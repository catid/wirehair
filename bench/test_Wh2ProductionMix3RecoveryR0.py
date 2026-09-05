#!/usr/bin/env python3
"""Synthetic recovery-protocol tests; never execute training or holdout roots."""
import copy
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2ProductionMix3RecoveryR0 as c


SOURCE_COMMIT = "a" * 40
LIBRARY_PATH = "/test/libwirehair.so.2.0.0"


def encoded(records):
    return "".join(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
                   for record in records)


def cell(phase, arm, K, attempt, root_index, schedule, outcomes):
    roots = c.TRAINING_ROOTS if phase == "train" else c.HOLDOUT_ROOTS
    root = roots[root_index]
    ids, attempted = c.trace(K, root, schedule)
    profile = c.profile_bytes(K, attempt)
    source_hash = c.digest(c.source_bytes(K))
    # Systematics are source bytes; arbitrary repair payloads remain stable for
    # each exact encoder/packet ID.  No native codec executes these fixtures.
    packets = b"".join(c.source_bytes(K)[packet_id * 2:packet_id * 2 + 2]
                       if packet_id < K else
                       bytes(((K + attempt + packet_id) & 255,
                              (3 * K + attempt + packet_id) & 255))
                       for packet_id in ids)
    return {
        "type": "cell", "phase": phase, "arm": arm, "K": K,
        "attempt": attempt, "root_index": root_index, "root": root,
        "schedule": schedule, "profile_hex": profile.hex(),
        "profile_sha256": c.digest(profile), "source_sha256": source_hash,
        "attempted_candidates": attempted,
        "trace_sha256": c.digest(b"".join(struct.pack("<I", i) for i in ids)),
        "ids": ids, "packets_hex": None if outcomes[0] == 10 else packets.hex(),
        "outcomes": list(outcomes),
        "recovered_sha256": [source_hash if outcome == 0 else None for outcome in outcomes],
    }


def transcript(phase="train", selected=(1, 2), baseline=(0, 0), failures=None):
    """Complete synthetic phases, with real-codec execution deliberately absent."""
    roots = c.TRAINING_ROOTS if phase == "train" else c.HOLDOUT_ROOTS
    records = [{"type": "begin", "protocol": c.PROTOCOL, "phase": phase,
                "source_commit": SOURCE_COMMIT, "library_path": LIBRARY_PATH}]
    for index, K in enumerate(c.KS):
        attempts = (list(range(selected[index] + 1)) if selected[index] >= 0
                    else list(range(256))) if phase == "train" else [selected[index]]
        arms = [("baseline", baseline[index])]
        arms.extend(("candidate", attempt) for attempt in attempts)
        for arm, attempt in arms:
            for root_index in range(len(roots)):
                for schedule in c.SCHEDULES:
                    outcomes = (0, 0, 0, 0)
                    if (phase == "train" and
                            (selected[index] == -1 or attempt < selected[index]) and
                            root_index == 0 and schedule == c.SCHEDULES[0]):
                        outcomes = (1, 0, 0, 0)
                    if failures is not None:
                        overridden = failures(phase, arm, K, attempt, root_index, schedule)
                        if overridden is not None:
                            outcomes = overridden
                    records.append(cell(phase, arm, K, attempt, root_index, schedule, outcomes))
    records.append({"type": "terminal", "phase": phase, "rows": len(records) - 1,
                    "selected_attempts": list(selected)})
    return records


def parse(records, phase="train", selected=None, **kwargs):
    return c.parse_phase(encoded(records), phase, source_commit=SOURCE_COMMIT,
                         library_path=LIBRARY_PATH, selected_attempts=selected, **kwargs)


def bundle_fixture(selected=(1, 2)):
    train_text = encoded(transcript(selected=selected))
    train = c.parse_phase(train_text, "train", SOURCE_COMMIT, LIBRARY_PATH)
    files = {"freeze.json": c.frozen_protocol({"source_commit": SOURCE_COMMIT,
                                               "library": LIBRARY_PATH}),
             "train.jsonl": train_text.encode(), "train.stderr": b"",
             "selection.json": {"protocol": c.PROTOCOL, "selected_attempts": list(selected),
                                "train_sha256": c.digest(train_text.encode())}}
    holdout = None
    if -1 not in selected:
        holdout_text = encoded(transcript("holdout", selected=selected))
        holdout = c.parse_phase(holdout_text, "holdout", SOURCE_COMMIT, LIBRARY_PATH, list(selected))
        files.update({"holdout.jsonl": holdout_text.encode(), "holdout.stderr": b""})
    summary = c.summarize(train, holdout)
    summary.update({"protocol": c.PROTOCOL, "elapsed_seconds": 0.25,
                    "promotion_claimed": False, "all_K_claimed": False, "timing_claimed": False})
    files["summary.json"] = summary
    return files


def publish_bundle(root, files):
    data = {name: value if type(value) is bytes else c.canonical(value) + b"\n"
            for name, value in files.items()}
    for name, value in data.items():
        (root / name).write_bytes(value)
    # The mutations deliberately recompute this manifest: scientific replay
    # must validate the contract and observations, not only file integrity.
    (root / "COMPLETE").write_bytes(c.canonical({name: c.digest(value)
                                                for name, value in data.items()}) + b"\n")


class ProductionMix3RecoveryTest(unittest.TestCase):
    def test_frozen_constants_roots_and_profile(self):
        self.assertEqual(c.KS, (3, 6))
        self.assertEqual(c.OVERHEADS, (0, 1, 2, 4))
        self.assertEqual(c.SCHEDULES, ("burst", "adversarial", "repair-only"))
        self.assertEqual(c.PROTOCOL, "wirehair.wh2.production-mix3-k3k6-r0")
        expected_training = ("0x7ccd510f122fc160", "0xb889883a79549774",
                             "0xb5666de0987896af") + tuple(
            "0x" + hashlib.sha256((c.PROTOCOL + ":train/" + str(i)).encode()).hexdigest()[:16]
            for i in range(3))
        expected_holdout = tuple(
            "0x" + hashlib.sha256((c.PROTOCOL + ":holdout/" + str(i)).encode()).hexdigest()[:16]
            for i in range(9))
        self.assertEqual(c.TRAINING_ROOTS, expected_training)
        self.assertEqual(c.HOLDOUT_ROOTS, expected_holdout)
        self.assertEqual(len(set(c.TRAINING_ROOTS + c.HOLDOUT_ROOTS)), 15)
        for K in c.KS:
            self.assertEqual(c.source_bytes(K), bytes((73 * i + 19 * K + 11) & 255
                                                     for i in range(2 * K)))
            for attempt in (0, 1, 255):
                expected = (b"WHV2" + struct.pack("<HHQQI", 1, 32,
                            0x4B295BBB47F4F9C9, K * 2, 2) + bytes((attempt, 0, 0, 0)))
                self.assertEqual(c.profile_bytes(K, attempt), expected)

    def test_complete_training_and_holdout(self):
        train = parse(transcript())
        self.assertEqual(list(train["selected_attempts"]), [1, 2])
        self.assertEqual(list(train["baseline_attempts"]), [0, 0])
        self.assertEqual(len(train["cells"]), 126)
        holdout = parse(transcript("holdout"), "holdout", [1, 2])
        self.assertEqual(len(holdout["cells"]), 108)
        self.assertEqual(c.summarize(train, holdout)["outcome"], "PASS")

    def test_unchanged_success_and_unchanged_holdout_failure_are_distinct(self):
        train = parse(transcript(selected=(0, 0)))
        holdout = parse(transcript("holdout", selected=(0, 0)), "holdout", [0, 0])
        self.assertEqual(c.summarize(train, holdout)["outcome"], "NO_CHANGE")

        def failures(phase, arm, K, attempt, root_index, schedule):
            if K == 3 and root_index == 0 and schedule == "repair-only":
                return (1, 1, 0, 0)

        failed = parse(transcript("holdout", selected=(0, 0), failures=failures),
                       "holdout", [0, 0])
        self.assertEqual(c.summarize(train, failed)["outcome"], "FAIL")

    def test_actual_baseline_attempt_is_not_assumed_zero(self):
        def bad_lower(phase, arm, K, attempt, root_index, schedule):
            if attempt == 0:
                return (10, 10, 10, 10)

        train = parse(transcript(selected=(1, 1), baseline=(1, 1), failures=bad_lower))
        holdout = parse(transcript("holdout", selected=(1, 1), baseline=(1, 1)),
                        "holdout", [1, 1])
        summary = c.summarize(train, holdout)
        self.assertEqual(summary["outcome"], "NO_CHANGE")
        self.assertFalse(summary["map_changed"])
        self.assertEqual(list(summary["baseline_attempts"]), [1, 1])

    def test_repairs_introductions_and_all_overhead_tails_are_paired(self):
        def failures(phase, arm, K, attempt, root_index, schedule):
            if K == 3 and root_index == 0 and schedule == "burst" and arm == "baseline":
                return (1, 1, 0, 0)
            if K == 6 and root_index == 1 and schedule == "adversarial":
                return (1, 0, 0, 0) if arm == "baseline" else (1, 1, 1, 0)

        train = parse(transcript())
        holdout = parse(transcript("holdout", failures=failures), "holdout", [1, 2])
        summary = c.summarize(train, holdout)
        self.assertEqual(summary["outcome"], "FAIL")
        self.assertEqual(summary["training"], {
            "baseline_failures": [2, 0, 0, 0], "candidate_failures": [0, 0, 0, 0],
            "repairs": [2, 0, 0, 0], "introductions": [0, 0, 0, 0], "cells_per_arm": 36,
        })
        self.assertEqual(summary["holdout"], {
            "baseline_failures": [2, 1, 0, 0], "candidate_failures": [1, 1, 1, 0],
            "repairs": [1, 1, 0, 0], "introductions": [0, 1, 1, 0], "cells_per_arm": 54,
        })

    def test_candidate_holdout_failure_does_not_reselect(self):
        train = parse(transcript())

        def failures(phase, arm, K, attempt, root_index, schedule):
            if arm == "candidate" and K == 6 and root_index == 8 and schedule == "adversarial":
                return (1, 0, 0, 0)

        holdout = parse(transcript("holdout", failures=failures), "holdout", [1, 2])
        self.assertEqual(c.summarize(train, holdout)["outcome"], "FAIL")
        self.assertEqual(list(holdout["selected_attempts"]), [1, 2])
        with self.assertRaises(c.ValidationError):
            parse(transcript("holdout"), "holdout", [2, 1])
        with self.assertRaises(c.ValidationError):
            parse(transcript("holdout"), "holdout")

    def test_cross_phase_maps_cannot_change(self):
        train = parse(transcript())
        altered_baseline = parse(transcript("holdout", baseline=(1, 0)),
                                 "holdout", [1, 2])
        with self.assertRaises(c.ValidationError):
            c.summarize(train, altered_baseline)
        altered_selection = parse(transcript("holdout", selected=(2, 1)),
                                  "holdout", [2, 1])
        with self.assertRaises(c.ValidationError):
            c.summarize(train, altered_selection)

    def test_bad_seed_is_whole_construction_failure_not_internal_error(self):
        def bad_seed(phase, arm, K, attempt, root_index, schedule):
            if arm == "candidate" and attempt == 1 and K == 6:
                return (10, 10, 10, 10)

        train = parse(transcript(failures=bad_seed))
        self.assertEqual(list(train["selected_attempts"]), [1, 2])
        rows = transcript(failures=bad_seed)
        row = next(row for row in rows if row.get("outcomes") == [10, 10, 10, 10])
        row["outcomes"][0] = 2
        with self.assertRaises(c.ValidationError):
            parse(rows)

    def test_full_attempt_exhaustion_is_not_success_or_missing_witness(self):
        exhausted = transcript(selected=(-1, 0))
        train = parse(exhausted)
        self.assertEqual(list(train["selected_attempts"]), [-1, 0])
        self.assertEqual(c.summarize(train, None)["outcome"], "EXHAUSTED")
        # A terminal -1 cannot excuse stopping the attempt search at 254.
        truncated = [row for row in exhausted if not (
            row.get("arm") == "candidate" and row.get("K") == 3 and row.get("attempt") == 255)]
        truncated[-1]["rows"] = len(truncated) - 2
        with self.assertRaises(c.ValidationError):
            parse(truncated)

    def test_lower_attempt_and_complete_coordinate_witnesses_required(self):
        for predicate in (
                lambda row: row.get("arm") == "candidate" and row.get("K") == 6 and row.get("attempt") == 1,
                lambda row: row.get("arm") == "candidate" and row.get("K") == 3 and row.get("root_index") == 5):
            records = [row for row in transcript() if not predicate(row)]
            records[-1]["rows"] = len(records) - 2
            with self.assertRaises(c.ValidationError):
                parse(records)

    def test_first_success_cannot_be_skipped_or_followed_by_more_attempts(self):
        records = transcript()
        for row in records:
            if row.get("arm") == "candidate" and row.get("K") == 6 and row.get("attempt") == 1:
                row["outcomes"] = [0, 0, 0, 0]
                row["recovered_sha256"] = [row["source_sha256"]] * 4
        with self.assertRaises(c.ValidationError):
            parse(records)

    def test_transcript_structure_and_provenance_are_exact(self):
        mutations = (
            lambda rows: rows.pop(),
            lambda rows: rows.pop(1),
            lambda rows: rows.insert(1, copy.deepcopy(rows[1])),
            lambda rows: rows.__setitem__(slice(1, 3), list(reversed(rows[1:3]))),
            lambda rows: rows[0].__setitem__("source_commit", "b" * 40),
            lambda rows: rows[0].__setitem__("library_path", "/other/libwirehair.so"),
            lambda rows: rows[0].__setitem__("protocol", "retired-protocol"),
            lambda rows: rows[1].__setitem__("extra", 0),
            lambda rows: rows[1].__setitem__("phase", "holdout"),
            lambda rows: rows[-1].__setitem__("rows", True),
            lambda rows: rows[-1].__setitem__("selected_attempts", [1, True]),
        )
        for mutation in mutations:
            rows = transcript()
            mutation(rows)
            with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                parse(rows)

    def test_cell_coordinates_profiles_traces_and_payload_shape_are_exact(self):
        mutations = (
            ("K", True), ("attempt", True), ("root_index", True),
            ("root", c.HOLDOUT_ROOTS[0]), ("schedule", "iid"),
            ("profile_hex", "00" * 32), ("profile_sha256", "b" * 64),
            ("source_sha256", "b" * 64), ("trace_sha256", "b" * 64),
            ("attempted_candidates", 0), ("ids", [0] * 7),
            ("packets_hex", "0000"), ("packets_hex", None),
            ("outcomes", [1, 0, 1, 0]), ("outcomes", [False, 0, 0, 0]),
            ("outcomes", [0, 0, 0]), ("outcomes", [1, 10, 10, 10]),
            ("recovered_sha256", ["b" * 64] * 4),
        )
        for key, value in mutations:
            rows = transcript()
            rows[1][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(c.ValidationError):
                parse(rows)
        rows = transcript()
        rows[1]["trace_sha256"] = c.digest(b"".join(struct.pack(">I", i) for i in rows[1]["ids"]))
        with self.assertRaises(c.ValidationError):
            parse(rows)

    def test_same_attempt_semantics_must_agree_between_arms(self):
        rows = transcript(selected=(0, 0))
        row = next(row for row in rows if row.get("arm") == "candidate" and
                   row.get("schedule") == "repair-only")
        row["packets_hex"] = "ffff" + row["packets_hex"][4:]
        with self.assertRaises(c.ValidationError):
            parse(rows)

    def test_bad_seed_cannot_vary_by_loss_coordinate(self):
        rows = transcript()
        row = next(row for row in rows if row.get("arm") == "candidate" and
                   row.get("K") == 6 and row.get("attempt") == 1)
        row["outcomes"] = [10, 10, 10, 10]
        row["packets_hex"] = None
        row["recovered_sha256"] = [None] * 4
        with self.assertRaises(c.ValidationError):
            parse(rows)

    def test_malformed_json_is_rejected(self):
        for text in ("", "{}", '{"type":"begin","type":"cell"}\n',
                     '{"x":NaN}\n', encoded(transcript())[:-1]):
            with self.subTest(text=text[:60]), self.assertRaises(c.ValidationError):
                c.parse_phase(text, "train")

    def test_offline_replay_never_calls_live_tools(self):
        for selected, outcome in (((1, 2), "PASS"), ((0, 0), "NO_CHANGE"), ((-1, 0), "EXHAUSTED")):
            with self.subTest(outcome=outcome), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                publish_bundle(root, bundle_fixture(selected))
                with mock.patch.object(c, "command", side_effect=AssertionError("live command in replay")):
                    self.assertEqual(c.replay(root)["outcome"], outcome)

    def test_replay_rejects_rehashed_contract_mutations(self):
        for key, value in (("K", [3]), ("block_bytes", 64), ("schedules", ["iid"]),
                           ("loss_ppm", 0), ("overheads", [0]), ("worker_budget_seconds", 61),
                           ("outer_budget_seconds", 71), ("max_prefix_decodes", 37441),
                           ("training_roots", list(c.HOLDOUT_ROOTS)), ("unexpected", True)):
            files = bundle_fixture()
            files["freeze.json"][key] = value
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                publish_bundle(root, files)
                with self.assertRaises(c.ValidationError):
                    c.replay(root)

    def test_replay_rejects_rehashed_errors_rosters_maps_and_summaries(self):
        mutations = (
            lambda files: files.__setitem__("train.stderr", b"worker error\n"),
            lambda files: files.__setitem__("holdout.stderr", b"worker error\n"),
            lambda files: files.pop("holdout.jsonl"),
            lambda files: files.__setitem__("unexpected.json", {}),
            lambda files: files["selection.json"].__setitem__("selected_attempts", [2, 1]),
            lambda files: files["selection.json"].__setitem__("train_sha256", "b" * 64),
            lambda files: files["summary.json"].__setitem__("elapsed_seconds", 70),
            lambda files: files["summary.json"].__setitem__("elapsed_seconds", -1),
            lambda files: files["summary.json"].__setitem__("elapsed_seconds", True),
            lambda files: files["summary.json"].__setitem__("promotion_claimed", True),
            lambda files: files["summary.json"].__setitem__("outcome", "NO_CHANGE"),
            lambda files: files["summary.json"].__setitem__("universal_success", True),
            lambda files: files["summary.json"]["training"]["repairs"].__setitem__(0, 99),
        )
        for mutation in mutations:
            files = bundle_fixture()
            mutation(files)
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                publish_bundle(root, files)
                with self.assertRaises(c.ValidationError):
                    c.replay(root)

    def test_publication_hashing_cannot_seal_success_after_deadline(self):
        for expire_during_hashing in (False, True):
            files = bundle_fixture()
            summary = files.pop("summary.json")
            clock = {"now": 69.0}
            original_digest = c.file_digest
            original_write = c.write_json

            def hash_file(path):
                if expire_during_hashing:
                    clock["now"] = 71.0
                return original_digest(path)

            def write_file(path, value):
                if Path(path).name == "COMPLETE" and expire_during_hashing:
                    # Inspect the bytes being sealed, not only the return value.
                    published = json.loads((Path(path).parent / "summary.json").read_bytes())
                    self.assertEqual(published["outcome"], "INVALID")
                original_write(path, value)

            with self.subTest(expired=expire_during_hashing), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                for name, value in files.items():
                    (root / name).write_bytes(value if type(value) is bytes else c.canonical(value) + b"\n")
                with mock.patch.object(c.time, "monotonic", side_effect=lambda: clock["now"]), \
                        mock.patch.object(c, "file_digest", side_effect=hash_file), \
                        mock.patch.object(c, "write_json", side_effect=write_file):
                    published = c.publish_result(root, summary, started=0.0)
                expected = "INVALID" if expire_during_hashing else "PASS"
                self.assertEqual(published["outcome"], expected)
                self.assertEqual(c.replay(root)["outcome"], expected)
                complete = json.loads((root / "COMPLETE").read_bytes())
                for name, expected_hash in complete.items():
                    self.assertEqual(original_digest(root / name), expected_hash)


class UntimedWorkerSelftest(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("WH2_PRODUCTION_MIX3_TEST_WORKER"),
                         "set WH2_PRODUCTION_MIX3_TEST_WORKER for root-free worker selftest")
    def test_root_free_public_api_selftest_only(self):
        worker = Path(os.environ["WH2_PRODUCTION_MIX3_TEST_WORKER"]).resolve(strict=True)
        result = subprocess.run([str(worker), "--selftest"], capture_output=True,
                                text=True, timeout=15, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        records = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(records[-1], {"type": "selftest", "pass": True})
        self.assertEqual(len(records), 13)
        for row in records[:-1]:
            self.assertEqual(row["phase"], "selftest")
            self.assertEqual(row["root"], "0x1234567890abcdef")
            self.assertNotIn(row["root"], c.TRAINING_ROOTS + c.HOLDOUT_ROOTS)
        for index, K in enumerate(c.KS):
            offset = index * 6
            for baseline, candidate in zip(records[offset:offset + 3],
                                            records[offset + 3:offset + 6]):
                self.assertEqual(baseline["K"], K)
                baseline = {key: value for key, value in baseline.items() if key != "arm"}
                candidate = {key: value for key, value in candidate.items() if key != "arm"}
                self.assertEqual(baseline, candidate)


if __name__ == "__main__":
    unittest.main()
