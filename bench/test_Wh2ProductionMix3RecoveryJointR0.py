#!/usr/bin/env python3
"""Joint-width protocol tests with synthetic outcomes, never codec campaigns."""
import copy
import contextlib
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2ProductionMix3RecoveryPortabilityR0 as previous

PREVIOUS_STATE = (previous.c.PROTOCOL, previous.c.HOLDOUT_ROOTS, previous.c.trace,
                  previous.c.source_bytes, previous.c.profile_bytes, previous.c.validate_cell)
import Wh2ProductionMix3RecoveryJointR0 as j

c = j.c
HISTORY_BUNDLES = tuple(Path("/var/tmp/" + name) for name in (
    "wh2-production-mix3-k3k6-r0", "wh2-production-mix3-k3k6-broad-r0",
    "wh2-production-mix3-k3k6-portability-r0"))


def encoded(rows):
    return b"".join(c.canonical(row) + b"\n" for row in rows)


def synthetic_receipt():
    sources = j.source_receipt()
    basis = c.parse_json(j.b.BASIS_FREEZE.read_bytes())["build"]
    return dict(protocol=j.PROTOCOL, source_commit="a" * 40, source_files=sources,
                worker=str(j.WORKER_SOURCE), worker_sha256=sources[str(j.WORKER_SOURCE.relative_to(j.ROOT))],
                interpreter=str(Path(sys.executable).resolve()), interpreter_sha256="b" * 64,
                interpreter_version="synthetic", interpreter_flags=["-I", "-B", "-S"],
                library=str(j.b.PINNED_LIBRARY), library_sha256=j.b.LIBRARY_SHA256,
                production_source_commit=j.b.PRODUCTION_COMMIT, production_build=basis)


class FakeNative:
    """No CDLL: deterministic packet bytes and controlled constructor/decoder outcomes."""
    def __init__(self, history, selected=(5, 6)):
        self.history, self.selected, self.opened, self.closed, self.calls = history, selected, [], [], []
        self.regression_targets, self.fresh_targets = {}, {}
        for K in j.KS:
            specs = j.specs_for(K, "search", history)
            self.regression_targets[K] = specs[2]
            regression = {(x["B"], tuple(x["ids"])) for x in specs if x["kind"] == "regression"}
            self.fresh_targets[K] = next(x for x in specs if x["kind"] == "fresh" and
                                         (x["B"], tuple(x["ids"][:K])) not in regression)

    @contextlib.contextmanager
    def encoder(self, K, B, attempt):
        actual = 9 if attempt is None else attempt
        choice = self.selected[j.KS.index(K)]
        bad = attempt is not None and (choice == -1 or [K, B, attempt] in j.BAD_PROFILES or
                                      (K == 6 and attempt in (0, 2, 4)))
        source = j.p.source_bytes(K, B)
        handle = (K, B, actual)
        self.opened.append(handle)
        try:
            yield (handle, j.p.profile_bytes(K, B, actual), j.p.ct.create_string_buffer(source, len(source)),
                   source, bad, actual)
        finally:
            self.closed.append(handle)

    def encode(self, handle, packet_id, B):
        K, width, attempt = handle
        assert width == B
        if packet_id < K:
            return j.p.source_bytes(K, B)[packet_id * B:(packet_id + 1) * B]
        return bytes((K + attempt + packet_id + i) & 255 for i in range(B))

    def decode(self, profile, source, B, ids, packets, count):
        K, attempt = len(source) // B, profile[28]
        self.calls.append((K, B, attempt, count))
        choice = self.selected[j.KS.index(K)]
        if attempt < choice:
            # Exercise both regression and fresh-training failure witnesses.
            failed = (self.fresh_targets if attempt == 5 else self.regression_targets)[K]
            if B == failed["B"] and ids[:count] == failed["ids"][:count]:
                return 1
        return 0


def fake_stream(phase, freeze, selection=None, selected=(5, 6)):
    native, rows = FakeNative(freeze["history"], selected), []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        c.write_json(root / "freeze.json", freeze)
        selection_path = root / "selection.json" if selection is not None else None
        if selection is not None:
            c.write_json(selection_path, selection)
        with mock.patch.object(j.p, "Native", return_value=native), \
                mock.patch.object(j, "check_build"), mock.patch.object(j.p, "emit", side_effect=rows.append):
            j.worker(phase, root / "freeze.json", 60, selection_path,
                     c.file_digest(selection_path) if selection_path else None)
    assert sorted(native.opened) == sorted(native.closed)
    return rows, native


def write_bundle(root, files):
    data = {name: value if isinstance(value, bytes) else c.canonical(value) + b"\n"
            for name, value in files.items()}
    for name, value in data.items():
        (root / name).write_bytes(value)
    (root / "COMPLETE").write_bytes(c.canonical({name: c.digest(value) for name, value in data.items()}) + b"\n")


def scoped(result):
    phases = ["search"] if result["outcome"] == "EXHAUSTED" else ["search", "validate", "holdout"]
    if result["outcome"] == "INVALID":
        phases = list(result["completed_phases"])
        if len(phases) < 3:
            phases.append(("search", "validate", "holdout")[len(phases)])
    return dict(result, protocol=j.PROTOCOL, elapsed_seconds=0.25,
                workers=[dict(phase=phase, elapsed_seconds=0.01) for phase in phases],
                promotion_claimed=False, all_K_claimed=False, timing_claimed=False)


class JointProtocolTest(unittest.TestCase):
    def test_private_import_and_frozen_root_roster(self):
        self.assertIsNot(j.p, previous)
        self.assertIsNot(c, previous.c)
        self.assertEqual((previous.c.PROTOCOL, previous.c.HOLDOUT_ROOTS, previous.c.trace,
                          previous.c.source_bytes, previous.c.profile_bytes, previous.c.validate_cell),
                         PREVIOUS_STATE)
        self.assertEqual(j.PROTOCOL, "wirehair.wh2.production-mix3-k3k6-joint-r0")
        roots = {phase: tuple("0x" + hashlib.sha256((j.PROTOCOL + ":" + phase + "/" + str(i))
                                                  .encode("ascii")).hexdigest()[:16]
                              for i in range(count)) for phase, count in (("train", 16), ("holdout", 64))}
        self.assertEqual(j.TRAINING_ROOTS, roots["train"])
        self.assertEqual(j.HOLDOUT_ROOTS, roots["holdout"])
        self.assertEqual(c.digest(c.canonical(roots)),
                         "8245376b0c91593e099f35207b436a0ed6343fd212bc6119f6357d74c1fd2470")
        self.assertEqual(len(set(roots["train"] + roots["holdout"])), 80)
        self.assertEqual((j.STREAM_LIMIT, j.TOTAL_LIMIT), (33554432, 67108864))

    @unittest.skipUnless(all(path.exists() for path in HISTORY_BUNDLES), "historical input bundles unavailable")
    def test_existing_history_ledgers_are_exact_read_only_inputs(self):
        before = {str(path / "COMPLETE"): c.file_digest(path / "COMPLETE") for path in HISTORY_BUNDLES}
        with mock.patch.object(c, "command", side_effect=AssertionError("live history tool")), \
                mock.patch.object(j.p, "Native", side_effect=AssertionError("native history call")):
            history = j.load_history()
            j.validate_history(history)
        self.assertEqual((len(history["prefixes"]), len(history["origins"]), len(history["bad_origins"])),
                         (73, 126, 360))
        self.assertEqual(c.digest(c.canonical(history["prefixes"])),
                         "e12d079b994aa04d7dd47656cf42895c1910f4239434ee8f3a2aae9bb29e0b1b")
        self.assertEqual(c.digest(c.canonical(history["origins"])),
                         "5362613e8ade02cec8932fad0fc1fd0e290eaa343b2eb11471ebd79a094034ff")
        self.assertEqual(c.digest(c.canonical(history["bad_origins"])),
                         "617c0db46b2a719d14f780f3ed44872c5ecf65d6f499e2500d3438d89dcf4613")
        self.assertEqual(len(set(history["prior_roots"])), 127)
        self.assertFalse(set(history["prior_roots"]) & set(j.TRAINING_ROOTS + j.HOLDOUT_ROOTS))
        self.assertEqual(before, {path: c.file_digest(path) for path in before})
        for key in ("prefixes", "origins", "bad_origins"):
            changed = copy.deepcopy(history)
            changed[key] = changed[key][:-1]
            with self.subTest(ledger=key), self.assertRaises(c.ValidationError):
                j.validate_history(changed)

    def test_supervisor_enforces_stream_and_cross_phase_aggregate_caps(self):
        def child(stdout, stderr=0):
            return [sys.executable, "-I", "-B", "-S", "-c",
                    "import os; os.write(1,b'x'*%d); os.write(2,b'e'*%d)" % (stdout, stderr)]

        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(j, "STREAM_LIMIT", 16), mock.patch.object(j, "TOTAL_LIMIT", 32):
            root = Path(directory)
            prior = 0
            # Three separate phases reach the inclusive aggregate boundary.
            for phase, size in (("search", 15), ("validate", 15), ("holdout", 2)):
                prior = j.execute_worker(child(size), root / (phase + ".jsonl"),
                                         root / (phase + ".stderr"), 5, prior)
            self.assertEqual(prior, 32)
            for name, size, stderr, prior in (("stream", 16, 0, 0), ("aggregate", 1, 0, 32),
                                              ("stderr", 0, 2, 31)):
                with self.subTest(cap=name), self.assertRaisesRegex(c.ValidationError, "raw byte cap"):
                    j.execute_worker(child(size, stderr), root / name, root / (name + ".err"), 5, prior)
                actual = sum((root / file).stat().st_size for file in (name, name + ".err"))
                self.assertLess((root / name).stat().st_size, 16)
                self.assertLess((root / (name + ".err")).stat().st_size, 16)
                self.assertLessEqual(actual + prior, 32)
            with self.assertRaisesRegex(c.ValidationError, "worker failed/stderr"):
                j.execute_worker(child(1, 1), root / "diagnostic", root / "diagnostic.err", 5)

    def test_supervisor_timeout_kills_and_reaps_synthetic_child(self):
        children, real_popen = [], j.subprocess.Popen

        def spawn(*args, **kwargs):
            process = real_popen(*args, **kwargs)
            children.append(process)
            return process

        with tempfile.TemporaryDirectory() as directory, mock.patch.object(j.subprocess, "Popen", side_effect=spawn):
            root = Path(directory)
            with self.assertRaisesRegex(c.ValidationError, "worker deadline"):
                j.execute_worker([sys.executable, "-I", "-B", "-S", "-c", "import time; time.sleep(30)"],
                                 root / "raw", root / "error", 0.05)
        self.assertEqual(len(children), 1)
        self.assertLess(children[0].returncode, 0)
        with self.assertRaises(ChildProcessError):
            os.waitpid(children[0].pid, os.WNOHANG)


@unittest.skipUnless(all(path.exists() for path in HISTORY_BUNDLES), "historical fixture inputs unavailable")
class JointStreamTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.history = j.load_history()
        cls.freeze = j.frozen_protocol(synthetic_receipt(), cls.history)
        cls.search_rows, cls.fake = fake_stream("search", cls.freeze)
        cls.search_text = encoded(cls.search_rows)
        cls.search = j.parse_search(cls.search_text.decode("ascii"), cls.freeze)
        cls.selection = dict(protocol=j.PROTOCOL, selected_attempts=[5, 6], search_sha256=c.digest(cls.search_text))
        cls.validation_rows, _ = fake_stream("validate", cls.freeze, cls.selection)
        cls.holdout_rows, _ = fake_stream("holdout", cls.freeze, cls.selection)
        cls.validation = cls.parse(cls.validation_rows, "validate")
        cls.holdout = cls.parse(cls.holdout_rows, "holdout")
        cls.files = {"freeze.json": cls.freeze, "selection.json": cls.selection,
                     "search.jsonl": cls.search_text, "search.stderr": b"",
                     "validate.jsonl": encoded(cls.validation_rows), "validate.stderr": b"",
                     "holdout.jsonl": encoded(cls.holdout_rows), "holdout.stderr": b"",
                     "summary.json": scoped(j.summarize(cls.search, cls.validation, cls.holdout))}

    @classmethod
    def parse(cls, rows, phase):
        return j.parse_validation(encoded(rows).decode("ascii"), phase, cls.freeze, cls.selection, cls.search)

    def test_complete_fake_search_and_every_validation_constraint(self):
        self.assertEqual(self.search["selected_attempts"], [5, 6])
        self.assertEqual(len(self.search["records"]), 13)
        self.assertEqual(len(self.validation["rows"]), 1014)
        self.assertEqual(len(self.holdout["rows"]), 2304)
        self.assertEqual(self.files["summary.json"]["outcome"], "PASS")
        self.assertEqual(self.files["summary.json"]["total_prefix_cases"], self.search["prefix_cases"] + 11958)
        self.assertEqual(self.validation["baseline_attempts"], {"3": [9, 9, 9], "6": [9, 9, 9]})
        for K, count in ((3, 210), (6, 297)):
            specs = j.specs_for(K, "search", self.history)
            self.assertEqual(len(specs), count)
            self.assertEqual(sum(x["kind"] == "fresh" for x in specs), 144)
            kinds = [x["kind"] for x in specs]
            self.assertEqual(kinds, sorted(kinds, key=lambda kind: kind == "fresh"))
            accepted = next(row for row in self.search["records"] if row["K"] == K and row["accepted"])
            self.assertEqual(accepted["passed"], count)
        self.assertTrue(any(row["witness"] and row["witness"]["kind"] == "fresh" for row in self.search["records"]))
        self.assertTrue(any(row["witness"] and row["passed"] > 0 for row in self.search["records"]))

    def test_search_rejects_skipped_attempts_wrong_first_witness_and_late_success(self):
        for mutation in ("skip", "late", "constructor_order", "constructor_after_bad", "passed",
                         "witness", "success_hash", "map"):
            rows = copy.deepcopy(self.search_rows)
            if mutation == "skip":
                del rows[1]
            elif mutation == "late":
                rows.insert(7, copy.deepcopy(rows[6]))
            elif mutation == "constructor_order":
                rows[1]["constructors"].reverse()
            elif mutation == "constructor_after_bad":
                target = next(row for row in rows[1:-1] if row["constructors"][-1]["outcome"] == 10)
                target["constructors"].append(copy.deepcopy(target["constructors"][-1]))
            elif mutation in ("passed", "witness"):
                target = next(row for row in rows[1:-1] if row["witness"])
                if mutation == "passed":
                    target["passed"] += 1
                else:
                    target["witness"]["outcomes"] = [0]
            elif mutation == "success_hash":
                target = next(row for row in rows[1:-1] if row["passed"] == 0)
                target["pass_sha256"] = "f" * 64
            else:
                rows[-1]["selected_attempts"] = [6, 6]
            rows[-1]["rows"] = len(rows) - 2
            with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                j.parse_search(encoded(rows).decode("ascii"), self.freeze)
        rows = copy.deepcopy(self.search_rows)
        known = next(row for row in rows[1:-1] if row["K"] == 3 and row["attempt"] == 1)
        known["constructors"][0]["outcome"] = 0
        with self.assertRaisesRegex(c.ValidationError, "known BadSeed"):
            j.parse_search(encoded(rows).decode("ascii"), self.freeze)

    def test_selected_byte_hash_revalidation_and_map_cannot_drift(self):
        for mutation in ("selected_failure", "repair_byte", "systematic_tail", "map", "search_hash", "selection_hash"):
            rows = copy.deepcopy(self.validation_rows)
            if mutation == "selected_failure":
                target = next(row for row in rows if row.get("arm") == "candidate")
                target["outcomes"][0], target["recovered_sha256"][0] = 1, None
            elif mutation in ("repair_byte", "systematic_tail"):
                systematic = mutation == "systematic_tail"
                target = next(row for row in rows if row.get("arm") == "candidate" and row["B"] == 1280 and
                              any((packet_id < row["K"]) == systematic for packet_id in row["ids"]))
                index = next(i for i, packet_id in enumerate(target["ids"]) if (packet_id < target["K"]) == systematic)
                payload = bytearray.fromhex(target["packets_hex"])
                payload[(index + 1) * 1280 - 1] ^= 1
                target["packets_hex"] = payload.hex()
            elif mutation == "map":
                rows[-1]["selected_attempts"] = [5, 7]
            elif mutation == "selection_hash":
                rows[0]["selection_sha256"] = "f" * 64
            else:
                search = copy.deepcopy(self.search)
                next(row for row in search["records"] if row["accepted"])["pass_sha256"] = "f" * 64
                with self.assertRaises(c.ValidationError):
                    j.parse_validation(encoded(rows).decode(), "validate", self.freeze, self.selection, search)
                continue
            with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                self.parse(rows, "validate")
        for path in j.IMMUTABLE_HELPERS:
            receipt = copy.deepcopy(self.freeze["build"])
            receipt["source_files"][path] = "f" * 64
            with self.subTest(helper=path), self.assertRaisesRegex(c.ValidationError, "immutable helper"):
                j.frozen_protocol(receipt, self.history)

    def test_offline_replay_pass_fail_and_rehashed_contract_mutations(self):
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(c, "command", side_effect=AssertionError("live replay tool")), \
                mock.patch.object(j.p, "Native", side_effect=AssertionError("native replay")), \
                mock.patch.object(j, "load_history", side_effect=AssertionError("live historical replay")):
            root = Path(directory)
            write_bundle(root, self.files)
            self.assertEqual(j.replay(root)["outcome"], "PASS")
            failed = copy.deepcopy(self.holdout_rows)
            target = next(row for row in failed if row.get("arm") == "candidate" and row["K"] == 6 and row["B"] == 64)
            target["outcomes"][0], target["recovered_sha256"][0] = 1, None
            files = dict(self.files, **{"holdout.jsonl": encoded(failed),
                                      "summary.json": scoped(j.summarize(self.search, self.validation, self.parse(failed, "holdout")))})
            write_bundle(root, files)
            self.assertEqual(j.replay(root)["outcome"], "FAIL")
            for mutation in ("prior_roots", "widths", "map", "stderr", "summary", "late_time",
                             "worker_order", "worker_total", "worker_type"):
                files = copy.deepcopy(self.files)
                if mutation == "prior_roots":
                    files["freeze.json"]["history"]["prior_roots"][0] = "0x0000000000000000"
                elif mutation == "widths":
                    files["freeze.json"]["widths"].reverse()
                elif mutation == "map":
                    files["selection.json"]["selected_attempts"] = [5, 7]
                elif mutation == "stderr":
                    files["validate.stderr"] = b"diagnostic\n"
                elif mutation == "summary":
                    files["summary.json"]["total_prefix_cases"] += 1
                elif mutation == "late_time":
                    files["summary.json"]["elapsed_seconds"] = 70
                elif mutation == "worker_order":
                    files["summary.json"]["workers"].reverse()
                elif mutation == "worker_total":
                    files["summary.json"]["elapsed_seconds"] = 60
                    for worker in files["summary.json"]["workers"]:
                        worker["elapsed_seconds"] = 20
                else:
                    files["summary.json"]["workers"][0]["elapsed_seconds"] = True
                write_bundle(root, files)
                with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                    j.replay(root)

    def test_exhaustion_has_full_lower_witnesses_and_no_following_phase(self):
        rows, native = fake_stream("search", self.freeze, selected=(-1, -1))
        text = encoded(rows)
        search = j.parse_search(text.decode("ascii"), self.freeze)
        self.assertEqual((len(search["records"]), search["selected_attempts"]), (512, [-1, -1]))
        self.assertEqual(native.calls, [])
        selection = dict(protocol=j.PROTOCOL, selected_attempts=[-1, -1], search_sha256=c.digest(text))
        files = {"freeze.json": self.freeze, "search.jsonl": text, "search.stderr": b"",
                 "selection.json": selection, "summary.json": scoped(j.summarize(search))}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_bundle(root, files)
            self.assertEqual(j.replay(root)["outcome"], "EXHAUSTED")
            files["validate.jsonl"] = b"partial\n"
            write_bundle(root, files)
            with self.assertRaises(c.ValidationError):
                j.replay(root)

    def test_controller_seals_selection_and_carries_bytes_across_phases(self):
        receipt = self.freeze["build"]

        def command(argv, *unused, **kwargs):
            return "" if argv[:2] in (["git", "status"], ["git", "ls-files"]) else receipt["source_commit"]

        real_search, real_validation = j.parse_search, j.parse_validation
        for invalid_validation, parse_cost in ((False, 0), (True, 0), (False, 11), (False, 31)):
            streams = {phase: self.files[phase + ".jsonl"] for phase in ("search", "validate", "holdout")}
            if invalid_validation:
                rows = copy.deepcopy(self.validation_rows)
                target = next(row for row in rows if row.get("arm") == "candidate")
                target["outcomes"][0], target["recovered_sha256"][0] = 1, None
                streams["validate"] = encoded(rows)
            seen, cumulative, now, timeouts = [], [0], [0.0], []

            def parse_search(*args):
                result = real_search(*args)
                now[0] += parse_cost
                return result

            def parse_validation(*args):
                result = real_validation(*args)
                now[0] += parse_cost
                return result

            def transport(argv, raw, error, timeout, prior_bytes=0):
                phase = argv[argv.index("--worker") + 1]
                self.assertEqual(phase, ("search", "validate", "holdout")[len(seen)])
                self.assertEqual(prior_bytes, cumulative[0])
                self.assertTrue(0 < timeout <= 60)
                timeouts.append(timeout)
                now[0] += 10 if parse_cost else 0
                if phase == "search":
                    self.assertNotIn("--selection", argv)
                else:
                    selection_path = Path(argv[argv.index("--selection") + 1])
                    self.assertEqual(c.parse_json(selection_path.read_bytes()), self.selection)
                    self.assertEqual(c.file_digest(selection_path), argv[argv.index("--selection-sha256") + 1])
                self.assertEqual(argv[1:4], ["-I", "-B", "-S"])
                Path(raw).write_bytes(streams[phase])
                Path(error).write_bytes(b"")
                seen.append(phase)
                cumulative[0] += len(streams[phase])
                return cumulative[0]

            with self.subTest(invalid_validation=invalid_validation, parse_cost=parse_cost), tempfile.TemporaryDirectory() as directory:
                root, bundle = Path(directory), Path(directory) / "bundle"
                c.write_json(root / "build.json", receipt)
                c.write_json(root / "inputs.json", self.history)
                with mock.patch.object(j, "NAMESPACE", bundle), mock.patch.object(j, "check_build"), \
                        mock.patch.object(j, "load_history", return_value=self.history), \
                        mock.patch.object(c, "command", side_effect=command), \
                        mock.patch.object(j, "execute_worker", side_effect=transport), \
                        mock.patch.object(j.time, "monotonic", side_effect=lambda: now[0]), \
                        mock.patch.object(j, "parse_search", side_effect=parse_search), \
                        mock.patch.object(j, "parse_validation", side_effect=parse_validation), \
                        mock.patch.object(j.p, "Native", side_effect=AssertionError("native controller call")):
                    result = j.run_once(root / "build.json", bundle)
                self.assertEqual(j.replay(bundle), result)
                invalid = invalid_validation or parse_cost == 31
                self.assertEqual(result["outcome"], "INVALID" if invalid else "PASS")
                self.assertEqual(seen, ["search", "validate"] if invalid else ["search", "validate", "holdout"])
                if invalid:
                    self.assertEqual(result["completed_phases"], ["search"] if invalid_validation else ["search", "validate"])
                    self.assertFalse((bundle / "holdout.jsonl").exists())
                if parse_cost:
                    self.assertEqual(timeouts, [60, 49, 28] if parse_cost == 11 else [60, 29])
                    self.assertEqual([worker["elapsed_seconds"] for worker in result["workers"]], [10] * len(seen))

    def test_invalid_partial_phase_and_publication_deadline_keep_progress(self):
        files = {name: value for name, value in self.files.items() if not name.startswith(("holdout.", "validate."))}
        files["validate.jsonl"] = b"truncated current phase"
        files["validate.stderr"] = b"interrupted\n"
        files["summary.json"] = scoped(dict(outcome="INVALID", error="worker interrupted", completed_phases=["search"]))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_bundle(root, files)
            self.assertEqual(j.replay(root)["completed_phases"], ["search"])
            files["holdout.jsonl"] = b"forbidden later phase\n"
            write_bundle(root, files)
            with self.assertRaises(c.ValidationError):
                j.replay(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, value in self.files.items():
                if name != "summary.json":
                    (root / name).write_bytes(value if isinstance(value, bytes) else c.canonical(value) + b"\n")
            now, real_digest = [0.0], c.file_digest

            def delayed_digest(path):
                now[0] = 70.0
                return real_digest(path)

            with mock.patch.object(j.time, "monotonic", side_effect=lambda: now[0]), \
                    mock.patch.object(c, "file_digest", side_effect=delayed_digest):
                result = j.publish_result(root, self.files["summary.json"], 0.0)
            self.assertEqual(result["outcome"], "INVALID")
            self.assertEqual(result["completed_phases"], ["search", "validate", "holdout"])
            self.assertEqual(j.replay(root), result)


if __name__ == "__main__":
    unittest.main()
