#!/usr/bin/env python3
"""Width-local receipt tests: synthetic outcomes only, never codec campaigns."""
import contextlib
import copy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2ProductionMix3RecoveryJointR0 as previous

PREVIOUS_STATE = (previous.PROTOCOL, previous.WIDTHS, previous.TRAINING_ROOTS,
                  previous.HOLDOUT_ROOTS, previous.c.PROTOCOL)
import Wh2ProductionMix3RecoveryWidthLocalR0 as w

p, c = w.p, w.c
CELLS = ((3, 2), (3, 64), (3, 1280), (6, 2), (6, 64), (6, 1280))
HISTORY_BUNDLES = tuple(Path("/var/tmp/" + name) for name in (
    "wh2-production-mix3-k3k6-r0", "wh2-production-mix3-k3k6-broad-r0",
    "wh2-production-mix3-k3k6-portability-r0", "wh2-production-mix3-k3k6-joint-r0"))


def history_available():
    """These are machine-local immutable fixtures, not portable test assets."""
    for path in HISTORY_BUNDLES:
        if not (path / "COMPLETE").is_file():
            return False
        manifest = c.parse_json((path / "COMPLETE").read_bytes())
        if not all((path / name).is_file() for name in manifest):
            return False
    return w.j.b.BASIS_FREEZE.is_file()


def encoded(rows):
    return b"".join(c.canonical(row) + b"\n" for row in rows)


def write_bundle(root, files):
    data = {name: value if isinstance(value, bytes) else c.canonical(value) + b"\n"
            for name, value in files.items()}
    for name, value in data.items():
        (root / name).write_bytes(value)
    (root / "COMPLETE").write_bytes(c.canonical({name: c.digest(value) for name, value in data.items()}) + b"\n")


def synthetic_receipt():
    sources = w.j.source_receipt()
    basis = c.parse_json(w.j.b.BASIS_FREEZE.read_bytes())["build"]
    return dict(protocol=w.PROTOCOL, source_commit="a" * 40, source_files=sources,
                worker=str(w.j.WORKER_SOURCE),
                worker_sha256=sources[str(w.j.WORKER_SOURCE.relative_to(w.j.ROOT))],
                interpreter=str(Path(sys.executable).resolve()), interpreter_sha256="b" * 64,
                interpreter_version="synthetic", interpreter_flags=["-I", "-B", "-S"],
                library=str(w.j.b.PINNED_LIBRARY), library_sha256=w.j.b.LIBRARY_SHA256,
                production_source_commit=w.j.b.PRODUCTION_COMMIT, production_build=basis)


def scoped(result):
    phases = ["search"] if result["outcome"] == "EXHAUSTED" else ["search", "validate", "holdout"]
    return dict(result, protocol=w.PROTOCOL, elapsed_seconds=0.25,
                workers=[dict(phase=phase, elapsed_seconds=0.01) for phase in phases],
                promotion_claimed=False, all_K_claimed=False, timing_claimed=False)


class FakeNative:
    """Exact-width handles and byte oracles, with synthetic rank outcomes."""
    def __init__(self, history, selected=None):
        self.bad = {(row["K"], row["B"], row["attempt"]) for row in history["bad_origins"]}
        self.selected = tuple(selected) if selected is not None else tuple(
            next(a for a in range(6 + index, 256) if (K, B, a) not in self.bad)
            for index, (K, B) in enumerate(CELLS))
        self.baselines = tuple(next(a for a in range(100 + index, 256) if (K, B, a) not in self.bad)
                               for index, (K, B) in enumerate(CELLS))
        self.opened, self.closed, self.calls = [], [], []

    @contextlib.contextmanager
    def encoder(self, K, B, attempt):
        index = CELLS.index((K, B))
        actual = self.baselines[index] if attempt is None else attempt
        source = p.source_bytes(K, B)
        handle = (K, B, actual)
        bad = attempt is not None and (self.selected[index] == -1 or (K, B, actual) in self.bad)
        self.opened.append((K, B, attempt))
        try:
            yield (handle, p.profile_bytes(K, B, actual), p.ct.create_string_buffer(source, len(source)),
                   source, bad, actual)
        finally:
            self.closed.append((K, B, attempt))

    def encode(self, handle, packet_id, B):
        K, width, attempt = handle
        assert width == B
        if packet_id < K:
            return p.source_bytes(K, B)[packet_id * B:(packet_id + 1) * B]
        return bytes((K + attempt + packet_id + i) & 255 for i in range(B))

    def decode(self, profile, source, B, ids, packets, count):
        K, attempt = len(source) // B, profile[28]
        self.calls.append((K, B, attempt, count))
        return int(attempt < self.selected[CELLS.index((K, B))])


def fake_stream(phase, freeze, selection=None, selected=None):
    native, rows = FakeNative(freeze["history"], selected), []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        c.write_json(root / "freeze.json", freeze)
        selection_path = root / "selection.json" if selection is not None else None
        if selection is not None:
            c.write_json(selection_path, selection)
        with mock.patch.object(p, "Native", return_value=native), \
                mock.patch.object(w.j, "check_build"), mock.patch.object(p, "emit", side_effect=rows.append):
            w.worker(phase, root / "freeze.json", 60, selection_path,
                     c.file_digest(selection_path) if selection_path else None)
    assert native.opened == native.closed
    return rows, native


class WidthLocalProtocolTest(unittest.TestCase):
    def test_private_import_frozen_roots_and_six_entry_selection(self):
        self.assertIsNot(w.j, previous)
        self.assertIsNot(c, previous.c)
        self.assertEqual((previous.PROTOCOL, previous.WIDTHS, previous.TRAINING_ROOTS,
                          previous.HOLDOUT_ROOTS, previous.c.PROTOCOL), PREVIOUS_STATE)
        self.assertEqual(tuple(map(tuple, w.CELLS)), CELLS)
        roots = {phase: tuple("0x" + hashlib.sha256((w.PROTOCOL + ":" + phase + "/" + str(i))
                                                  .encode("ascii")).hexdigest()[:16]
                              for i in range(count)) for phase, count in (("train", 16), ("holdout", 64))}
        self.assertEqual(w.TRAINING_ROOTS, roots["train"])
        self.assertEqual(w.HOLDOUT_ROOTS, roots["holdout"])
        self.assertEqual(c.digest(c.canonical(roots)),
                         "2297d69d6b7a5a83f715589f981deea0d48a676599e02e5697d8ecff6dcfa89f")
        self.assertEqual(len(set(roots["train"] + roots["holdout"])), 80)
        selection = dict(protocol=w.PROTOCOL, selected_attempts=[0, 1, 2, 3, 4, 255], search_sha256="a" * 64)
        self.assertEqual(w.j.validate_selection(selection), selection["selected_attempts"])
        for bad in ([5, 3], [0] * 5, [0] * 7, [0, 0, 0, 0, 0, True], [0, 0, 0, 0, 0, 256]):
            with self.subTest(map=bad), self.assertRaises(c.ValidationError):
                w.j.validate_selection(dict(selection, selected_attempts=bad))

    @unittest.skipUnless(history_available(), "machine-local immutable history fixtures unavailable")
    def test_history_extension_no_omission_and_exact_reserved_roots(self):
        before = {str(path / "COMPLETE"): c.file_digest(path / "COMPLETE") for path in HISTORY_BUNDLES}
        with mock.patch.object(p, "Native", side_effect=AssertionError("native history call")), \
                mock.patch.object(c, "command", side_effect=AssertionError("history subprocess")):
            history = w.load_history()
            w.validate_history(history)
        expected = {
            "prefixes": (75, "6f1fac82e97e72c124b83c495316eca69acb4180097f8ea60e05674404fca33b"),
            "origins": (212, "89e514a5b6a673c0fe5a5188503c0f125717eb35962f48cd7bfa55a3646d8d39"),
            "bad_origins": (786, "2cca4560caa79b35f7d5103a694231e22d2797a4f466109fb1b3f9dce3241a54"),
            "known_bad": (426, "1438930c6732e27b8f9a34ce7a8b788f441d9a4f25400fb109ce9b38cfb218dd"),
            "prior_roots": (207, "5b89e3a62123156bd8a4918f34d1364e52dd9e0d091300a302d944c344c3997d"),
            "consumed_roots": (143, "ada764d2cb8afa51e9ea13c130e8536ae0b2d23a4ad1209b532987739c1b291e"),
            "unconsumed_holdout_roots": (64, "85e942184cb461bd91914c9f97e293d13daf8564413af97a9440ff86af641a65")}
        for key, (count, digest) in expected.items():
            self.assertEqual((len(history[key]), c.digest(c.canonical(history[key]))), (count, digest))
            changed = copy.deepcopy(history)
            changed[key] = changed[key][:-1]
            with self.subTest(ledger=key), self.assertRaises(c.ValidationError):
                w.validate_history(changed)
        self.assertEqual([sum(row["K"] == K for row in history["prefixes"]) for K in (3, 6)], [24, 51])
        self.assertFalse(set(history["prior_roots"]) & set(w.TRAINING_ROOTS + w.HOLDOUT_ROOTS))
        self.assertEqual(before, {path: c.file_digest(path) for path in before})

    def test_inherited_supervisor_enforces_stream_and_combined_phase_caps(self):
        def child(size):
            return [sys.executable, "-I", "-B", "-S", "-c", "import os; os.write(1,b'x'*%d)" % size]

        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(w.j, "STREAM_LIMIT", 16), mock.patch.object(w.j, "TOTAL_LIMIT", 32):
            root, total = Path(directory), 0
            for phase, size in (("search", 15), ("validate", 15), ("holdout", 2)):
                total = w.j.execute_worker(child(size), root / phase, root / (phase + ".err"), 5, total)
            self.assertEqual(total, 32)
            for name, size, prior in (("stream", 16, 0), ("aggregate", 1, 32)):
                with self.subTest(cap=name), self.assertRaisesRegex(c.ValidationError, "raw byte cap"):
                    w.j.execute_worker(child(size), root / name, root / (name + ".err"), 5, prior)
                self.assertLess((root / name).stat().st_size, 16)
                self.assertLessEqual((root / name).stat().st_size + prior, 32)


@unittest.skipUnless(history_available(), "machine-local immutable history fixtures unavailable")
class WidthLocalStreamTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.history = w.load_history()
        cls.freeze = w.frozen_protocol(synthetic_receipt(), cls.history)
        cls.search_rows, cls.native = fake_stream("search", cls.freeze)
        cls.search_text = encoded(cls.search_rows)
        cls.search = w.parse_search(cls.search_text.decode("ascii"), cls.freeze)
        cls.selection = dict(protocol=w.PROTOCOL, selected_attempts=list(cls.native.selected),
                             search_sha256=c.digest(cls.search_text))
        cls.validation_rows, _ = fake_stream("validate", cls.freeze, cls.selection)
        cls.holdout_rows, _ = fake_stream("holdout", cls.freeze, cls.selection)
        cls.validation, cls.holdout = cls.parse(cls.validation_rows, "validate"), cls.parse(cls.holdout_rows, "holdout")
        cls.files = {"freeze.json": cls.freeze, "selection.json": cls.selection,
                     "search.jsonl": cls.search_text, "search.stderr": b"",
                     "validate.jsonl": encoded(cls.validation_rows), "validate.stderr": b"",
                     "holdout.jsonl": encoded(cls.holdout_rows), "holdout.stderr": b"",
                     "summary.json": scoped(w.summarize(cls.search, cls.validation, cls.holdout))}

    @classmethod
    def parse(cls, rows, phase):
        return w.parse_validation(encoded(rows).decode("ascii"), phase, cls.freeze, cls.selection, cls.search)

    def test_all_six_choices_revalidate_with_same_kb_baselines(self):
        self.assertEqual(self.search["selected_attempts"], list(self.native.selected))
        self.assertEqual(self.validation["baseline_attempts"], list(self.native.baselines))
        self.assertEqual(self.holdout["baseline_attempts"], list(self.native.baselines))
        self.assertEqual((len(self.validation["rows"]), len(self.holdout["rows"])), (1026, 2304))
        self.assertEqual(self.files["summary.json"]["outcome"], "PASS")
        self.assertEqual(self.files["summary.json"]["total_prefix_cases"], self.search["prefix_cases"] + 11970)
        self.assertEqual((self.freeze["max_prefix_cases"], self.freeze["max_constructors"]), (143298, 1560))
        self.assertEqual(self.search["constructors"], sum(x + 1 for x in self.native.selected))
        expected = [(K, B, attempt) for (K, B), choice in zip(CELLS, self.native.selected)
                    for attempt in range(choice + 1)]
        self.assertEqual(self.native.opened, expected)
        self.assertTrue(set(expected) & self.native.bad)
        accepted = [row for row in self.search["records"] if row["accepted"]]
        self.assertEqual([row["passed"] for row in accepted], [72, 72, 72, 99, 99, 99])
        self.assertEqual((previous.PROTOCOL, previous.WIDTHS, previous.TRAINING_ROOTS,
                          previous.HOLDOUT_ROOTS, previous.c.PROTOCOL), PREVIOUS_STATE)
        state = (w.j.WIDTHS, w.j.TRAINING_ROOTS, w.j.HOLDOUT_ROOTS, w.j.BAD_PROFILES)
        with self.assertRaisesRegex(RuntimeError, "synthetic oracle"):
            with w.at_width(64, self.history):
                self.assertEqual(w.j.WIDTHS, (64,))
                raise RuntimeError("synthetic oracle")
        self.assertEqual((w.j.WIDTHS, w.j.TRAINING_ROOTS, w.j.HOLDOUT_ROOTS, w.j.BAD_PROFILES), state)

    def test_first_cell_exhaustion_still_searches_later_five_and_seals_full_map(self):
        choices = (-1,) + self.native.selected[1:]
        rows, native = fake_stream("search", self.freeze, selected=choices)
        raw = encoded(rows)
        search = w.parse_search(raw.decode("ascii"), self.freeze)
        self.assertEqual(search["selected_attempts"], list(choices))
        self.assertEqual(len(search["records"]), 256 + sum(x + 1 for x in choices[1:]))
        self.assertEqual(native.opened[:256], [(3, 2, x) for x in range(256)])
        self.assertEqual({(K, B) for K, B, _ in native.opened}, set(CELLS))
        self.assertFalse(any(K == 3 and B == 2 for K, B, _, _ in native.calls))
        files = {"freeze.json": self.freeze, "search.jsonl": raw, "search.stderr": b"",
                 "selection.json": dict(self.selection, selected_attempts=list(choices), search_sha256=c.digest(raw)),
                 "summary.json": scoped(w.summarize(search))}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_bundle(root, files)
            self.assertEqual(w.replay(root)["outcome"], "EXHAUSTED")
            files["validate.jsonl"] = b"forbidden phase\n"
            write_bundle(root, files)
            with self.assertRaises(c.ValidationError):
                w.replay(root)

    def test_search_rejects_missing_cells_attempts_wrong_width_and_bad_profile_contradiction(self):
        for mutation in ("skip", "width", "passed", "known_bad", "map"):
            rows = copy.deepcopy(self.search_rows)
            if mutation == "skip":
                del rows[1]
            elif mutation == "width":
                rows[1]["B"] = 64
            elif mutation == "passed":
                next(row for row in rows[1:-1] if row["witness"])["passed"] += 1
            elif mutation == "known_bad":
                row = next(row for row in rows[1:-1] if (row["K"], row["B"], row["attempt"]) in self.native.bad)
                row["constructors"][0]["outcome"] = 0
            else:
                rows[-1]["selected_attempts"] = list(self.native.selected[:-1])
            rows[-1]["rows"] = len(rows) - 2
            with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                w.parse_search(encoded(rows).decode("ascii"), self.freeze)

    def test_revalidation_profiles_bytes_and_same_cell_baselines_cannot_drift(self):
        for mutation in ("selected_failure", "repair_byte", "selection_hash", "wrong_width"):
            rows = copy.deepcopy(self.validation_rows)
            row = next(row for row in rows if row.get("arm") == "candidate")
            if mutation == "selected_failure":
                row["outcomes"][0], row["recovered_sha256"][0] = 1, None
            elif mutation == "repair_byte":
                index = next(i for i, packet_id in enumerate(row["ids"]) if packet_id >= row["K"])
                payload = bytearray.fromhex(row["packets_hex"])
                payload[index * row["B"]] ^= 1
                row["packets_hex"] = payload.hex()
            elif mutation == "selection_hash":
                rows[0]["selection_sha256"] = "f" * 64
            else:
                row["B"] = 1280
            with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                self.parse(rows, "validate")
        holdout = copy.deepcopy(self.holdout)
        holdout["baseline_attempts"][0] = (holdout["baseline_attempts"][0] + 1) % 256
        with self.assertRaisesRegex(c.ValidationError, "baseline"):
            w.summarize(self.search, self.validation, holdout)

    def test_offline_pass_fail_and_rehashed_omission_helper_map_and_budget_mutations(self):
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(w.j, "load_history", side_effect=AssertionError("live history during replay")), \
                mock.patch.object(c, "command", side_effect=AssertionError("replay subprocess")), \
                mock.patch.object(p, "Native", side_effect=AssertionError("native replay")):
            root = Path(directory)
            write_bundle(root, self.files)
            self.assertEqual(w.replay(root)["outcome"], "PASS")
            rows = copy.deepcopy(self.holdout_rows)
            row = next(row for row in rows if row.get("arm") == "candidate" and row["K"] == 6 and row["B"] == 64)
            row["outcomes"][0], row["recovered_sha256"][0] = 1, None
            failed = dict(self.files, **{"holdout.jsonl": encoded(rows),
                                        "summary.json": scoped(w.summarize(self.search, self.validation,
                                                                           self.parse(rows, "holdout")))})
            write_bundle(root, failed)
            self.assertEqual(w.replay(root)["outcome"], "FAIL")
            for mutation in ("omission", "helper", "map", "summary", "workers", "outer"):
                files = copy.deepcopy(self.files)
                if mutation == "omission":
                    basis = files["freeze.json"]["history"]["basis"]
                    del basis["search"][1]
                    basis["complete"]["search.jsonl"] = c.digest(encoded(basis["search"]))
                elif mutation == "helper":
                    files["freeze.json"]["build"]["source_files"]["bench/Wh2ProductionMix3RecoveryJointR0.py"] = "f" * 64
                elif mutation == "map":
                    files["selection.json"]["selected_attempts"][0] += 1
                elif mutation == "summary":
                    files["summary.json"]["total_prefix_cases"] += 1
                elif mutation == "workers":
                    files["summary.json"]["elapsed_seconds"] = 60
                    for worker in files["summary.json"]["workers"]:
                        worker["elapsed_seconds"] = 20
                else:
                    files["summary.json"]["elapsed_seconds"] = 70
                write_bundle(root, files)
                with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                    w.replay(root)

    def test_controller_seals_full_map_and_separates_worker_and_outer_budgets(self):
        receipt = self.freeze["build"]

        def command(argv, *unused, **kwargs):
            return "" if argv[:2] in (["git", "status"], ["git", "ls-files"]) else receipt["source_commit"]

        real_search, real_validation = w.j.parse_search, w.j.parse_validation
        for invalid_validation, parse_cost in ((False, 0), (True, 0), (False, 11), (False, 31)):
            streams = {phase: self.files[phase + ".jsonl"] for phase in ("search", "validate", "holdout")}
            if invalid_validation:
                rows = copy.deepcopy(self.validation_rows)
                row = next(row for row in rows if row.get("arm") == "candidate")
                row["outcomes"][0], row["recovered_sha256"][0] = 1, None
                streams["validate"] = encoded(rows)
            seen, total, now, timeouts = [], [0], [0.0], []

            def timed(parser):
                def invoke(*args):
                    result = parser(*args)
                    now[0] += parse_cost
                    return result
                return invoke

            def transport(argv, raw, error, timeout, prior_bytes=0):
                phase = argv[argv.index("--worker") + 1]
                self.assertEqual(phase, ("search", "validate", "holdout")[len(seen)])
                self.assertEqual(prior_bytes, total[0])
                self.assertEqual(argv[1:4], ["-I", "-B", "-S"])
                timeouts.append(timeout)
                now[0] += 10 if parse_cost else 0
                if phase == "search":
                    self.assertNotIn("--selection", argv)
                else:
                    path = Path(argv[argv.index("--selection") + 1])
                    self.assertEqual(c.parse_json(path.read_bytes()), self.selection)
                    self.assertEqual(c.file_digest(path), argv[argv.index("--selection-sha256") + 1])
                Path(raw).write_bytes(streams[phase])
                Path(error).write_bytes(b"")
                seen.append(phase)
                total[0] += len(streams[phase])
                return total[0]

            with self.subTest(invalid=invalid_validation, parsing=parse_cost), tempfile.TemporaryDirectory() as directory:
                root, bundle = Path(directory), Path(directory) / "bundle"
                c.write_json(root / "build.json", receipt)
                c.write_json(root / "inputs.json", self.history)
                with mock.patch.object(w.j, "NAMESPACE", bundle), mock.patch.object(w.j, "check_build"), \
                        mock.patch.object(w.j, "load_history", return_value=self.history), \
                        mock.patch.object(c, "command", side_effect=command), \
                        mock.patch.object(w.j, "execute_worker", side_effect=transport), \
                        mock.patch.object(w.j.time, "monotonic", side_effect=lambda: now[0]), \
                        mock.patch.object(w.j, "parse_search", side_effect=timed(real_search)), \
                        mock.patch.object(w.j, "parse_validation", side_effect=timed(real_validation)), \
                        mock.patch.object(p, "Native", side_effect=AssertionError("native controller call")):
                    result = w.run_once(root / "build.json", bundle)
                self.assertEqual(w.replay(bundle), result)
                invalid = invalid_validation or parse_cost == 31
                self.assertEqual(result["outcome"], "INVALID" if invalid else "PASS")
                self.assertEqual(seen, ["search", "validate"] if invalid else ["search", "validate", "holdout"])
                if invalid:
                    self.assertFalse((bundle / "holdout.jsonl").exists())
                if parse_cost:
                    self.assertEqual(timeouts, [60, 49, 28] if parse_cost == 11 else [60, 29])
                    self.assertEqual([row["elapsed_seconds"] for row in result["workers"]], [10] * len(seen))

    def test_publication_deadline_seals_only_invalid_and_retains_completed_phases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, value in self.files.items():
                if name != "summary.json":
                    (root / name).write_bytes(value if isinstance(value, bytes) else c.canonical(value) + b"\n")
            now, real_digest = [0.0], c.file_digest

            def delayed_digest(path):
                now[0] = 70.0
                return real_digest(path)

            with mock.patch.object(w.j.time, "monotonic", side_effect=lambda: now[0]), \
                    mock.patch.object(c, "file_digest", side_effect=delayed_digest):
                result = w.j.publish_result(root, self.files["summary.json"], 0.0)
            self.assertEqual(result["outcome"], "INVALID")
            self.assertEqual(result["completed_phases"], ["search", "validate", "holdout"])
            self.assertEqual(w.replay(root), result)


if __name__ == "__main__":
    unittest.main()
