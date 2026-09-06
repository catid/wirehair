#!/usr/bin/env python3
"""Neutral controller tests. No compiler, candidate, or fixed namespace run."""
import contextlib
import copy
import importlib.util
import io
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("tm_native_launcher_test", HERE / "Wh2ThueMorseNativeRunR0.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class ParsingTests(unittest.TestCase):
    def test_cache_is_typed_by_exact_text_and_rejects_duplicates(self):
        self.assertEqual(M.parse_cache(b"# comment\n//x\nA:BOOL=OFF\nB:STRING=\n"), {"A": "OFF", "B": ""})
        for invalid in (b"A:BOOL=ON\nA:BOOL=OFF", b"=x", b"bad"):
            with self.assertRaises(ValueError):
                M.parse_cache(invalid)

    def test_compiler_dependency_closure_rejects_stale_or_missing_entries(self):
        valid = b"file.o: #deps 2, deps mtime 123 (VALID)\n    /usr/include/a.h\n    source.cpp\n\n"
        self.assertEqual(M.dependency_paths(valid), ["/usr/include/a.h", "source.cpp"])
        for invalid in (valid.replace(b"VALID", b"STALE"), b"", b"    \n", b"garbage"):
            with self.assertRaises(ValueError):
                M.dependency_paths(invalid)

    def test_single_shared_runtime_symbol_definitions(self):
        names = ("GF256Ctx", "gf256_init_", "gf256_mulset_multi_mem", "gf256_muladd_mem",
                 "wirehair_encode", "wirehair_decode", "wirehair_v2_encode", "wirehair_v2_decode")
        valid = "\n".join("00001 T " + name for name in names).encode()
        M.verify_symbols(valid)
        for invalid in (valid + b"\n00002 T GF256Ctx", valid.replace(b"gf256_init_", b"missing")):
            with self.assertRaises(ValueError):
                M.verify_symbols(invalid)

    def test_unresolved_runtime_rejected(self):
        with self.assertRaises(ValueError):
            M.runtime_paths(b"libbad.so => not found\n")

    def test_exact_compiled_commit_and_language_standard(self):
        commit = "a" * 40
        text = '/usr/bin/g++ -DWIREHAIR_WH2_SOURCE_GIT_COMMIT=\\"' + commit + '\\" -std=gnu++11 -c file.cpp'
        M.verify_compile_commands(text, commit)
        for invalid in (text.replace(commit, "b" * 40), text.replace("gnu++11", "gnu++17"),
                        text.replace("-std=gnu++11", ""), "no compiler"):
            with self.assertRaises(ValueError):
                M.verify_compile_commands(invalid, commit)


class LifecycleTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="wh2-native-controller-neutral-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.receipt_path = self.root / "launch.json"
        self.receipt = {"protocol": M.PROTOCOL, "build_directory": str(self.root),
                        "worker_argv": ["NEVER_EXECUTE_NATIVE"], "exact_integer": 1}
        M.C.write_new(self.receipt_path, M.C.canonical(self.receipt))
        self.output = self.root / "result"
        self.patches = (
            mock.patch.object(M, "OUTPUT", self.output),
            mock.patch.object(M, "current_receipt", return_value=self.receipt),
            mock.patch.object(M.C, "capture_worker", return_value={
                "stdout": M.C.canonical({"protocol": M.PROTOCOL}) + b"\n", "stderr": b"",
                "returncode": 0, "failure": None, "elapsed_seconds": 0.001}),
            mock.patch.object(M, "sibling", return_value=SimpleNamespace(
                reduce_raw=lambda *args, **kwargs: {"outcome": "PASS"})),
        )
        self.mocks = [patch.start() for patch in self.patches]
        for patch in self.patches:
            self.addCleanup(patch.stop)

    def launch(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return M.run(self.receipt_path, time.monotonic())

    def summary(self):
        return M.C.strict_json((self.output / "summary.json").read_bytes())

    def test_complete_publication_has_native_claim_and_cannot_rerun(self):
        self.assertEqual(self.launch(), 0)
        summary = self.summary()
        self.assertEqual(summary["outcome"], "PASS")
        self.assertTrue(summary["WH1_compared"])
        self.assertFalse(summary["public_candidate_api_claimed"])
        self.assertFalse(summary["promotion_claimed"])
        complete = M.C.strict_json((self.output / "COMPLETE.json").read_bytes())
        for name, pin in complete["files"].items():
            data = (self.output / name).read_bytes()
            self.assertEqual(pin, {"bytes": len(data), "sha256": M.C.sha(data)})
            self.assertEqual((self.output / name).stat().st_mode & 0o777, 0o400)
        with self.assertRaises(FileExistsError):
            self.launch()
        self.assertEqual(self.mocks[2].call_count, 1)

    def test_preflight_float_integer_alias_never_claims_or_runs(self):
        forged = dict(self.receipt, exact_integer=1.0)
        self.mocks[1].return_value = forged
        with self.assertRaisesRegex(ValueError, "launch receipt changed"):
            self.launch()
        self.assertFalse(self.output.exists())
        self.mocks[2].assert_not_called()

    def test_post_pin_bool_alias_is_invalid_and_retains_evidence(self):
        self.mocks[1].side_effect = [self.receipt, dict(self.receipt, exact_integer=True)]
        self.assertEqual(self.launch(), 1)
        self.assertEqual(self.summary()["outcome"], "INVALID")
        self.assertEqual(self.summary()["worker"]["post_pins"], "failed_or_unobserved")
        self.assertTrue((self.output / "raw.json").exists())

    def test_worker_failure_and_malformed_json_are_terminal(self):
        self.mocks[2].return_value["returncode"] = 7
        self.assertEqual(self.launch(), 1)
        self.assertEqual(self.summary()["outcome"], "INVALID")
        self.assertFalse(self.summary()["WH1_compared"])

    def test_control_fail_is_reported_not_relaxed(self):
        self.mocks[3].return_value.reduce_raw = lambda *a, **kw: {"outcome": "CONTROL_FAIL"}
        self.assertEqual(self.launch(), 0)
        self.assertEqual(self.summary()["outcome"], "CONTROL_FAIL")

    def test_duplicate_json_key_cannot_reach_reducer(self):
        self.mocks[2].return_value["stdout"] = b'{"protocol":"a","protocol":"b"}\n'
        self.assertEqual(self.launch(), 1)
        self.assertIn("duplicate JSON", self.summary()["failure"])


if __name__ == "__main__":
    unittest.main()
