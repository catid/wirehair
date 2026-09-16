"""Synthetic mutation tests for the isolated K4 R1 qualifier.

These tests never invoke CMake, a compiler, a codec worker, or a scientific
namespace.  They exercise the fail-closed parsers and receipt bindings.
"""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import Wh2K4FreshQualifierR1 as Q


class QualifierTest(unittest.TestCase):
    def test_pin_map_rejects_schema_and_conflicts(self):
        record = {"path": "/synthetic/a", "bytes": 1, "sha256": "a" * 64}
        self.assertEqual(Q.pin_map([record, copy.deepcopy(record)]),
                         {record["path"]: record})
        for changed in (dict(record, bytes=True), dict(record, sha256="A" * 64),
                        dict(record, path="relative")):
            with self.assertRaises(ValueError):
                Q.pin_map([changed])
        with self.assertRaises(ValueError):
            Q.pin_map([record, dict(record, bytes=2)])

    def test_depfile_requires_exact_target_and_absolute_inputs(self):
        with tempfile.TemporaryDirectory(prefix="wh2-k4-r1-dep-") as directory:
            root = Path(directory)
            target = root / "x.o"
            source = root / "x.cpp"
            target.write_bytes(b"o")
            source.write_bytes(b"s")
            raw = (str(target) + ": " + str(target) + " " + str(source) + "\n").encode()
            self.assertEqual(Q._dep_paths(raw, target), {target, source})
            for mutated in (raw.replace(str(target).encode(), b"wrong", 1),
                            raw.replace(str(source).encode(), b"relative", 1)):
                with self.assertRaises(ValueError):
                    Q._dep_paths(mutated, target)

    def test_compile_recipe_rejects_lto_scalar_and_extra_flags(self):
        source = Q.ROOT / Q.PRODUCERS[0]
        library = Path("/var/tmp/wh2-k4-fresh-neutral-r1.synthetic/native/library")
        prefix = "CMakeFiles/wirehair_objects.dir/"
        output = prefix + source.relative_to(Q.ROOT).as_posix() + ".o"
        entry = {"directory": str(library), "file": str(source), "output": output,
                 "command": " ".join(["/usr/bin/c++"] + Q._compiler_flags("native") +
                                      ["-o", output, "-c", str(source)])}
        # The path is synthetic, so only recipe parsing is exercised here.
        with patch.object(Q, "exact", side_effect=Q.exact):
            Q._validate_compile_command(entry, source, library / output, library,
                                        "native", prefix)
        for extra in ("-flto", "-DANDROID", "-fno-lto"):
            changed = copy.deepcopy(entry)
            changed["command"] = changed["command"].replace(
                "-o " + output, extra + " -o " + output)
            with self.assertRaises(ValueError):
                Q._validate_compile_command(changed, source, library / output,
                                            library, "native", prefix)

    def test_artifact_roster_is_exact_and_rehashed(self):
        with tempfile.TemporaryDirectory(prefix="wh2-k4-r1-artifacts-") as directory:
            root = Path(directory) / "native"
            root.mkdir()
            (root / "SOURCE.json").write_bytes(b"source")
            (root / "command-000.stdout").write_bytes(b"stdout")
            records = [Q.pin(root / "SOURCE.json"), Q.pin(root / "command-000.stdout")]
            result = {"artifacts": records}
            Q._validate_artifacts(root, result)
            (root / "extra").write_bytes(b"unexpected")
            with self.assertRaises(ValueError):
                Q._validate_artifacts(root, result)
            (root / "extra").unlink()
            (root / "command-000.stdout").write_bytes(b"changed")
            with self.assertRaises(ValueError):
                Q._validate_artifacts(root, result)

    def test_proof_write_is_after_validation(self):
        # A rejected mode root must not create a destination or any evidence.
        with tempfile.TemporaryDirectory(prefix="wh2-k4-r1-proof-") as directory:
            mode = Path(directory) / "native"
            mode.mkdir()
            proof = Path(directory) / "proof.json"
            with self.assertRaises(ValueError):
                Q.qualify(mode, "native", proof)
            self.assertFalse(proof.exists())


if __name__ == "__main__":
    unittest.main()
