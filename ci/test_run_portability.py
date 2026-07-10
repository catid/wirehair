#!/usr/bin/env python3
"""Unit tests for the architecture-specific portability runner."""

from pathlib import Path
import tempfile
import unittest

try:
    from . import run_portability
except ImportError:
    import run_portability


class PortabilityRunnerTests(unittest.TestCase):
    def test_cross_targets_pin_toolchains_emulators_and_runtime_contracts(self):
        self.assertEqual(
            {"aarch64", "i686", "riscv64", "s390x"},
            set(run_portability.CROSS_TARGETS))
        self.assertEqual(
            "qemu-s390x",
            run_portability.CROSS_TARGETS["s390x"]["emulator"])
        self.assertEqual(
            "big endian",
            run_portability.CROSS_TARGETS["s390x"]["elf_data"])
        self.assertEqual(
            "ELF32", run_portability.CROSS_TARGETS["i686"]["elf_class"])
        self.assertIn(
            "allocation_overflow_test",
            run_portability.CROSS_TARGETS["i686"]["targets"])
        self.assertEqual(
            "scalar", run_portability.CROSS_TARGETS["riscv64"]["backend"])
        self.assertIn(
            "gf256_inplace_test",
            run_portability.CROSS_TARGETS["aarch64"]["targets"])

    def test_pinned_version_families_are_bounded(self):
        self.assertIsNotNone(run_portability.GCC_VERSION_RE.search("13.3.0"))
        self.assertIsNone(run_portability.GCC_VERSION_RE.search("14.2.0"))
        self.assertIsNotNone(
            run_portability.QEMU_VERSION_RE.search("qemu-x version 8.2.2"))
        self.assertIsNone(
            run_portability.QEMU_VERSION_RE.search("qemu-x version 9.0.0"))

    def test_validate_gf256_backend_distinguishes_neon_scalar_and_auxv(self):
        run_portability.validate_gf256_backend(
            {"GF256_TARGET_MOBILE": "", "GF256_TRY_NEON": "",
             "LINUX_ARM": ""},
            "neon", linux_arm=True)
        run_portability.validate_gf256_backend(
            {"GF256_TARGET_MOBILE": ""},
            "scalar", linux_arm=False)
        with self.assertRaisesRegex(RuntimeError, "backend mismatch"):
            run_portability.validate_gf256_backend(
                {"GF256_TARGET_MOBILE": ""},
                "neon", linux_arm=False)
        with self.assertRaisesRegex(RuntimeError, "auxv guard"):
            run_portability.validate_gf256_backend(
                {"GF256_TARGET_MOBILE": "", "LINUX_ARM": ""},
                "scalar", linux_arm=False)

    def test_validate_elf_header_checks_class_endian_and_machine(self):
        header = """
          Class:                             ELF32
          Data:                              2's complement, little endian
          Machine:                           Intel 80386
        """
        run_portability.validate_elf_header(
            header, run_portability.CROSS_TARGETS["i686"])
        with self.assertRaisesRegex(RuntimeError, "architecture mismatch"):
            run_portability.validate_elf_header(
                header, run_portability.CROSS_TARGETS["s390x"])

    def test_macho_exports_strip_one_abi_underscore_and_reject_extras(self):
        parsed = run_portability.parse_macho_exports(
            "_wirehair_decode\n_wirehair_free\n")
        self.assertEqual(
            frozenset({"wirehair_decode", "wirehair_free"}), parsed)
        self.assertIn(
            "_accidental_cpp_export",
            run_portability.parse_macho_exports(
                "_wirehair_decode\n_accidental_cpp_export\n"))

    def test_export_map_parser_accepts_versioned_policy_and_rejects_wildcard(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "wirehair.map"
            path.write_text(
                "WIREHAIR_2.0 { global: wirehair_decode; "
                "wirehair_free; local: *; };\n",
                encoding="ascii")
            self.assertEqual(
                frozenset({"wirehair_decode", "wirehair_free"}),
                run_portability.load_public_exports(path))
            path.write_text(
                "WIREHAIR_2.0 { global: *; local: *; };\n",
                encoding="ascii")
            with self.assertRaisesRegex(RuntimeError, "malformed"):
                run_portability.load_public_exports(path)

    def test_parser_exposes_native_and_cross_lanes_with_xcode_pin(self):
        cross = run_portability.parse_args([
            "cross", "--target", "riscv64", "--jobs", "3"])
        self.assertIs(cross.function, run_portability.run_cross)
        self.assertEqual("riscv64", cross.target)
        self.assertEqual(3, cross.jobs)
        macos = run_portability.parse_args(["macos-arm64"])
        self.assertIs(macos.function, run_portability.run_macos_arm64)
        self.assertEqual("16.4", macos.xcode_version)
        self.assertEqual("16F6", macos.xcode_build)

    def test_workflows_route_fast_and_scheduled_architecture_budgets(self):
        per_change = (
            run_portability.ROOT / ".github" / "workflows" / "ci.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("runs-on: ubuntu-24.04-arm", per_change)
        self.assertIn("run_portability.py native-arm64", per_change)
        self.assertIn("run_portability.py cross --target i686", per_change)
        self.assertIn("timeout-minutes: 25", per_change)
        self.assertIn("timeout-minutes: 20", per_change)

        scheduled = (
            run_portability.ROOT / ".github" / "workflows" / "scheduled.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("run_portability.py cross --target s390x", scheduled)
        self.assertIn("run_portability.py cross --target riscv64", scheduled)
        self.assertIn("runs-on: macos-15", scheduled)
        self.assertIn("/Applications/Xcode_16.4.app", scheduled)


if __name__ == "__main__":
    unittest.main()
