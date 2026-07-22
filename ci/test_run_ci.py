#!/usr/bin/env python3
"""Unit tests for the cross-platform CI runner's path and environment logic."""

from contextlib import redirect_stderr
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

try:
    from . import python_package_test, run_ci
except ImportError:
    import python_package_test
    import run_ci


class RunnerTests(unittest.TestCase):
    def test_python_package_work_directory_cannot_remove_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            source = root / "source"
            prefix = root / "native-prefix"
            library = prefix / "lib" / "libwirehair.so"
            source.mkdir()
            library.parent.mkdir(parents=True)
            library.write_bytes(b"native")
            args = mock.Mock(
                source_dir=source,
                library=library,
                native_prefix=prefix,
                work_dir=root,
            )
            with mock.patch.object(
                    python_package_test, "parse_args", return_value=args), \
                    mock.patch.object(
                        python_package_test.shutil, "rmtree") as remove, \
                    self.assertRaisesRegex(
                        RuntimeError, "refusing to remove work directory"):
                python_package_test.main()
            remove.assert_not_called()

    def test_find_executable_supports_single_and_multi_config_layouts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            suffix = ".exe" if os.name == "nt" else ""
            single = root / ("unit_test" + suffix)
            single.write_bytes(b"")
            self.assertEqual(single.resolve(), run_ci.find_executable(root, "unit_test"))
            single.unlink()

            multi = root / "Debug" / ("unit_test" + suffix)
            multi.parent.mkdir()
            multi.write_bytes(b"")
            self.assertEqual(
                multi.resolve(), run_ci.find_executable(root, "unit_test", "Debug"))

    def test_find_installed_python_and_shared_libraries(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            module = prefix / "python" / "whirehair.py"
            module.parent.mkdir(parents=True)
            module.write_text("", encoding="ascii")
            linux_library = prefix / "lib" / "libwirehair.so.2"
            linux_library.parent.mkdir()
            linux_library.write_bytes(b"")
            windows_library = prefix / "bin" / "wirehair.dll"
            windows_library.parent.mkdir()
            windows_library.write_bytes(b"")

            self.assertEqual(module.resolve(), run_ci.find_python_module(prefix))
            self.assertEqual(
                linux_library.resolve(),
                run_ci.find_shared_library(prefix, platform="linux"),
            )
            self.assertEqual(
                windows_library.resolve(),
                run_ci.find_shared_library(prefix, platform="win32"),
            )

    def test_runtime_environment_prepends_installed_library_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary).resolve()
            with mock.patch.dict(os.environ, {"LD_LIBRARY_PATH": "/prior"}, clear=False):
                linux = run_ci.runtime_env(prefix, platform="linux")
            self.assertEqual(
                [str(prefix / "lib"), str(prefix / "lib64"), "/prior"],
                linux["LD_LIBRARY_PATH"].split(os.pathsep),
            )

            with mock.patch.dict(os.environ, {"PATH": "prior"}, clear=False):
                windows = run_ci.runtime_env(prefix, platform="win32")
            self.assertEqual(
                [str(prefix / "bin"), str(prefix / "lib"), "prior"],
                windows["PATH"].split(os.pathsep),
            )

    def test_parser_rejects_unknown_linkage_and_accepts_extra_cmake_args(self):
        parsed = run_ci.parse_args([
            "matrix", "--linkage", "shared", "--cmake-arg=-DTESTING=ON",
            "--repository-gates", "--tool-smoke",
        ])
        self.assertEqual("shared", parsed.linkage)
        self.assertEqual(["-DTESTING=ON"], parsed.cmake_arg)
        self.assertTrue(parsed.repository_gates)
        self.assertTrue(parsed.tool_smoke)
        defaults = run_ci.parse_args(["matrix", "--linkage", "static"])
        self.assertFalse(defaults.repository_gates)
        self.assertFalse(defaults.tool_smoke)
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            run_ci.parse_args(["matrix", "--linkage", "both"])

    def test_consumer_forwards_only_selected_compilers(self):
        args = mock.Mock(
            cmake_arg=[
                "-DFEATURE=ON",
                "-DCMAKE_C_COMPILER=clang",
                "-DCMAKE_CXX_COMPILER=clang++",
            ],
            strict=False,
            generator="Ninja",
        )
        self.assertEqual(
            ["-DCMAKE_C_COMPILER=clang", "-DCMAKE_CXX_COMPILER=clang++"],
            run_ci.consumer_toolchain_args(args),
        )

    def test_strict_flags_cover_unix_and_msvc_c_and_cxx(self):
        unix = mock.Mock(strict=True, generator="Ninja")
        self.assertEqual(
            ["-DWIREHAIR_STRICT_WARNINGS=ON"],
            run_ci.strict_cmake_args(unix),
        )
        msvc = mock.Mock(strict=True, generator="Visual Studio 17 2022")
        self.assertEqual(
            ["-DWIREHAIR_STRICT_WARNINGS=ON"],
            run_ci.strict_cmake_args(msvc),
        )
        self.assertFalse(any(
            argument.startswith(("-DCMAKE_C_FLAGS=", "-DCMAKE_CXX_FLAGS="))
            for argument in run_ci.strict_cmake_args(msvc)
        ))

    def test_build_captures_output_in_failure_artifact_log(self):
        args = mock.Mock(
            build_dir=Path("build/ci-test"),
            config="Release",
            jobs=2,
        )
        with mock.patch.object(run_ci, "run") as run:
            run_ci.build(args, targets=("unit_test",))
        run.assert_called_once_with(
            [
                "cmake", "--build", args.build_dir,
                "--config", "Release", "--parallel", "2",
                "--target", "unit_test",
            ],
            log_path=args.build_dir / "ci-logs" / "build.log",
        )

    def test_ctest_parallelism_and_exclusion_keep_failure_output(self):
        args = mock.Mock(
            build_dir=Path("build/ci-test"),
            config="Release",
            jobs=7,
        )
        with mock.patch.object(run_ci, "run") as run:
            run_ci.ctest(args, exclude_label=run_ci.CI_ONCE_LABEL)
        run.assert_called_once_with([
            "ctest", "--test-dir", args.build_dir,
            "--build-config", "Release",
            "--output-on-failure", "--no-tests=error",
            "--parallel", "7",
            "--label-exclude", run_ci.CI_ONCE_LABEL,
        ], env=None)

    def test_ci_once_inventory_requires_exact_owned_gate_set(self):
        expected = set(run_ci.CI_ONCE_PACKAGE_TESTS)
        if os.name != "nt":
            expected.add("build_policy_e2e")

        def inventory(names):
            return json.dumps({
                "tests": [
                    {
                        "name": name,
                        "properties": [{
                            "name": "LABELS",
                            "value": ["ci-once"],
                        }],
                    }
                    for name in names
                ],
            })

        args = mock.Mock(build_dir=Path("build/test"), config="Release")
        with mock.patch.object(
                run_ci, "capture_output", return_value=inventory(expected)):
            run_ci.validate_ci_once_inventory(args)

        with mock.patch.object(run_ci.os, "name", "nt"), \
                mock.patch.object(
                    run_ci, "capture_output",
                    return_value=inventory(run_ci.CI_ONCE_PACKAGE_TESTS)):
            run_ci.validate_ci_once_inventory(args)

        for names in (expected - {next(iter(expected))},
                      expected | {"accidental_duplicate_gate"}):
            with self.subTest(names=names), mock.patch.object(
                    run_ci, "capture_output", return_value=inventory(names)), \
                    self.assertRaisesRegex(RuntimeError, "ownership mismatch"):
                run_ci.validate_ci_once_inventory(args)

        with mock.patch.object(
                run_ci, "capture_output", return_value="not JSON"), \
                self.assertRaisesRegex(RuntimeError, "malformed JSON"):
            run_ci.validate_ci_once_inventory(args)

    def test_data_integrity_lane_pins_every_repository_gate(self):
        bash = "/test/bin/bash"
        with mock.patch.object(run_ci.shutil, "which", return_value=bash), \
                mock.patch.object(run_ci, "run") as run:
            run_ci.run_data_integrity(mock.sentinel.args)

        self.assertEqual([
            mock.call([
                sys.executable, "-m", "unittest", "-v",
                "tables.test_dense_count_validate",
                "experiments.precode.test_validate_results",
                "experiments.precode.test_rank_total",
                "experiments.test_validate_byte_metrics",
            ]),
            mock.call([
                sys.executable,
                str(run_ci.ROOT / "tables" / "test_heavy_matrix_docs.py"),
            ]),
            mock.call([
                bash,
                str(run_ci.ROOT / "experiments" / "precode" /
                    "check_results_integrity.sh"),
            ]),
            mock.call([
                bash,
                str(run_ci.ROOT / "experiments" / "test_byte_metrics.sh"),
            ]),
            mock.call([
                bash,
                str(run_ci.ROOT / "experiments" /
                    "test_wh2_thermal_runner.sh"),
            ]),
        ], run.call_args_list)

    def test_data_integrity_parser_is_a_standalone_non_matrix_lane(self):
        parsed = run_ci.parse_args(["data-integrity"])
        self.assertIs(parsed.function, run_ci.run_data_integrity)
        self.assertFalse(hasattr(parsed, "linkage"))
        self.assertFalse(hasattr(parsed, "build_dir"))

    def test_data_integrity_lane_requires_bash_before_running_tests(self):
        with mock.patch.object(run_ci.shutil, "which", return_value=None), \
                mock.patch.object(run_ci, "run") as run, \
                self.assertRaisesRegex(RuntimeError, "bash is required"):
            run_ci.run_data_integrity(mock.sentinel.args)
        run.assert_not_called()

    def test_msvc_strict_matrix_enables_and_smokes_explicit_tools(self):
        args = mock.Mock(
            linkage="static",
            strict=True,
            generator="Visual Studio 17 2022",
            repository_gates=False,
            tool_smoke=True,
        )
        with mock.patch.object(run_ci, "python_unit_tests") as python_tests, \
                mock.patch.object(run_ci, "configure") as configure, \
                mock.patch.object(run_ci, "build"), \
                mock.patch.object(run_ci, "validate_ci_once_inventory"), \
                mock.patch.object(run_ci, "explicit_tools_smoke") as smoke, \
                mock.patch.object(run_ci, "large_message_profiles"), \
                mock.patch.object(run_ci, "ctest") as ctest, \
                mock.patch.object(run_ci, "install"), \
                mock.patch.object(run_ci, "validate_msvc_exports") as exports, \
                mock.patch.object(run_ci, "package_consumer"):
            run_ci.run_matrix(args)
        configure.assert_called_once_with(
            args, shared=False, explicit_tools=True)
        smoke.assert_called_once_with(args)
        python_tests.assert_not_called()
        ctest.assert_called_once_with(
            args, exclude_label=run_ci.CI_ONCE_LABEL)
        exports.assert_not_called()

    def test_repository_gate_owner_runs_python_and_nested_ctests_once(self):
        args = mock.Mock(
            linkage="shared",
            strict=True,
            generator="Ninja",
            repository_gates=True,
            tool_smoke=False,
        )
        with mock.patch.object(run_ci, "python_unit_tests") as python_tests, \
                mock.patch.object(run_ci, "configure") as configure, \
                mock.patch.object(run_ci, "build"), \
                mock.patch.object(run_ci, "validate_ci_once_inventory"), \
                mock.patch.object(run_ci, "explicit_tools_smoke") as smoke, \
                mock.patch.object(run_ci, "large_message_profiles"), \
                mock.patch.object(run_ci, "ctest") as ctest, \
                mock.patch.object(run_ci, "install"), \
                mock.patch.object(run_ci, "validate_msvc_exports"), \
                mock.patch.object(run_ci, "package_consumer"), \
                mock.patch.object(
                    run_ci, "python_distribution_test") as distribution, \
                mock.patch.object(run_ci, "python_native_test") as native:
            run_ci.run_matrix(args)
        python_tests.assert_called_once_with()
        configure.assert_called_once_with(
            args, shared=True, explicit_tools=False)
        smoke.assert_not_called()
        ctest.assert_called_once_with(args)
        distribution.assert_called_once_with(args)
        native.assert_not_called()

    def test_repository_gate_owner_requires_loadable_shared_library(self):
        args = mock.Mock(
            linkage="static",
            repository_gates=True,
            tool_smoke=False,
        )
        with mock.patch.object(run_ci, "python_unit_tests") as python_tests, \
                self.assertRaisesRegex(RuntimeError, "shared-library cell"):
            run_ci.run_matrix(args)
        python_tests.assert_not_called()

    def test_python_distribution_gate_uses_installed_library_and_logs(self):
        args = mock.Mock(
            build_dir=Path("build/ci-gcc-shared"),
            install_dir=Path("build/install-gcc-shared"),
        )
        library = Path("build/install-gcc-shared/lib/libwirehair.so")
        with mock.patch.object(
                run_ci, "find_shared_library", return_value=library), \
                mock.patch.object(run_ci, "run") as run:
            run_ci.python_distribution_test(args)
        run.assert_called_once_with([
            sys.executable,
            str(run_ci.ROOT / "ci" / "python_package_test.py"),
            "--source-dir", run_ci.ROOT,
            "--library", library,
            "--native-prefix", args.install_dir,
            "--work-dir", args.build_dir / "python-package-e2e",
        ], log_path=args.build_dir / "ci-logs" / "python-package-e2e.log")

    def test_python_native_gate_checks_automatic_then_explicit_discovery(self):
        args = mock.Mock(install_dir=Path("build/install-gcc-shared"))
        module = args.install_dir / "python" / "whirehair.py"
        library = args.install_dir / "lib" / "libwirehair.so"
        environment = {"LD_LIBRARY_PATH": str(args.install_dir / "lib")}
        with mock.patch.object(
                run_ci, "find_python_module", return_value=module), \
                mock.patch.object(
                    run_ci, "find_shared_library", return_value=library), \
                mock.patch.object(
                    run_ci, "runtime_env", return_value=environment), \
                mock.patch.object(run_ci, "run") as run:
            run_ci.python_native_test(args)
        self.assertEqual(2, run.call_count)
        automatic, explicit = run.call_args_list
        self.assertEqual([sys.executable, "-I", "-c"], automatic.args[0][:3])
        self.assertEqual(module.parent, automatic.args[0][-2])
        self.assertEqual(library, automatic.args[0][-1])
        self.assertNotIn("WIREHAIR_LIBRARY", automatic.kwargs["env"])
        self.assertEqual(
            str(library), explicit.kwargs["env"]["WIREHAIR_LIBRARY"])

    def test_msvc_strict_shared_matrix_validates_installed_dll_exports(self):
        args = mock.Mock(
            linkage="shared",
            strict=True,
            generator="Visual Studio 17 2022",
            repository_gates=False,
            tool_smoke=False,
        )
        with mock.patch.object(run_ci, "python_unit_tests"), \
                mock.patch.object(run_ci, "configure"), \
                mock.patch.object(run_ci, "build"), \
                mock.patch.object(run_ci, "validate_ci_once_inventory"), \
                mock.patch.object(run_ci, "explicit_tools_smoke"), \
                mock.patch.object(run_ci, "large_message_profiles"), \
                mock.patch.object(run_ci, "ctest"), \
                mock.patch.object(run_ci, "install"), \
                mock.patch.object(run_ci, "validate_msvc_exports") as exports, \
                mock.patch.object(run_ci, "package_consumer"), \
                mock.patch.object(run_ci, "python_native_test"):
            run_ci.run_matrix(args)
        exports.assert_called_once_with(args)

    def test_msvc_export_gate_selection_is_strict_native_shared_only(self):
        for linkage, strict, generator, expected in (
                ("shared", True, "Visual Studio 17 2022", True),
                ("static", True, "Visual Studio 17 2022", False),
                ("shared", False, "Visual Studio 17 2022", False),
                ("shared", True, "Ninja", False)):
            with self.subTest(
                    linkage=linkage, strict=strict, generator=generator):
                args = mock.Mock(
                    linkage=linkage, strict=strict, generator=generator)
                self.assertEqual(
                    expected, run_ci.should_validate_msvc_exports(args))

    def test_non_msvc_matrix_keeps_explicit_tools_out_of_default_build(self):
        args = mock.Mock(
            linkage="static",
            strict=True,
            generator="Ninja",
            repository_gates=False,
            tool_smoke=False,
        )
        with mock.patch.object(run_ci, "python_unit_tests"), \
                mock.patch.object(run_ci, "configure") as configure, \
                mock.patch.object(run_ci, "build"), \
                mock.patch.object(run_ci, "validate_ci_once_inventory"), \
                mock.patch.object(run_ci, "explicit_tools_smoke") as smoke, \
                mock.patch.object(run_ci, "large_message_profiles"), \
                mock.patch.object(run_ci, "ctest"), \
                mock.patch.object(run_ci, "install"), \
                mock.patch.object(run_ci, "validate_msvc_exports") as exports, \
                mock.patch.object(run_ci, "package_consumer"):
            run_ci.run_matrix(args)
        configure.assert_called_once_with(
            args, shared=False, explicit_tools=False)
        smoke.assert_not_called()
        exports.assert_not_called()

    def test_sanitizer_excludes_uninstrumented_nested_builds(self):
        args = mock.Mock()
        with mock.patch.object(run_ci, "configure"), \
                mock.patch.object(run_ci, "build"), \
                mock.patch.object(run_ci, "validate_ci_once_inventory"), \
                mock.patch.object(run_ci, "large_message_profiles"), \
                mock.patch.object(run_ci, "ctest") as ctest:
            run_ci.run_sanitizer(args)
        ctest.assert_called_once_with(
            args,
            exclude_label=run_ci.CI_ONCE_LABEL,
            env=os.environ.copy(),
        )

    def test_workflow_has_single_gate_owners_and_nonduplicate_push_trigger(self):
        workflow = (
            run_ci.ROOT / ".github" / "workflows" / "ci.yml"
        ).read_text(encoding="utf-8")
        self.assertEqual(
            1, workflow.count("repository_args: --repository-gates"))
        self.assertEqual(1, workflow.count("tool_args: --tool-smoke"))
        self.assertEqual(5, workflow.count('repository_args: ""'))
        self.assertEqual(5, workflow.count('tool_args: ""'))
        self.assertIn("      - master\n", workflow)
        self.assertIn("  pull_request:\n", workflow)
        self.assertIn("${{ matrix.repository_args }}", workflow)
        self.assertIn("${{ matrix.tool_args }}", workflow)

    def test_negative_path_requires_failure_and_expected_diagnostic(self):
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "negative.log"
            run_ci.run_expect_failure(
                [sys.executable, "-c", "print('expected'); raise SystemExit(7)"],
                env=os.environ.copy(),
                log_path=log,
                expected_returncode=7,
                required_text="expected",
            )
            self.assertIn("expected", log.read_text(encoding="utf-8"))
            with self.assertRaises(RuntimeError):
                run_ci.run_expect_failure(
                    [sys.executable, "-c", "print('unexpected success')"],
                    env=os.environ.copy(),
                    log_path=log,
                    expected_returncode=7,
                    required_text="expected",
                )
            with self.assertRaises(RuntimeError):
                run_ci.run_expect_failure(
                    [sys.executable, "-c", "print('expected'); raise SystemExit(8)"],
                    env=os.environ.copy(),
                    log_path=log,
                    expected_returncode=7,
                    required_text="expected",
                )
            with self.assertRaises(RuntimeError):
                run_ci.run_expect_failure(
                    [sys.executable, "-c",
                     "print('expected AddressSanitizer'); raise SystemExit(7)"],
                    env=os.environ.copy(),
                    log_path=log,
                    expected_returncode=7,
                    required_text="expected",
                    forbidden_text=("AddressSanitizer",),
                )

    def test_log_line_validation_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "qemu.log"
            log.write_text("prefix expected suffix\nexpected\n", encoding="ascii")
            run_ci.require_log_line(log, "expected")
            with self.assertRaises(RuntimeError):
                run_ci.require_log_line(log, "prefix expected")

    def test_mingw_export_and_import_parsers(self):
        nm_output = "\n".join(
            "00000000 T " + symbol
            for symbol in sorted(run_ci.WIREHAIR_C_API_EXPORTS)
        )
        self.assertEqual(
            run_ci.WIREHAIR_C_API_EXPORTS,
            run_ci.parse_mingw_exports(
                nm_output +
                "\n0000 I __imp_wirehair_decode\n0000 T unrelated"),
        )
        self.assertTrue(run_ci.imports_wirehair_dll(
            "\tDLL Name: KERNEL32.dll\n\tDLL Name: libwirehair.dll\n"))
        self.assertFalse(run_ci.imports_wirehair_dll(
            "\tDLL Name: KERNEL32.dll\n\tDLL Name: libstdc++-6.dll\n"))

    def test_export_map_is_versioned_header_abi_source_of_truth(self):
        version, symbols = run_ci.load_wirehair_export_map()
        self.assertEqual("WIREHAIR_2.0", version)
        self.assertEqual(run_ci.WIREHAIR_ABI_VERSION, version)
        self.assertEqual(run_ci.WIREHAIR_C_API_EXPORTS, symbols)
        self.assertEqual(symbols, run_ci.load_wirehair_header_exports())

    def test_export_map_parser_rejects_policy_broadening_and_drift(self):
        malformed_maps = (
            "WIREHAIR_2.0 { global: wirehair_decode; };",
            "WIREHAIR_2.0 { global: *; local: *; };",
            ("WIREHAIR_2.0 { global: wirehair_decode; "
             "wirehair_decode; local: *; };"),
            ("WIREHAIR_2.0 { global: wirehair_encode; "
             "wirehair_decode; local: *; };"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "wirehair.map"
            for contents in malformed_maps:
                with self.subTest(contents=contents):
                    manifest.write_text(contents, encoding="ascii")
                    with self.assertRaises(RuntimeError):
                        run_ci.load_wirehair_export_map(manifest)

    def test_dumpbin_parser_accepts_only_export_table_rows(self):
        rows = "\n".join(
            "%4d %4X %08X %s" % (ordinal, ordinal - 1, 4096 + ordinal, symbol)
            for ordinal, symbol in enumerate(
                sorted(run_ci.WIREHAIR_C_API_EXPORTS), start=1)
        )
        output = """Microsoft (R) COFF/PE Dumper
wirehair_not_an_export appears in a diagnostic
  ordinal hint RVA      name
%s
  99  63 00009999 unrelated_export
  Summary
""" % rows
        self.assertEqual(
            run_ci.WIREHAIR_C_API_EXPORTS | {"unrelated_export"},
            run_ci.parse_dumpbin_exports(output),
        )
        self.assertEqual(
            frozenset({"wirehair_decode"}),
            run_ci.parse_dumpbin_exports(
                "  1 0 00001000 wirehair_decode = forwarded.decode\n"),
        )
        self.assertEqual(
            frozenset({"wirehair_decode", "wirehair_unexpected"}),
            run_ci.parse_dumpbin_exports(
                "  1 0 00001000 wirehair_decode\n"
                "  2 1 00001010 wirehair_unexpected\n"),
        )

    def test_dumpbin_selection_supports_vs_environment_without_windows(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "VC" / "Tools" / "MSVC" / "14.42.34433"
            dumpbin = tools / "bin" / "Hostx64" / "x64" / "dumpbin.exe"
            dumpbin.parent.mkdir(parents=True)
            dumpbin.write_bytes(b"")
            selected = run_ci.find_dumpbin(
                environment={"VCToolsInstallDir": str(tools)},
                which=lambda unused: None,
            )
            self.assertEqual(dumpbin.resolve(), selected)

            library = root / "wirehair.dll"
            self.assertEqual(
                [str(dumpbin.resolve()), "/nologo", "/exports", str(library)],
                run_ci.dumpbin_exports_command(selected, library),
            )

    def test_dumpbin_selection_prefers_path_and_newest_vs_toolset(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            direct = root / "path" / "dumpbin.exe"
            direct.parent.mkdir()
            direct.write_bytes(b"")
            self.assertEqual(
                direct.resolve(),
                run_ci.find_dumpbin(
                    environment={},
                    which=lambda name: str(direct) if name == "dumpbin" else None,
                ),
            )

            installation = root / "Visual Studio" / "2022" / "Enterprise"
            older = (installation / "VC" / "Tools" / "MSVC" / "14.9.0" /
                     "bin" / "Hostx64" / "x64" / "dumpbin.exe")
            newer = (installation / "VC" / "Tools" / "MSVC" / "14.42.0" /
                     "bin" / "Hostx64" / "x64" / "dumpbin.exe")
            older.parent.mkdir(parents=True)
            newer.parent.mkdir(parents=True)
            older.write_bytes(b"")
            newer.write_bytes(b"")
            self.assertEqual(
                newer.resolve(),
                run_ci.find_dumpbin(
                    environment={"VSINSTALLDIR": str(installation)},
                    which=lambda unused: None,
                ),
            )

    def test_msvc_export_validator_requires_exact_public_api(self):
        rows = "\n".join(
            "%d %X %08X %s" % (ordinal, ordinal - 1, ordinal, symbol)
            for ordinal, symbol in enumerate(
                sorted(run_ci.WIREHAIR_C_API_EXPORTS), start=1)
        )
        args = mock.Mock(
            build_dir=Path("build/ci-msvc-shared"),
            install_dir=Path("build/install-msvc-shared"),
        )
        library = Path("C:/install/bin/wirehair.dll")
        dumpbin = Path("C:/VS/VC/Tools/dumpbin.exe")
        with mock.patch.object(
                run_ci, "find_shared_library", return_value=library), \
                mock.patch.object(run_ci, "find_dumpbin", return_value=dumpbin), \
                mock.patch.object(
                    run_ci, "capture_output", return_value=rows) as capture:
            run_ci.validate_msvc_exports(args)
        capture.assert_called_once_with(
            [str(dumpbin), "/nologo", "/exports", str(library)],
            log_path=args.build_dir / "ci-logs" / "msvc-dumpbin-exports.txt",
        )

        incomplete = rows.replace("wirehair_decode", "unrelated_export", 1)
        with mock.patch.object(
                run_ci, "find_shared_library", return_value=library), \
                mock.patch.object(run_ci, "find_dumpbin", return_value=dumpbin), \
                mock.patch.object(
                    run_ci, "capture_output", return_value=incomplete), \
                self.assertRaisesRegex(RuntimeError, "missing=.*wirehair_decode"):
            run_ci.validate_msvc_exports(args)

        extra = rows + "\n16 F 00000010 accidental_cpp_export\n"
        with mock.patch.object(
                run_ci, "find_shared_library", return_value=library), \
                mock.patch.object(run_ci, "find_dumpbin", return_value=dumpbin), \
                mock.patch.object(
                    run_ci, "capture_output", return_value=extra), \
                self.assertRaisesRegex(RuntimeError, "extra=.*accidental_cpp_export"):
            run_ci.validate_msvc_exports(args)


if __name__ == "__main__":
    unittest.main()
