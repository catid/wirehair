#!/usr/bin/env python3
"""Tests for the source-pinned WH2 hook build matrix."""

import contextlib
import importlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
matrix = importlib.import_module("wh2_hook_build_matrix")


def _run(command, cwd):
    return subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=120,
    )


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _init_repository(path, files):
    path.mkdir()
    for relative, contents in files.items():
        _write(path / relative, contents)
    _run(["git", "init", "-q"], path)
    _run(["git", "config", "user.email", "matrix@example.invalid"], path)
    _run(["git", "config", "user.name", "Matrix Test"], path)
    _run(["git", "add", "--", *sorted(files)], path)
    _run(["git", "commit", "-qm", "fixture"], path)


class SourcePinTests(unittest.TestCase):
    def test_clean_commit_binds_tree_and_worktree_sha256(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "source"
            _init_repository(repository, {
                "main.cpp": "int main() { return 0; }\n",
                "nested/header.h": "#define VALUE 7\n",
            })
            receipt = matrix.source_receipt(repository)
            self.assertRegex(receipt["commit"], r"^[0-9a-f]{40}$")
            self.assertRegex(receipt["tree"], r"^[0-9a-f]{40}$")
            self.assertEqual(receipt["tracked_file_count"], 2)
            self.assertEqual(
                receipt["tracked_manifest_sha256"],
                matrix.canonical_sha256(receipt["tracked_files"]),
            )
            paths = {
                row["path"]: row for row in receipt["tracked_files"]
            }
            self.assertEqual(
                paths["main.cpp"]["sha256"],
                __import__("hashlib").sha256(
                    b"int main() { return 0; }\n").hexdigest(),
            )

    def test_modified_and_untracked_sources_fail_closed(self):
        cases = ("modified", "untracked")
        for case in cases:
            with self.subTest(case=case):
                with tempfile.TemporaryDirectory() as temporary:
                    repository = Path(temporary) / "source"
                    _init_repository(repository, {"main.cpp": "old\n"})
                    if case == "modified":
                        _write(repository / "main.cpp", "new\n")
                    else:
                        _write(repository / "untracked.cpp", "new\n")
                    with self.assertRaisesRegex(
                            matrix.MatrixError, "dirty|differs"):
                        matrix.source_receipt(repository)

    def test_assume_unchanged_file_still_fails_object_check(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "source"
            _init_repository(repository, {"main.cpp": "old\n"})
            _run(["git", "update-index", "--assume-unchanged", "main.cpp"],
                 repository)
            _write(repository / "main.cpp", "new\n")
            with self.assertRaisesRegex(matrix.MatrixError, "differs"):
                matrix.source_receipt(repository)

    def test_ignored_file_and_tracked_symlink_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "ignored"
            _init_repository(repository, {
                ".gitignore": "ignored.h\n",
                "main.cpp": "int main() { return 0; }\n",
            })
            _write(repository / "ignored.h", "#define INJECTED 1\n")
            with self.assertRaisesRegex(matrix.MatrixError, "ignored"):
                matrix.source_receipt(repository)

        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "symlink"
            _init_repository(repository, {"target": "contents\n"})
            os.symlink("target", repository / "link")
            _run(["git", "add", "link"], repository)
            _run(["git", "commit", "-qm", "tracked symlink"], repository)
            with self.assertRaisesRegex(matrix.MatrixError, "unsupported"):
                matrix.source_receipt(repository)


class EvidenceTests(unittest.TestCase):
    def _compile_database(self, root, token):
        source = root / "source"
        build = root / "build"
        source.mkdir()
        build.mkdir()
        _write(source / "probe.cpp", "int main() { return 0; }\n")
        path = build / "compile_commands.json"
        path.write_text(json.dumps([{
            "directory": str(build),
            "command":
                f"/usr/bin/c++ {token} -O3 -c "
                f"{source / 'probe.cpp'} -o CMakeFiles/probe.o",
            "file": str(source / "probe.cpp"),
            "output": "CMakeFiles/probe.o",
        }]), encoding="utf-8")
        return source, build, path

    def test_compile_projection_strips_only_expected_hook_macro(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, build, path = self._compile_database(
                root, f"-D{matrix.HOOK_MACRO}=1")
            evidence = matrix.compile_commands_evidence(
                path, state="on", source=source, build=build)
            self.assertEqual(evidence["entries"], 1)
            tokens = evidence["normalized"][0]["tokens"]
            self.assertNotIn(matrix.HOOK_MACRO, " ".join(tokens))
            self.assertIn("-O3", tokens)

    def test_macro_wrong_state_value_and_response_file_fail_closed(self):
        cases = (
            ("on", "", "exactly once"),
            ("off", f"-D{matrix.HOOK_MACRO}=1", "hook-off"),
            ("on", f"-D{matrix.HOOK_MACRO}=2", "unexpected value"),
            ("on", "@flags.rsp", "response-file"),
            ("off", "-Wl,@link.rsp", "response-file"),
            ("off", "-Wp,@preprocessor.rsp", "response-file"),
            ("off", "-Wa,@assembler.rsp", "response-file"),
        )
        for state, token, error in cases:
            with self.subTest(state=state, token=token):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    source, build, path = self._compile_database(
                        root, token)
                    with self.assertRaisesRegex(matrix.MatrixError, error):
                        matrix.compile_commands_evidence(
                            path, state=state,
                            source=source, build=build)

    def test_reserved_cmake_and_flag_injection_are_rejected(self):
        for value in (
                "BUILD_TESTS=ON",
                "CMAKE_CXX_FLAGS=-O3",
                "WH_PGO_DIR=/external/profile",
                f"EXTRA=-D{matrix.HOOK_MACRO}=1",
                "OTHER=ON"):
            with self.subTest(value=value):
                with self.assertRaises(matrix.MatrixError):
                    matrix._parse_definition(value)
        for value in (f"-D{matrix.HOOK_MACRO}=1", "-O3", "-Wl,-z,now"):
            with self.subTest(value=value):
                with self.assertRaises(matrix.MatrixError):
                    matrix._validate_flag_string(value, "flags")

    def test_build_environment_removes_unbound_external_inputs(self):
        injected = {
            "PATH": os.environ.get("PATH", "/usr/bin"),
            "CPATH": "/external/include",
            "CPLUS_INCLUDE_PATH": "/external/cxx",
            "LIBRARY_PATH": "/external/lib",
            "LD_PRELOAD": "/external/preload.so",
            "GCC_EXEC_PREFIX": "/external/gcc",
            "COMPILER_PATH": "/external/compiler",
            "SOURCE_DATE_EPOCH": "1",
            "CCACHE_DIR": "/external/ccache",
            "CMAKE_TOOLCHAIN_FILE": "/external/toolchain",
            "WIREHAIR_V2_PEEL_DEGREES": "1,1",
        }
        with mock.patch.dict(os.environ, injected, clear=True):
            environment, evidence = matrix._sanitized_environment()
        for name in injected:
            if name == "PATH":
                self.assertEqual(environment[name], injected[name])
            else:
                self.assertNotIn(name, environment)
        self.assertEqual(
            evidence["policy"], "minimal-allowlist-v1")
        self.assertEqual(evidence["removed_count"], len(injected) - 1)
        self.assertEqual(environment["LC_ALL"], "C")
        self.assertEqual(environment["TZ"], "UTC")

    def test_git_subprocess_environment_is_minimal_and_pinned(self):
        captured = {}

        def fake_run(command, **kwargs):
            captured.update(kwargs["environment"])
            return {"stdout": b"", "stderr": b"", "command": command}

        injected = {
            "PATH": "/usr/bin",
            "HOME": "/tmp/home",
            "BASH_ENV": "/secret/bash-env",
            "LD_PRELOAD": "/secret/preload",
            "API_TOKEN": "secret",
            "GIT_DIR": "/forged/repository",
        }
        with mock.patch.dict(os.environ, injected, clear=True):
            with mock.patch.object(matrix, "_run", side_effect=fake_run):
                matrix._git(Path("/tmp"), ["version"])
        self.assertEqual(captured["PATH"], "/usr/bin")
        self.assertEqual(captured["HOME"], "/tmp/home")
        for name in ("BASH_ENV", "LD_PRELOAD", "API_TOKEN", "GIT_DIR"):
            self.assertNotIn(name, captured)
        self.assertEqual(captured["GIT_CONFIG_NOSYSTEM"], "1")
        self.assertEqual(captured["GIT_NO_REPLACE_OBJECTS"], "1")

    def test_protocol_parser_rejects_pgo_profiles(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                matrix._parser().parse_args([
                    "--build-root", "/tmp/build",
                    "--output", "/tmp/output",
                    "--target", "target",
                    "--binary-relative", "target",
                    "--pgo", "GENERATE",
                ])

    def test_cache_rejects_project_include_and_linker_launch_seams(self):
        required = {
            "BUILD_CODEC_V2": "ON",
            "BUILD_TESTS": "OFF",
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_C_FLAGS": "",
            "CMAKE_CXX_FLAGS": "",
            "CMAKE_EXE_LINKER_FLAGS": "",
            "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
            "CMAKE_FIND_USE_PACKAGE_REGISTRY": "OFF",
            "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY": "OFF",
            "CMAKE_MODULE_LINKER_FLAGS": "",
            "CMAKE_SHARED_LINKER_FLAGS": "",
            "WH_PGO_MODE": "OFF",
            "WIREHAIR_BUILD_BENCHMARKS": "ON",
        }
        seams = (
            "CMAKE_PROJECT_TOP_LEVEL_INCLUDES",
            "CMAKE_PROJECT_matrix_probe_INCLUDE",
            "CMAKE_PROJECT_matrix_probe_INCLUDE_BEFORE",
            "CMAKE_C_LINKER_LAUNCHER",
            "CMAKE_CXX_LINKER_LAUNCHER",
        )
        for seam in seams:
            with self.subTest(seam=seam):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    source = root / "source"
                    build = root / "build"
                    source.mkdir()
                    build.mkdir()
                    cache = dict(required)
                    cache[seam] = "/external/injector"
                    cache_path = build / "CMakeCache.txt"
                    _write(cache_path, "".join(
                        f"{name}:STRING={value}\n"
                        for name, value in sorted(cache.items())
                    ))
                    with self.assertRaisesRegex(
                            matrix.MatrixError, seam):
                        matrix.cache_evidence(
                            cache_path, state="off",
                            source=source, build=build)

    def test_receipt_publication_is_atomic_and_no_replace(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "receipt.json"
            matrix._atomic_write_receipt(path, b"first\n")
            self.assertEqual(path.read_bytes(), b"first\n")
            with self.assertRaises(matrix.MatrixError):
                matrix._atomic_write_receipt(path, b"second\n")
            self.assertEqual(path.read_bytes(), b"first\n")
            self.assertEqual(
                list(Path(temporary).glob(".*.tmp")), [])

    def test_tool_version_receipt_rejects_executable_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            tool = Path(temporary) / "tool"
            _write(tool, "#!/bin/sh\nexit 0\n")
            tool.chmod(0o755)

            def mutate_tool(unused_command, **unused_kwargs):
                _write(tool, "#!/bin/sh\nexit 1\n")
                tool.chmod(0o755)
                return {
                    "command": [str(tool), "--version"],
                    "stdout": b"version\n",
                    "stderr": b"",
                }

            with mock.patch.object(
                    matrix, "_run", side_effect=mutate_tool):
                with self.assertRaisesRegex(
                        matrix.MatrixError, "changed"):
                    matrix._tool_receipt(
                        tool,
                        cwd=Path(temporary),
                        environment={"PATH": temporary},
                        timeout=10,
                    )

    def test_tool_recheck_rejects_receipt_change(self):
        expected = {"sha256": "old"}
        with mock.patch.object(
                matrix, "_tool_receipt",
                return_value={"sha256": "new"}):
            with self.assertRaisesRegex(
                    matrix.MatrixError, "Python changed"):
                matrix._require_tool_unchanged(
                    "Python", "/python", expected,
                    cwd=Path("/tmp"),
                    environment={"PATH": "/usr/bin"},
                    timeout=10,
                )


class ComparabilityTests(unittest.TestCase):
    def test_build_order_is_interleaved_by_macro_state(self):
        self.assertEqual(matrix.MATRIX, (
            ("hook-on-a", "on", "a"),
            ("hook-off-a", "off", "a"),
            ("hook-on-b", "on", "b"),
            ("hook-off-b", "off", "b"),
        ))

    def _builds(self):
        builds = {}
        for name, state, replicate in matrix.MATRIX:
            state_digest = "a" * 64 if state == "on" else "b" * 64
            builds[name] = {
                "name": name,
                "cache": {"normalized_sha256": "1" * 64},
                "compile_commands": {"normalized_sha256": "2" * 64},
                "link": {"normalized_sha256": "3" * 64},
                "tools": {"compiler": {"sha256": "4" * 64}},
                "binary": {
                    "bytes": 10,
                    "sha256": state_digest,
                    "build_id": {
                        "value": "aa" if state == "on" else "bb",
                    },
                },
                "hook_state": state,
                "replicate": replicate,
            }
        return builds

    def test_exact_matrix_is_explicitly_comparable(self):
        result = matrix.validate_comparability(self._builds())
        self.assertTrue(result["cache_normalized_sha256_identical"])
        self.assertTrue(result["compile_commands_normalized_sha256_identical"])
        self.assertTrue(result["link_normalized_sha256_identical"])
        self.assertTrue(result["hook_on_ab_binary_identical"])
        self.assertTrue(result["hook_off_ab_binary_identical"])
        self.assertFalse(result["hook_on_off_binary_identical"])

    def test_flag_or_ab_binary_difference_fails_closed(self):
        cases = (
            lambda builds:
                builds["hook-off-b"]["compile_commands"].__setitem__(
                    "normalized_sha256", "9" * 64),
            lambda builds:
                builds["hook-on-b"]["binary"].__setitem__(
                    "sha256", "9" * 64),
        )
        for mutate in cases:
            with self.subTest(mutate=mutate):
                builds = self._builds()
                mutate(builds)
                with self.assertRaises(matrix.MatrixError):
                    matrix.validate_comparability(builds)

    def test_identical_hook_on_off_binaries_fail_closed(self):
        builds = self._builds()
        on_binary = builds["hook-on-a"]["binary"]
        for name in ("hook-off-a", "hook-off-b"):
            builds[name]["binary"]["bytes"] = on_binary["bytes"]
            builds[name]["binary"]["sha256"] = on_binary["sha256"]
            builds[name]["binary"]["build_id"]["value"] = \
                on_binary["build_id"]["value"]
        with self.assertRaisesRegex(matrix.MatrixError, "identical"):
            matrix.validate_comparability(builds)


@unittest.skipUnless(
    shutil.which("cmake") and shutil.which("readelf") and
    shutil.which("c++"),
    "requires the native CMake/ELF toolchain",
)
class TinyEndToEndTests(unittest.TestCase):
    def test_four_fresh_release_trees_publish_authenticated_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "source"
            helper = Path(matrix.__file__).read_text(encoding="utf-8")
            files = {
                "bench/wh2_hook_build_matrix.py": helper,
                "CMakeLists.txt": (
                    "cmake_minimum_required(VERSION 3.15)\n"
                    "project(matrix_probe LANGUAGES C CXX)\n"
                    "add_executable(probe main.cpp helper.c)\n"
                ),
                "main.cpp": (
                    "extern \"C\" int helper(void);\n"
                    "int main() {\n"
                    "  const int value = helper();\n"
                    "  return (value == 1 || value == 2) ? 0 : 1;\n"
                    "}\n"
                ),
                "helper.c": (
                    "#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)\n"
                    "int helper(void) { return 1; }\n"
                    "#else\n"
                    "int helper(void) { return 2; }\n"
                    "#endif\n"
                ),
            }
            _init_repository(repository, files)
            build_root = root / "builds"
            output = root / "receipt.json"
            command = [
                sys.executable,
                str(repository / "bench/wh2_hook_build_matrix.py"),
                "--source", str(repository),
                "--build-root", str(build_root),
                "--output", str(output),
                "--target", "probe",
                "--binary-relative", "probe",
                "--jobs", "2",
                "--timeout", "120",
            ]
            result = _run(command, root)
            announcement = json.loads(result.stdout)
            receipt = json.loads(output.read_bytes())
            self.assertEqual(receipt["schema"], matrix.SCHEMA)
            self.assertEqual(
                announcement["sha256"],
                __import__("hashlib").sha256(output.read_bytes()).hexdigest(),
            )
            self.assertEqual(set(receipt["builds"]), {
                name for name, unused_state, unused_replicate
                in matrix.MATRIX
            })
            self.assertTrue(
                receipt["comparability"]["hook_on_ab_binary_identical"])
            self.assertTrue(
                receipt["comparability"]["hook_off_ab_binary_identical"])
            self.assertFalse(
                receipt["comparability"]["hook_on_off_binary_identical"])
            self.assertEqual(
                receipt["configuration"]["build_order"],
                [name for name, unused_state, unused_replicate
                 in matrix.MATRIX],
            )
            self.assertFalse(receipt["timing_executed"])
            python = receipt["tools"]["python"]
            self.assertEqual(
                Path(python["configured_path"]),
                Path(os.path.abspath(sys.executable)),
            )
            self.assertEqual(
                Path(python["resolved_path"]),
                Path(sys.executable).resolve(),
            )
            self.assertRegex(python["sha256"], r"^[0-9a-f]{64}$")
            self.assertIn("version", python)
            for name, state, unused_replicate in matrix.MATRIX:
                build = receipt["builds"][name]
                self.assertEqual(build["hook_state"], state)
                self.assertRegex(
                    build["binary"]["build_id"]["value"],
                    r"^[0-9a-f]+$",
                )
                self.assertTrue(
                    (build_root / name / "compile_commands.json").is_file())
                for tool in (
                        "CMAKE_AR", "CMAKE_CXX_COMPILER", "CMAKE_LINKER",
                        "CMAKE_MAKE_PROGRAM", "CMAKE_NM", "CMAKE_OBJDUMP",
                        "CMAKE_RANLIB"):
                    self.assertIn(tool, build["tools"])


if __name__ == "__main__":
    unittest.main()
