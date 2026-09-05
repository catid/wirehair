#!/usr/bin/env python3
"""Synthetic tests for the public borrowed-source r1 supervisor/reducer."""

import contextlib
import copy
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock


sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2PublicBorrowedCurrentVsWh1R1 as screen  # noqa: E402


def health_receipt(non_idle_ticks=10):
    ticks = [non_idle_ticks, 0, 0, 100, 0, 0, 0, 0]
    return {
        "available_cpu_count": 1,
        "controller_affinity": [screen.CONTROLLER_CPU],
        "controller_cpu": screen.CONTROLLER_CPU,
        "edac": [{"counter": "ce_count",
                  "path": "/sys/devices/system/edac/mc/mc0/ce_count", "value": 0}],
        "governor": "powersave",
        "loadavg": "0.00 0.00 0.00 1/100 123",
        "sibling_cpu": screen.TARGET_SIBLING,
        "sibling_non_idle_ticks": non_idle_ticks,
        "sibling_proc_stat": "cpu{} {}".format(
            screen.TARGET_SIBLING, " ".join(str(value) for value in ticks)),
        "sibling_tick_fields": ticks,
        "target_cpu": screen.TARGET_CPU,
        "thermal": [{"path": "/sys/class/hwmon/hwmon0/temp1_input",
                     "value_millic": 45000}],
        "thermal_endpoint_max_millic": screen.THERMAL_ENDPOINT_MAX_MILLIC,
        "thread_siblings_list": "56,120",
    }


def artifact(path, digest=None):
    return {
        "device": 1, "inode": int(screen.synthetic_hash(path)[:8], 16),
        "mode": 0o755, "mtime_ns": 1,
        "path": path, "sha256": digest or screen.synthetic_hash(path),
        "size": 128,
    }


def text_hash(text):
    return screen.sha256_bytes(text.encode("utf-8"))


def command(tokens):
    return " ".join(token if token == "&&" else shlex.quote(token)
                    for token in tokens)


def build_fixture(expected, interpreter, source_root, build_root, worker):
    compiler = expected["compiler_path"]
    definitions = [
        "-DWIREHAIR_DLL=1",
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_ID="GNU"',
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION="{}"'.format(
            expected["compiler_version"]),
        '-DWIREHAIR_WH2_CXX_COMPILER_PATH="{}"'.format(compiler),
        '-DWIREHAIR_WH2_CXX_COMPILER_SHA256="{}"'.format(
            expected["compiler_sha256"]),
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(expected["source_commit"]),
    ]
    library_definitions = [
        "-DWIREHAIR_BUILDING=1", "-DWIREHAIR_DLL=1", "-Dwirehair_EXPORTS"]

    def compile_rows(sources, prefix, macros, pic):
        rows = {}
        objects = []
        for relative in sources:
            obj = prefix + relative + ".o"
            objects.append(obj)
            rows[relative] = command(
                [compiler] + macros + ["-I" + source_root + "/include",
                "-O3", "-DNDEBUG", "-std=gnu++11"] + pic +
                ["-Wall", "-Wextra", "-march=native", "-MD", "-MT", obj,
                 "-MF", obj + ".d", "-o", obj, "-c", source_root + "/" + relative])
        return rows, objects

    compile_commands, objects = compile_rows(
        screen.R1_TARGET_SOURCES,
        "CMakeFiles/wirehair_wh2_public_borrowed_current_vs_wh1_r1.dir/",
        definitions, [])
    library_commands, library_objects = compile_rows(
        screen.R1_LIBRARY_SOURCES, "CMakeFiles/wirehair.dir/",
        library_definitions, ["-fPIC"])
    library_link = command([
        ":", "&&", compiler, "-fPIC", "-O3", "-DNDEBUG",
        "-Wl,--version-script=" + source_root + "/abi/wirehair.map",
        "-shared", "-Wl,-soname,libwirehair.so.2", "-o", "libwirehair.so.2.0.0",
    ] + library_objects + ["-lm", "&&", ":"])
    link = command([
        ":", "&&", compiler, "-O3", "-DNDEBUG",
    ] + objects + [
        "-o", Path(worker).name, "-Wl,-rpath," + build_root,
        "libwirehair.so.2.0.0", "-ldl", "&&", ":",
    ])
    commands = "\n".join(list(library_commands.values()) + [library_link] +
                         list(compile_commands.values()) + [link]) + "\n"
    dry_run = "ninja: no work to do.\n"
    return {
        "cache": {"BUILD_SHARED_LIBS": "ON", "CMAKE_BUILD_TYPE": "Release",
                  "CMAKE_GENERATOR": "Ninja", "MARCH_NATIVE": "ON",
                  "WH_LTO": "OFF", "WH_PGO_MODE": "OFF"},
        "cmake_build_root": build_root, "cmake_compiler": compiler,
        "cmake_compiler_resolved": compiler,
        "cmake_interpreter": interpreter["path"],
        "cmake_interpreter_resolved": interpreter["path"],
        "cmake_source_root": source_root, "compiler_invocation": compiler,
        "compile_command_sha256": {name: text_hash(row)
                                   for name, row in compile_commands.items()},
        "compile_commands": compile_commands, "compile_definitions": definitions,
        "dry_run_sha256": text_hash(dry_run), "dry_run_text": dry_run,
        "library_compile_command_sha256": {name: text_hash(row)
                                           for name, row in library_commands.items()},
        "library_compile_commands": library_commands,
        "library_compile_definitions": library_definitions,
        "library_link_command": library_link,
        "library_link_command_sha256": text_hash(library_link),
        "library_link_object_roster": library_objects,
        "library_sources": list(screen.R1_LIBRARY_SOURCES),
        "link_command": link, "link_command_sha256": text_hash(link),
        "link_object_roster": objects, "ninja_commands_sha256": text_hash(commands),
        "ninja_commands_text": commands,
        "target_sources": list(screen.R1_TARGET_SOURCES),
    }


def dynamic_fixture(worker, library, library_artifact):
    symbols = sorted(screen.REQUIRED_DYNAMIC_SYMBOLS)
    imports = [name + "@" + screen.WIREHAIR_SYMBOL_VERSION for name in symbols]
    library_nm = "".join(name + "@@WIREHAIR_2.0 T 10 10\n" for name in symbols)
    worker_nm = "".join(name + " U\n" for name in imports)
    full_nm = "main T 10 10\nwirehair::wh2_benchmark::Run() T 20 10\n"
    ldd = "libwirehair.so.2 => {} (0xADDR)\n".format(library)
    readelf = " 0x0000000000000001 (NEEDED) Shared library: [libwirehair.so.2]\n"
    nm_version = "GNU nm (GNU Binutils) 2.42\n"
    return {
        "forbidden_symbol_fragments": list(screen.FORBIDDEN_SYMBOL_FRAGMENTS),
        "ldd_library_path": library, "ldd_library_resolved_path": library,
        "ldd_library_stat": {key: value for key, value in library_artifact.items()
                             if key != "sha256"},
        "ldd_normalized_sha256": text_hash(ldd), "ldd_normalized_text": ldd,
        "library_nm_argv": ["/usr/bin/nm", "-D", "--format=posix",
                            "--defined-only", library],
        "library_symbol_table_text": library_nm,
        "nm_artifact": artifact("/usr/bin/x86_64-linux-gnu-nm"),
        "nm_invoked_path": "/usr/bin/nm",
        "nm_version_sha256": text_hash(nm_version), "nm_version_text": nm_version,
        "readelf_dynamic_sha256": text_hash(readelf),
        "readelf_dynamic_text": readelf,
        "required_dynamic_symbols": symbols, "required_dynamic_symbols_verified": True,
        "symbol_table_sha256": text_hash(library_nm),
        "worker_forbidden_symbol_fragments": list(screen.WORKER_FORBIDDEN_SYMBOL_FRAGMENTS),
        "worker_full_nm_argv": ["/usr/bin/nm", "-C", "--format=posix", worker],
        "worker_full_symbol_table_sha256": text_hash(full_nm),
        "worker_full_symbol_table_text": full_nm,
        "worker_import_symbol_version": screen.WIREHAIR_SYMBOL_VERSION,
        "worker_import_symbols": imports, "worker_import_symbols_verified": True,
        "worker_nm_argv": ["/usr/bin/nm", "-D", "--format=posix", worker],
        "worker_symbol_table_sha256": text_hash(worker_nm),
        "worker_symbol_table_text": worker_nm,
    }


def bundle_fixture(invalid_launch=False, **ratios):
    """Portable synthetic evidence: no host probes and no benchmark invocation."""
    raw, expected, maps_hash = screen.synthetic_transcript(**ratios)
    config, _, statistics, invalid, reject = screen.parse_transcript(
        raw, expected, maps_hash)
    expected.update({
        "native_config_identity_sha256": config["native_config_identity_sha256"],
        "target_identity_sha256": config["target_identity"]["canonical_sha256"],
    })
    source_root = "/tmp/r1-source"
    build_root = "/tmp/r1"
    worker = build_root + "/wirehair_wh2_public_borrowed_current_vs_wh1_r1"
    controller = source_root + "/bench/Wh2PublicBorrowedCurrentVsWh1R1.py"
    library = expected["library_path"]
    compiler_version = "g++ (GCC) {}\n".format(expected["compiler_version"])
    compiler = {
        "path": expected["compiler_path"], "sha256": expected["compiler_sha256"],
        "version": expected["compiler_version"], "version_text": compiler_version,
        "version_text_sha256": text_hash(compiler_version),
        "version_sha256": text_hash(expected["compiler_version"] + "\n"),
    }
    interpreter_path = "/usr/bin/python3.8"
    interpreter = {
        "artifact": artifact(interpreter_path), "dont_write_bytecode": True,
        "ignore_environment": True, "implementation": "cpython",
        "invoked_path": interpreter_path, "isolated": True, "no_site": True,
        "no_user_site": True, "path": interpreter_path,
        "version": "3.8.20 (synthetic) [GCC 13.3.0]",
        "version_info": [3, 8, 20, "final", 0],
    }
    artifacts = {
        "compiler": artifact(compiler["path"], compiler["sha256"]),
        "controller": artifact(controller),
        "controller_interpreter": interpreter["artifact"],
        "library": artifact(library), "worker": artifact(worker),
    }
    git_version = "git version 2.43.0\n"
    git = {
        "commit": expected["source_commit"], "git_artifact": artifact("/usr/bin/git"),
        "git_invoked_path": "/usr/bin/git", "git_version_sha256": text_hash(git_version),
        "git_version_text": git_version,
        "required_tracked_inputs": {path: screen.synthetic_hash(path)[:40]
                                   for path in screen.R1_TRACKED_INPUTS},
        "tracked_clean": True, "tree": "3" * 40,
        "upstream_commit": expected["source_commit"],
    }
    build = build_fixture(expected, interpreter, source_root, build_root, worker)
    dynamic = dynamic_fixture(worker, library, artifacts["library"])
    r0_entries = [{"exists": False, "path": str(path)} for path in screen.R0_ROOTS]
    r0 = {"entries": r0_entries,
          "snapshot_sha256": screen.sha256_bytes(screen.canonical_bytes(r0_entries))}
    before, after = health_receipt(), health_receipt(12)
    native = {
        "config_sha256": screen.sha256_bytes(screen.canonical_bytes(config) + b"\n"),
        "native_config_identity_sha256": expected["native_config_identity_sha256"],
        "runtime_library_maps_sha256": maps_hash, "runtime_library_path": library,
        "stderr_sha256": screen.sha256_bytes(b""),
        "target_identity_sha256": expected["target_identity_sha256"],
    }
    claim = {
        "campaign": screen.CAMPAIGN,
        "controller_interpreter_sha256": interpreter["artifact"]["sha256"],
        "controller_sha256": artifacts["controller"]["sha256"], "created_unix_ns": 1,
        "library_sha256": artifacts["library"]["sha256"],
        "native_config_identity_sha256": expected["native_config_identity_sha256"],
        "schema": screen.CLAIM_SCHEMA, "source_commit": expected["source_commit"],
        "worker_sha256": artifacts["worker"]["sha256"],
    }
    claim_hash = screen.sha256_bytes(screen.canonical_bytes(claim) + b"\n")
    marker = {
        "campaign": screen.CAMPAIGN, "claim_sha256": claim_hash,
        "native_config_identity_sha256": expected["native_config_identity_sha256"],
        "schema": screen.WORKER_STARTED_SCHEMA,
        "source_commit": expected["source_commit"], "status": "started",
        "worker_sha256": claim["worker_sha256"],
    }
    provenance = {
        "artifacts": artifacts, "artifacts_after": copy.deepcopy(artifacts),
        "build": build, "build_after": copy.deepcopy(build), "build_root": build_root,
        "campaign": screen.CAMPAIGN, "claim_sha256": claim_hash, "compiler": compiler,
        "controller_affinity": {
            "affinity_before": [0, 56, 120], "affinity_after": [0], "controller_cpu": 0,
            "sibling_available_before": True, "singleton_verified": True,
            "target_available_before": True,
        },
        "controller_path": controller, "dynamic": dynamic,
        "dynamic_after": copy.deepcopy(dynamic), "expected": expected,
        "git": git, "git_after": copy.deepcopy(git),
        "health_adjudication": screen.validate_health_pair(before, after),
        "health_after": after, "health_before": before,
        "interpreter": interpreter, "interpreter_after": copy.deepcopy(interpreter),
        "library_path": library, "native_config": native,
        "outer_deadline_seconds": int(screen.OUTER_DEADLINE_SECONDS),
        "r0_after": copy.deepcopy(r0), "r0_before": r0,
        "schema": screen.PROVENANCE_SCHEMA, "source_root": source_root,
        "worker_exit": 0, "worker_path": worker, "worker_timed_out": False,
    }
    infrastructure = []
    if invalid_launch:
        raw = b""
        marker = None
        statistics, invalid, reject = None, [], []
        infrastructure = ["execution:synthetic launch failed", "worker_started_missing"]
        provenance["worker_exit"] = None
    raw_hash = screen.sha256_bytes(raw)
    outcome = "invalid" if infrastructure or invalid else "reject" if reject else "pass"
    summary = screen.make_summary(
        outcome, statistics, infrastructure, invalid, reject, raw_hash, 0.125)
    return {"CLAIM": claim, "raw.jsonl": raw, "worker.stderr": b"",
            "summary.json": summary, "provenance.json": provenance,
            "WORKER_STARTED": marker}


def publish_bundle(root, fixture):
    """Refresh transport hashes, leaving semantic cross-links for replay to check."""
    values = copy.deepcopy(fixture)
    marker = values.pop("WORKER_STARTED")
    files = {name: value if type(value) is bytes else
             screen.canonical_bytes(value) + b"\n" for name, value in values.items()}
    if marker is not None:
        files["WORKER_STARTED"] = screen.canonical_bytes(marker) + b"\n"
    provenance = values["provenance.json"]
    provenance.update({
        "claim_sha256": screen.sha256_bytes(files["CLAIM"]),
        "raw_sha256": screen.sha256_bytes(files["raw.jsonl"]),
        "stderr_sha256": screen.sha256_bytes(files["worker.stderr"]),
        "worker_started_sha256": (screen.sha256_bytes(files["WORKER_STARTED"])
                                  if marker is not None else None),
    })
    files["provenance.json"] = screen.canonical_bytes(provenance) + b"\n"
    files["COMPLETE"] = screen.canonical_bytes({
        "campaign": screen.CAMPAIGN, "claim_sha256": screen.sha256_bytes(files["CLAIM"]),
        "outcome": values["summary.json"]["outcome"],
        "provenance_sha256": screen.sha256_bytes(files["provenance.json"]),
        "raw_sha256": screen.sha256_bytes(files["raw.jsonl"]),
        "schema": screen.COMPLETE_SCHEMA, "status": "complete",
        "summary_sha256": screen.sha256_bytes(files["summary.json"]),
    }) + b"\n"
    os.chmod(str(root), 0o700)
    for name, data in files.items():
        path = root / name
        if path.exists():
            os.chmod(str(path), 0o600)
        path.write_bytes(data)
        os.chmod(str(path), 0o400)


class PublicBorrowedR1Test(unittest.TestCase):
    def parse(self, raw, expected, maps_hash):
        return screen.parse_transcript(raw, expected, maps_hash)

    def test_exact_passing_roster_and_reduction(self):
        raw, expected, maps_hash = screen.synthetic_transcript()
        self.assertEqual(len(raw.splitlines()), 1442)
        config, terminal, statistics, invalid, reject = self.parse(
            raw, expected, maps_hash)
        self.assertEqual(invalid, [])
        self.assertEqual(reject, [])
        self.assertEqual(terminal["panel_count"], 1440)
        self.assertEqual(
            config["expected_measured_invocations"], 852480)
        result = statistics["aggregates"]["C/L"][
            "prebuilt-systematic"]["equal-cell"]
        self.assertAlmostEqual(result["geometric_mean"], 0.98, places=12)
        self.assertLess(result["upper95"], 0.99)

    def test_noncanonical_and_roster_mutations_are_invalid(self):
        with self.assertRaises(screen.ValidationError):
            screen.parse_canonical_line(b'{"a": 1}\n', "spacing")
        raw, expected, maps_hash = screen.synthetic_transcript()
        records = [json.loads(line) for line in raw.splitlines()]
        records[1]["sequence"] = 3
        altered = b"".join(
            screen.canonical_bytes(record) + b"\n" for record in records)
        with self.assertRaises(screen.ValidationError):
            self.parse(altered, expected, maps_hash)
        records = [json.loads(line) for line in raw.splitlines()]
        records[17]["runtime_library_maps_sha256"] = "f" * 64
        altered = b"".join(
            screen.canonical_bytes(record) + b"\n" for record in records)
        with self.assertRaises(screen.ValidationError):
            self.parse(altered, expected, maps_hash)
        records = [json.loads(line) for line in raw.splitlines()]
        records[1]["slots"][4]["invocation_count"] += 1
        altered = b"".join(
            screen.canonical_bytes(record) + b"\n" for record in records)
        with self.assertRaises(screen.ValidationError):
            self.parse(altered, expected, maps_hash)

    def test_partial_final_and_maps_receipts_are_mandatory(self):
        raw, expected, maps_hash = screen.synthetic_transcript()
        records = [json.loads(line) for line in raw.splitlines()]
        records[0]["cells"][0]["partial_final_oracles"][0][
            "arms"][0]["source_immutable_verified"] = False
        altered = b"".join(
            screen.canonical_bytes(record) + b"\n" for record in records)
        with self.assertRaises(screen.ValidationError):
            self.parse(altered, expected, maps_hash)

    def test_exact_source_profile_and_partial_cross_arm_contracts(self):
        raw, expected, maps_hash = screen.synthetic_transcript()
        for mutate in (
            lambda config: config.__setitem__(
                "source_seed_base", "0x0000000000000000"),
            lambda config: config["cells"][0].__setitem__(
                "source_seed", "0x0000000000000000"),
            lambda config: config["cells"][0]["partial_final_oracles"][0][
                "arms"][1].__setitem__("systematic_sha256", "f" * 64),
            lambda config: config["cells"][0]["oracles"][0].__setitem__(
                "profile_hex", " " +
                config["cells"][0]["oracles"][0]["profile_hex"][1:]),
        ):
            records = [json.loads(line) for line in raw.splitlines()]
            mutate(records[0])
            altered = b"".join(
                screen.canonical_bytes(record) + b"\n" for record in records)
            with self.assertRaises(screen.ValidationError):
                self.parse(altered, expected, maps_hash)

    def test_health_pair_enforces_sibling_and_edac_caps(self):
        before = health_receipt(10)
        after = health_receipt(15)
        self.assertEqual(
            screen.validate_health_pair(before, after)[
                "sibling_non_idle_tick_delta"], 5)
        after = health_receipt(16)
        with self.assertRaises(screen.ValidationError):
            screen.validate_health_pair(before, after)
        after = health_receipt(15)
        after["edac"][0]["value"] = 1
        with self.assertRaises(screen.ValidationError):
            screen.validate_health_pair(before, after)

    def test_health_receipts_bind_raw_ticks_topology_and_thermal_endpoints(self):
        mutations = {
            "empty thermal": lambda receipt: receipt.__setitem__("thermal", []),
            "too hot": lambda receipt: receipt["thermal"][0].__setitem__(
                "value_millic", screen.THERMAL_ENDPOINT_MAX_MILLIC + 1),
            "duplicate sensor": lambda receipt: receipt["thermal"].append(
                copy.deepcopy(receipt["thermal"][0])),
            "forged tick sum": lambda receipt: receipt.__setitem__(
                "sibling_non_idle_ticks", 9),
            "forged raw stat": lambda receipt: receipt.__setitem__(
                "sibling_proc_stat", "cpu56 9 0 0 100 0 0 0 0"),
            "bad topology": lambda receipt: receipt.__setitem__(
                "thread_siblings_list", "56,not-a-cpu"),
            "wrong controller": lambda receipt: receipt.__setitem__(
                "controller_affinity", [1]),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                receipt = health_receipt()
                mutate(receipt)
                with self.assertRaises(screen.ValidationError):
                    screen.validate_health_receipt(receipt, "synthetic health")

    def test_roundtrip_requires_repair_packets_and_source_hash_closure(self):
        raw, expected, maps_hash = screen.synthetic_transcript()
        mutations = {
            "full source closure": lambda config: config["cells"][0]["oracles"][0].__setitem__(
                "roundtrip_sha256", "f" * 64),
            "full systematic-only": lambda config: config["cells"][0]["oracles"][0].__setitem__(
                "roundtrip_first_id", 0),
            "full repair receipt": lambda config: config["cells"][0]["oracles"][0].__setitem__(
                "roundtrip_repair_only_verified", False),
            "full too few packets": lambda config: config["cells"][0]["oracles"][0].__setitem__(
                "roundtrip_packet_count", 1),
            "partial source closure": lambda config: config["cells"][0]["partial_final_oracles"][0][
                "arms"][0].__setitem__("roundtrip_sha256", "f" * 64),
            "partial systematic-only": lambda config: config["cells"][0]["partial_final_oracles"][0][
                "arms"][0].__setitem__("roundtrip_first_id", 0),
            "partial too many packets": lambda config: config["cells"][0]["partial_final_oracles"][0][
                "arms"][0].__setitem__("roundtrip_packet_count", 73),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                records = [json.loads(line) for line in raw.splitlines()]
                mutate(records[0])
                records[0]["native_config_identity_sha256"] = (
                    screen.native_config_identity_sha256(records[0]))
                altered = b"".join(screen.canonical_bytes(record) + b"\n"
                                   for record in records)
                with self.assertRaises(screen.ValidationError):
                    self.parse(altered, expected, maps_hash)

    def test_aa_is_invalid_and_speed_miss_is_reject(self):
        raw, expected, maps_hash = screen.synthetic_transcript(aa_ratio=1.03)
        _, _, _, invalid, reject = self.parse(raw, expected, maps_hash)
        self.assertTrue(invalid)
        self.assertTrue(all(item.startswith("aa_ci:") for item in invalid))
        raw, expected, maps_hash = screen.synthetic_transcript(
            cl_systematic_ratio=1.0)
        _, _, _, invalid, reject = self.parse(raw, expected, maps_hash)
        self.assertEqual(invalid, [])
        self.assertTrue(reject)
        self.assertTrue(any(item.startswith("systematic_group_upper95:")
                            for item in reject))

    def test_controller_selftest_does_not_touch_fixed_namespace(self):
        before = screen.fixed_namespace_snapshot()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(screen.selftest(), 0)
        self.assertEqual(screen.fixed_namespace_snapshot(), before)

    def test_selftest_cli_rejects_workload_arguments(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            screen.parse_args(["--selftest", "--worker", "/tmp/nope"])

    def replay_fixture(self, fixture):
        with tempfile.TemporaryDirectory(prefix="wh2-r1-replay-") as temporary:
            root = Path(temporary).resolve()
            publish_bundle(root, fixture)
            with contextlib.redirect_stdout(io.StringIO()):
                return screen.replay(SimpleNamespace(bundle=str(root)))

    def test_complete_synthetic_bundles_replay_each_disposition(self):
        before = screen.fixed_namespace_snapshot()
        self.assertEqual(self.replay_fixture(bundle_fixture()), 0)
        self.assertEqual(self.replay_fixture(bundle_fixture(
            cl_systematic_ratio=1.0)), 2)
        self.assertEqual(self.replay_fixture(bundle_fixture(aa_ratio=1.03)), 1)
        self.assertEqual(self.replay_fixture(bundle_fixture(invalid_launch=True)), 1)
        self.assertEqual(screen.fixed_namespace_snapshot(), before)

    def test_replay_preserves_invalid_postflight_prefix(self):
        fixture = bundle_fixture(invalid_launch=True)
        fixture["summary.json"]["infrastructure_failures"].append(
            "postflight:synthetic health probe failure")
        for field in screen.PROVENANCE_POSTFLIGHT_FIELDS:
            fixture["provenance.json"][field] = None
        self.assertEqual(self.replay_fixture(fixture), 1)
        fixture["provenance.json"]["artifacts_after"] = copy.deepcopy(
            fixture["provenance.json"]["artifacts"])
        with self.assertRaisesRegex(screen.ValidationError, "not a prefix"):
            self.replay_fixture(fixture)

    def test_native_config_replay_allows_only_maps_to_differ(self):
        fixture = bundle_fixture()
        config = json.loads(fixture["raw.jsonl"].splitlines()[0])
        preflight_maps = screen.synthetic_hash("different-preflight-maps")
        config["runtime_library_maps_sha256"] = preflight_maps
        native = fixture["provenance.json"]["native_config"]
        native["runtime_library_maps_sha256"] = preflight_maps
        native["config_sha256"] = screen.sha256_bytes(
            screen.canonical_bytes(config) + b"\n")
        self.assertEqual(self.replay_fixture(fixture), 0)
        native["config_sha256"] = "f" * 64
        with self.assertRaises(screen.ValidationError):
            self.replay_fixture(fixture)

    def test_invalid_bundle_replay_verifies_exact_file_hashes(self):
        fixture = bundle_fixture(invalid_launch=True)
        with tempfile.TemporaryDirectory(prefix="wh2-r1-replay-") as temporary:
            root = Path(temporary).resolve()
            publish_bundle(root, fixture)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(screen.replay(SimpleNamespace(bundle=str(root))), 1)
            raw_path = root / "raw.jsonl"
            os.chmod(str(raw_path), 0o600)
            raw_path.write_bytes(b"\n")
            os.chmod(str(raw_path), 0o400)
            with self.assertRaises(screen.ValidationError):
                screen.replay(SimpleNamespace(bundle=str(root)))

    def test_replay_rejects_rehashed_semantic_cross_link_tampering(self):
        fixture = bundle_fixture()
        mutations = {
            "missing postflight": lambda value: value["provenance.json"].__setitem__(
                "build_after", None),
            "forged health adjudication": lambda value: value["provenance.json"][
                "health_adjudication"].__setitem__("sibling_non_idle_tick_delta", 0),
            "wrong upstream": lambda value: value["provenance.json"]["git"].__setitem__(
                "upstream_commit", "f" * 40),
            "wrong tracked blob roster": lambda value: value["provenance.json"]["git"][
                "required_tracked_inputs"].pop(screen.R1_TRACKED_INPUTS[0]),
            "wrong claim commit": lambda value: value["CLAIM"].__setitem__(
                "source_commit", "f" * 40),
            "wrong worker artifact": lambda value: value["provenance.json"]["artifacts"][
                "worker"].__setitem__("sha256", "f" * 64),
            "wrong interpreter": lambda value: value["provenance.json"]["interpreter"].__setitem__(
                "no_site", False),
            "wrong native identity": lambda value: value["provenance.json"]["native_config"].__setitem__(
                "native_config_identity_sha256", "f" * 64),
            "wrong marker worker": lambda value: value["WORKER_STARTED"].__setitem__(
                "worker_sha256", "f" * 64),
            "missing marker": lambda value: value.__setitem__("WORKER_STARTED", None),
            "wrong exit": lambda value: value["provenance.json"].__setitem__("worker_exit", 1),
            "nonempty stderr": lambda value: value.__setitem__("worker.stderr", b"unexpected\n"),
            "wrong statistics": lambda value: value["summary.json"]["statistics"]["aggregates"][
                "C/L"]["prebuilt-systematic"]["equal-cell"].__setitem__("geometric_mean", 0.5),
            "wrong R0 hash": lambda value: value["provenance.json"]["r0_before"].__setitem__(
                "snapshot_sha256", "f" * 64),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                altered = copy.deepcopy(fixture)
                mutate(altered)
                with self.assertRaises(screen.ValidationError):
                    self.replay_fixture(altered)

    def test_replay_reparses_rehashed_compile_and_import_receipts(self):
        fixture = bundle_fixture()
        for kind in ("worker compiler flags", "library compiler flags", "link input", "unversioned import"):
            with self.subTest(kind=kind):
                altered = copy.deepcopy(fixture)
                provenance = altered["provenance.json"]
                if kind == "unversioned import":
                    dynamic = provenance["dynamic"]
                    dynamic["worker_symbol_table_text"] = dynamic["worker_symbol_table_text"].replace(
                        "@WIREHAIR_2.0", "", 1)
                    dynamic["worker_symbol_table_sha256"] = text_hash(dynamic["worker_symbol_table_text"])
                    provenance["dynamic_after"] = copy.deepcopy(dynamic)
                else:
                    build = provenance["build"]
                    if kind == "link input":
                        before = build["link_command"]
                        after = before.replace("-ldl", "-ldl unexpected.a")
                        build["link_command"] = after
                        build["link_command_sha256"] = text_hash(after)
                    else:
                        prefix = "library_" if kind.startswith("library") else ""
                        commands = build[prefix + "compile_commands"]
                        source = sorted(commands)[0]
                        before = commands[source]
                        after = before.replace("-O3", "-O0")
                        commands[source] = after
                        build[prefix + "compile_command_sha256"][source] = text_hash(after)
                    build["ninja_commands_text"] = build["ninja_commands_text"].replace(before, after)
                    build["ninja_commands_sha256"] = text_hash(build["ninja_commands_text"])
                    provenance["build_after"] = copy.deepcopy(build)
                with self.assertRaises(screen.ValidationError):
                    self.replay_fixture(altered)

    def test_replay_crosslinks_equal_before_after_receipts(self):
        fixture = bundle_fixture()
        for kind in ("worker hash", "git commit", "interpreter hash", "compiler hash"):
            with self.subTest(kind=kind):
                altered = copy.deepcopy(fixture)
                provenance = altered["provenance.json"]
                if kind == "git commit":
                    provenance["git"]["commit"] = "f" * 40
                    provenance["git"]["upstream_commit"] = "f" * 40
                    provenance["git_after"] = copy.deepcopy(provenance["git"])
                else:
                    name = {"worker hash": "worker", "interpreter hash": "controller_interpreter",
                            "compiler hash": "compiler"}[kind]
                    provenance["artifacts"][name]["sha256"] = "f" * 64
                    provenance["artifacts_after"] = copy.deepcopy(provenance["artifacts"])
                    if kind == "interpreter hash":
                        provenance["interpreter"]["artifact"]["sha256"] = "f" * 64
                        provenance["interpreter_after"] = copy.deepcopy(provenance["interpreter"])
                with self.assertRaises(screen.ValidationError):
                    self.replay_fixture(altered)

    def test_worker_marker_is_optional_only_before_start_and_binds_claim(self):
        fixture = bundle_fixture()
        claim = fixture["CLAIM"]
        claim_hash = screen.sha256_bytes(screen.canonical_bytes(claim) + b"\n")
        with tempfile.TemporaryDirectory(prefix="wh2-r1-marker-") as temporary:
            root = Path(temporary).resolve()
            fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
            try:
                self.assertIsNone(screen.read_worker_started(fd, claim, claim_hash))
                marker_path = root / screen.WORKER_STARTED_FILE
                marker_path.write_bytes(screen.canonical_bytes(fixture["WORKER_STARTED"]) + b"\n")
                os.chmod(str(marker_path), 0o400)
                marker, digest = screen.read_worker_started(fd, claim, claim_hash)
                self.assertEqual(marker, fixture["WORKER_STARTED"])
                self.assertEqual(digest, screen.file_sha256(marker_path))
                with self.assertRaises(screen.ValidationError):
                    screen.read_worker_started(fd, claim, "f" * 64)
                os.chmod(str(marker_path), 0o600)
                with self.assertRaises(screen.ValidationError):
                    screen.read_worker_started(fd, claim, claim_hash)
            finally:
                os.close(fd)

    def supervised_python(self, code, seconds=5.0, patches=()):
        """Substitute a harmless Python child for the real benchmark executable."""
        actual_popen = subprocess.Popen
        children = []

        def launch(argv, **kwargs):
            child = actual_popen([sys.executable, "-I", "-B", "-S", "-c", code], **kwargs)
            children.append(child)
            return child

        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(screen.subprocess, "Popen", side_effect=launch))
            for patch in patches:
                stack.enter_context(patch)
            try:
                return screen.run_worker(
                    "/synthetic/unused-worker", "/synthetic/libwirehair.so.2.0.0",
                    time.monotonic() + seconds, "1" * 64, "2" * 64, "3" * 64)
            finally:
                for child in children:
                    self.assertIsNotNone(child.poll(), "supervisor left a live child")
                    self.assertTrue(child.stdout.closed)
                    self.assertTrue(child.stderr.closed)

    def test_worker_supervision_drains_both_streams_and_reaps(self):
        raw, stderr, exit_code, _, timed_out = self.supervised_python(
            "import os; os.write(1,b'out'); os.write(2,b'err')")
        self.assertEqual((raw, stderr, exit_code, timed_out), (b"out", b"err", 0, False))

    def test_worker_supervision_bounds_streams_and_setup_failures(self):
        for fd, cap in ((1, "MAX_RAW_BYTES"), (2, "MAX_STDERR_BYTES")):
            with self.subTest(fd=fd):
                with self.assertRaisesRegex(screen.ValidationError, "exceeds its byte cap"):
                    self.supervised_python(
                        "import os; os.write({},b'x'*8192)".format(fd),
                        patches=(mock.patch.object(screen, cap, 1024),))
        with self.assertRaisesRegex(screen.ValidationError, "synthetic setup failure"):
            self.supervised_python("import time; time.sleep(30)", patches=(
                mock.patch.object(screen.os, "set_blocking", side_effect=OSError(
                    "synthetic setup failure")),))

    def test_worker_supervision_deadline_terminates_child(self):
        _, _, exit_code, elapsed, timed_out = self.supervised_python(
            "import time; time.sleep(30)", seconds=0.05)
        self.assertTrue(timed_out)
        self.assertNotEqual(exit_code, 0)
        self.assertLess(elapsed, 3.0)

    @unittest.skipUnless(hasattr(os, "fork") and Path("/proc").is_dir(), "Linux process groups")
    def test_worker_supervision_exception_kills_descendant_after_leader_exit(self):
        actual_read = os.read
        descendant = []

        def fail_after_descendant_ready(fd, count):
            data = actual_read(fd, count)
            if data.startswith(b"grandchild:"):
                descendant.append(int(data.split(b":", 1)[1]))
                raise OSError("synthetic pipe failure")
            return data

        code = (
            "import os, signal, time\n"
            "if os.fork() == 0:\n"
            "    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
            "    os.write(1, ('grandchild:%d\\n' % os.getpid()).encode())\n"
            "time.sleep(30)\n"
        )
        try:
            with self.assertRaisesRegex(screen.ValidationError, "synthetic pipe failure"):
                self.supervised_python(code, patches=(
                    mock.patch.object(screen.os, "read", side_effect=fail_after_descendant_ready),))
            self.assertEqual(len(descendant), 1)
            status_path = Path("/proc/{}/status".format(descendant[0]))
            for _ in range(100):
                try:
                    status = status_path.read_text(encoding="ascii")
                except FileNotFoundError:
                    break
                if "State:\tZ" in status:
                    break
                time.sleep(0.01)
            else:
                self.fail("supervisor left its descendant alive")
        finally:
            for pid in descendant:
                try:
                    os.kill(pid, 9)
                except ProcessLookupError:
                    pass


if __name__ == "__main__":
    unittest.main()
