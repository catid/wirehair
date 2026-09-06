#!/usr/bin/env python3
"""Neutral schema/statistics/publication tests; no codec, native clocks or ISA."""
import copy
import importlib.util
import io
import math
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location("repair_alignment_neutral_test", Path(__file__).with_name("Wh2RepairAlignmentR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
HISTORICAL_READER = M.historical_receipt


def identity():
    derived = dict(ccd_id=6, complex_id=6, core_id=50, family=26,
                   full_apic_id=100, initial_apic_id_8=100, logical_processors_per_package=128,
                   model=8, package_id=0, stepping=1, thread_id=0, threads_per_core=2)
    return dict(derived=derived, neutral_identity=True)


def address(value):
    return dict(address=value, mod32=value % 32, mod64=value % 64, mod4096=value % 4096)


def fixture(placement=1.1, noisy_AA=False):
    """Synthetic source/intermediates/column lists; never construct real rows."""
    source = bytes((37 * i + i // 11) & 255 for i in range(7680))
    master = bytes((i * 19 + i // 29 + 11) & 255 for i in range(46080))
    columns = [[j, 6 + j, 12 + j, 18 + j] for j in range(6)]
    packets = bytes(master[j * 1280 + offset] ^ master[(6 + j) * 1280 + offset] ^
                    master[(12 + j) * 1280 + offset] ^ master[(18 + j) * 1280 + offset]
                    for j in range(6) for offset in range(1280))
    profile = b"WHV2" + struct.pack("<HHQQI", 1, 32, 0x4b295bbb47f4f9c9, 7680, 1280) + bytes([23, 0, 0, 0])
    descriptor = dict(profile_hex=profile.hex(), profile_sha256=M.C.sha(profile),
                      source_sha256=M.C.sha(source), intermediate_sha256=M.C.sha(master),
                      message_bytes=7680, block_bytes=1280, source_count=6, precode_count=30,
                      intermediate_bytes=46080, source_policy=2, seed_attempt=23,
                      params=dict(block_count=6, staircase=6, dense_rows=12, heavy_rows=12, source_hits=2,
                                  dense_identity_corner=False, heavy_family=0, dense_anchors=0, seed=17),
                      config=dict(peel_seed=29, mix_count=3), runtime=dict(source_prime=7, precode_prime=31))
    metadata = dict(source_hex=source.hex(), source_sha256=M.C.sha(source), intermediate_hex=master.hex(),
                    intermediate_sha256=M.C.sha(master), expected_packets_hex=packets.hex(),
                    expected_packets_sha256=M.C.sha(packets), columns=columns, expected_operations=[4] * 6,
                    handles=[copy.deepcopy(descriptor), copy.deepcopy(descriptor)])
    addresses = dict(source=address(0x100000), master=address(0x200000),
                     handles=[address(0x300000), address(0x310000)],
                     intermediates=[address(0x400010), address(0x500030)],
                     carriers=[address(0x600000), address(0x700000)],
                     outputs=[address(0x800000), address(0x900000)], scratch=address(0xa00000),
                     public_function=4096, shadow_function=8192, runner_function=12288, prelude_function=16384)
    start, cpu = 10**12, 10**15
    prelude = dict(iterations=1 << 20, seed=M.PRELUDE_SEED, final_state=M.PRELUDE_FINAL,
                   m0=start + 100, c0=cpu, m1=start + 200, m2=start + 1000,
                   c1=cpu + 1000, m3=start + 1200, ru0=[0] * 4, ru1=[0] * 4)
    runtime = dict(polynomial=333, address=0xb00000, ssse3=1, avx2=1, gfni=1, avx512=1)
    raw = dict(protocol=M.PROTOCOL, schema=M.PROTOCOL + ".raw.v1", outcome="COMPLETE", failure=None,
               target_cpu=50, identity_before=identity(), identity_after=identity(),
               runtime_before=runtime, runtime_after=dict(runtime), worker_start_ns=start,
               monotonic_resolution_ns=1, thread_resolution_ns=1, metadata=metadata,
               handles_after=[copy.deepcopy(descriptor), copy.deepcopy(descriptor)],
               addresses=addresses, addresses_after=copy.deepcopy(addresses),
               preflight=dict(public=12, shadow=24, original_view=6, scalar=6), prelude=prelude,
               callbacks=3360, encode_calls=1290240, checked_packets=20160,
               checksum=M.expected_checksum(packets), sum_encode_wall_ns=0,
               elapsed_ns=500000000, records=[])
    side_by_position = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1)
    for coordinate in M.roster():
        index, rep, order, phase, comparison, position, cell = coordinate
        factor = (1.2 if cell >= 2 else 1) * (1.05 if cell % 2 else 1) if phase == 0 else (
            placement if cell % 2 else 1)
        if noisy_AA and phase == 0 and comparison == 0 and side_by_position[position] ^ order:
            factor *= 1.1
        wall = round(20000 * factor)
        m0, c0 = start + 10000 + index * 100000, cpu + 10000 + index * 100000
        usage = [index // 100, 0, index // 200, index // 300]
        after = [(index + 1) // 100, 0, (index + 1) // 200, (index + 1) // 300]
        raw["records"].append(list(coordinate) + [m0, c0, m0 + 100, m0 + 100 + wall,
            c0 + wall + 200, m0 + wall + 500, usage, after, 384, [0] * 6, [1280] * 6,
            [None] * 6 if phase == 0 else [256] * 6])
        raw["sum_encode_wall_ns"] += wall
    return raw


class AlignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = fixture()

    def setUp(self):
        for owner, name in ((M, "historical_receipt"), (M.D, "historical_receipt"),
                            (M.M, "authenticated_lookup"), (M.M, "run"), (M.M, "build_only"),
                            (M.D, "run"), (M.D, "build_only"), (M.N, "checked"), (M.C, "capture_worker")):
            patch = mock.patch.object(owner, name, side_effect=AssertionError("forbidden actual work"))
            patch.start()
            self.addCleanup(patch.stop)
        patch = mock.patch.object(M.R.H, "validate_target_identity", return_value=None)
        patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch.object(M.time, "monotonic", return_value=100.0)
        patch.start()
        self.addCleanup(patch.stop)

    def test_exact_roster_and_log_contrasts(self):
        roster = list(M.roster())
        self.assertEqual(len(roster), 3360)
        self.assertEqual(roster[0], (0, 0, 0, 0, 0, 0, 0))
        self.assertEqual(roster[-1], (3359, 11, 0, 1, 4, 9, 1))
        result = M.reduce_raw(self.raw)
        self.assertEqual(result["outcome"], "PLACEMENT_SUPPORTED")
        self.assertEqual(len(result["panels"]), 336)
        self.assertEqual(len(result["records"]), 3360)
        self.assertEqual(sum(row["warmup"] for row in result["records"]), 672)
        self.assertEqual(len(result["statistics"]), 28)
        self.assertEqual(sum("AA_pass" in row for row in result["statistics"]), 16)
        self.assertEqual(result["failed_controls"], [])
        for item in result["statistics"]:
            if item["phase"] == 1 and item["comparison"] >= 4:
                self.assertAlmostEqual(item["estimate"]["ratio"], 1.1)
            self.assertEqual(item["estimate"]["replicates"], 12)
        for item in result["handle_output_interaction"]:
            self.assertAlmostEqual(item["estimate"]["ratio"], 1)
        self.assertEqual(sum(row["counter_deltas"]["minflt"] for row in result["records"]), 33)
        self.assertEqual(result["startup_to_first_panel_ns"], 10000)
        self.assertEqual(result["conditioning"]["prelude_wall_ns"], 800)
        self.assertNotIn("encode_wall_ns", result["conditioning"])
        for name in ("speed_claimed", "WH1_compared", "production_promotion_claimed"):
            self.assertFalse(result[name])

    def test_global_control_failure_overrides_all_supported_treatments(self):
        result = M.reduce_raw(fixture(noisy_AA=True))
        self.assertEqual(result["outcome"], "CONTROL_FAIL")
        self.assertEqual(result["failed_controls"], [[0, 0, 0], [0, 0, 1]])
        self.assertTrue(all(item["placement_supported"] for item in result["statistics"] if "placement_supported" in item))
        self.assertEqual(len(result["records"]), 3360)

    def test_inconclusive_and_opposite_are_not_support(self):
        for factor in (1.0, .9):
            result = M.reduce_raw(fixture(placement=factor))
            self.assertEqual(result["outcome"], "PLACEMENT_INCONCLUSIVE")
            self.assertEqual(result["failed_controls"], [])

    def test_both_carriers_and_both_orders_required(self):
        raw = copy.deepcopy(self.raw)
        for row in raw["records"]:
            if row[2] == 1 and row[3] == 1 and row[4] == 5 and row[6] == 3:
                raw["sum_encode_wall_ns"] -= row[10] - row[9] - 20000
                row[10] = row[9] + 20000
        self.assertEqual(M.reduce_raw(raw)["outcome"], "PLACEMENT_INCONCLUSIVE")

    def test_replication_not_slot_pseudoreplication_and_strict_AA(self):
        values = [(-1)**i * .01 for i in range(12)]
        result = M.confidence(values)
        expected = M.T11 * math.sqrt(sum(x * x for x in values) / 11 / 12)
        self.assertAlmostEqual(result["lower95_log"], -expected)
        self.assertAlmostEqual(result["upper95_log"], expected)
        for count in (11, 13):
            with self.assertRaises(ValueError):
                M.confidence([0] * count)
        boundary = dict(replicates=12, mean_log=0, lower95_log=-math.log1p(.02),
                        upper95_log=0, ratio=1, lower95=1/1.02, upper95=1)
        with mock.patch.object(M, "confidence", return_value=boundary):
            result = M.reduce_raw(self.raw)
        self.assertEqual(len(result["failed_controls"]), 16)

    def test_warmup_tail_and_all_counter_samples_retained(self):
        raw = copy.deepcopy(self.raw)
        row = raw["records"][-1]
        for index in (10, 11, 12):
            row[index] += 1000000
        raw["sum_encode_wall_ns"] += 1000000
        result = M.reduce_raw(raw)
        self.assertEqual(len(result["records"]), 3360)
        self.assertEqual(result["records"][-1]["encode_wall_ns"], 1022000)
        self.assertTrue(all(row["thread_span_ns"] >= row["encode_wall_ns"] for row in result["records"]))
        self.assertNotIn("off_cpu_ns", result["records"][0])

    def test_scalar_oracle_and_compact_checksum_are_independent(self):
        packets = bytes.fromhex(self.raw["metadata"]["expected_packets_hex"])
        self.assertEqual(M.verify_metadata(self.raw["metadata"]), packets)
        inner = 1469598103934665603
        for byte in packets:
            inner = ((inner * 1099511628211) & ((1 << 64) - 1)) ^ byte
        outer = 1469598103934665603
        for byte in inner.to_bytes(8, "little") * 3360:
            outer = ((outer * 1099511628211) & ((1 << 64) - 1)) ^ byte
        self.assertEqual(M.expected_checksum(packets), outer)
        state = M.PRELUDE_SEED
        for expected in (0xdc1b77ae0bf34dad, 0x64f0eeb9026e6076, 0x7b07ce91e5906136):
            state ^= (state << 13) & M.U64
            state ^= state >> 7
            state ^= (state << 17) & M.U64
            self.assertEqual(state, expected)
        self.assertEqual(M.PRELUDE_FINAL, 4869338620102145051)

    def test_schema_boolean_float_missing_truncated_and_status_mutants(self):
        mutations = [(["callbacks"], True), (["encode_calls"], 1.0), (["checked_packets"], 0),
            (["failure"], "bad"), (["outcome"], "PASS"), (["records", 0, 0], False),
            (["records", 0, 3], 1), (["records", 0, 6], 1), (["records", 0, 7], True),
            (["records", 0, 8], 1.0), (["records", 0, 13, 0], False),
            (["records", 0, 15], 383), (["records", 0, 16, 0], 1),
            (["records", 0, 17, 0], 1279), (["records", 0, 18, 0], 256),
            (["sum_encode_wall_ns"], 200000001), (["sum_encode_wall_ns"], 1),
            (["elapsed_ns"], 10**10), (["elapsed_ns"], 1), (["worker_start_ns"], 0),
            (["prelude", "iterations"], 1024), (["prelude", "final_state"], 0),
            (["prelude", "m2"], None), (["prelude", "ru1", 0], True),
            (["preflight", "public"], 6), (["checksum"], 0),
            (["runtime_after", "address"], 0), (["identity_after", "derived", "core_id"], 51),
            (["addresses", "carriers", 0, "mod64"], 16),
            (["addresses_after", "source", "address"], 1),
            (["metadata", "columns", 0, 0], True), (["metadata", "columns", 0, 1], 0),
            (["metadata", "columns", 0, 2], 36), (["metadata", "expected_operations", 0], 4.0),
            (["metadata", "expected_packets_sha256"], "0" * 64),
            (["metadata", "source_sha256"], "0" * 64),
            (["handles_after", 0, "seed_attempt"], 1),
            (["monotonic_resolution_ns"], 0)]
        for path, value in mutations:
            raw = copy.deepcopy(self.raw)
            target = raw
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                M.reduce_raw(raw)
        for raw in (dict(self.raw, extra=1), dict(self.raw, records=self.raw["records"][:-1])):
            with self.assertRaises(ValueError):
                M.reduce_raw(raw)
        raw = dict(self.raw)
        del raw["records"]
        with self.assertRaises(ValueError):
            M.reduce_raw(raw)

    def test_profile_equation_and_packet_corruption_rejected(self):
        for path, value in ((["params", "dense_identity_corner"], 0), (["params", "dense_anchors"], 1),
                            (["config", "mix_count"], 2), (["runtime", "source_prime"], 5),
                            (["source_policy"], 1), (["profile_hex"], "00" * 32)):
            raw = copy.deepcopy(self.raw)
            for item in raw["metadata"]["handles"] + raw["handles_after"]:
                target = item
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                M.reduce_raw(raw)
        raw = copy.deepcopy(self.raw)
        raw["metadata"]["expected_packets_hex"] = "ff" + raw["metadata"]["expected_packets_hex"][2:]
        with self.assertRaisesRegex(ValueError, "scalar XOR"):
            M.reduce_raw(raw)

    def test_clock_counter_cross_record_order_and_independent_epochs(self):
        mutations = [(0, 9, self.raw["records"][0][7] - 1), (0, 10, self.raw["records"][0][9]),
                     (0, 11, self.raw["records"][0][8] - 1), (0, 12, self.raw["records"][0][10] - 1),
                     (1, 7, self.raw["records"][0][12] - 1), (1, 8, self.raw["records"][0][11] - 1),
                     (100, 13, [0] * 4), (100, 14, [0] * 4)]
        for index, column, value in mutations:
            raw = copy.deepcopy(self.raw)
            raw["records"][index][column] = value
            with self.subTest(index=index, column=column), self.assertRaises(ValueError):
                M.reduce_raw(raw)
        raw = copy.deepcopy(self.raw)
        raw["records"][-1][8] = raw["prelude"]["c0"] + raw["elapsed_ns"]
        raw["records"][-1][11] = raw["records"][-1][8] + 1
        with self.assertRaisesRegex(ValueError, "lifetime"):
            M.reduce_raw(raw)
        result = M.reduce_raw(self.raw)
        self.assertEqual(result["records"][0]["thread_span_ns"], result["records"][0]["encode_wall_ns"] + 200)

    def test_storage_crossovers_and_guards_have_distinct_spans(self):
        for key in ("carriers", "outputs", "intermediates", "handles"):
            addresses = copy.deepcopy(self.raw["addresses"])
            addresses[key][1] = copy.deepcopy(addresses[key][0])
            with self.subTest(key=key), self.assertRaises(ValueError):
                M.verify_addresses(addresses)
        addresses = copy.deepcopy(self.raw["addresses"])
        addresses["carriers"][0] = address(addresses["carriers"][0]["address"] + 16)
        with self.assertRaisesRegex(ValueError, "page-aligned"):
            M.verify_addresses(addresses)

    def test_deadline_before_work(self):
        with self.assertRaises(TimeoutError):
            M.reduce_raw(self.raw, deadline=0)
        with self.assertRaises(TimeoutError):
            M.expected_checksum(bytes(7680), deadline=0)

    def test_build_facade_flags_and_exact_link_inputs(self):
        commands = M.build_commands(Path("/tmp/neutral-alignment-build"))
        self.assertEqual(len(commands), 3)
        self.assertEqual([Path(command[command.index("-c") + 1]).name for command in commands[:2]],
                         ["Wh2RepairAlignmentR0Bridge.cpp", "Wh2RepairAlignmentR0.cpp"])
        for flag in ("-DWIREHAIR_BUILDING=1", "-fPIC"):
            self.assertIn(flag, commands[0])
            self.assertNotIn(flag, commands[1])
        for command in commands[:2]:
            for flag in ("-DWIREHAIR_STATIC=1", "-march=native", "-std=gnu++11", "-Werror"):
                self.assertIn(flag, command)
        self.assertEqual(sum(part.endswith("libwirehair.a") for part in commands[-1]), 1)
        for command in commands:
            self.assertFalse(any("NativeCodec" in part or "MulRowsR0.o" in part for part in command))
            self.assertNotIn("-flto", command)
            self.assertNotIn("--worker", command)

    def test_link_map_rejects_duplicate_facade_or_candidate(self):
        build = Path("/tmp/neutral-build")
        content = (str(build / "Wh2RepairAlignmentR0Bridge.o") + "\n" + str(M.M.BASE / "libwirehair.a") +
                   "\nwirehair_v2_encode\n").encode("ascii")
        M.validate_link_map(content, build)
        for forbidden in (b"libwirehair.a(WirehairV2Profile.cpp.o)", b"Wh2ThueMorseNativeCodec.o", b"Wh2ThueMorseMulRowsR0.o"):
            with self.assertRaises(ValueError):
                M.validate_link_map(content + forbidden, build)
        with self.assertRaises(ValueError):
            M.validate_link_map(b"empty", build)

    def test_stale_source_and_object_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source", Path(temp) / "object"
            M.C.write_new(source, b"neutral source")
            M.C.write_new(output, b"neutral object")
            manifest = dict(inputs=[M.N.pin(source)], outputs=[M.N.pin(output)])
            M.M.validate_build_manifest(manifest, {source}, {output})
            for key in ("inputs", "outputs"):
                forged = copy.deepcopy(manifest)
                forged[key][0]["sha256"] = "0" * 64
                with self.assertRaises(ValueError):
                    M.M.validate_build_manifest(forged, {source}, {output})

    def test_historical_source_and_dependency_pins_are_all_checked(self):
        source = b"synthetic historical source"
        pin = dict(path="/tmp/neutral-dependency", bytes=20, sha256="2" * 64)
        old = dict(sources={name: M.C.sha(source) for name in M.D.SOURCES},
                   dependencies=[pin], build_files=[pin], helpers=[pin], interpreter=pin)
        data = M.C.canonical(old)

        def read(path, *args, **kwargs):
            return data if path == M.HISTORICAL_RECEIPT else source

        with mock.patch.object(M, "HISTORICAL_SHA", M.C.sha(data)), \
                mock.patch.object(M.C, "read_regular", side_effect=read), \
                mock.patch.object(M.D, "historical_receipt", return_value={}) as transitive, \
                mock.patch.object(M.N, "pin", return_value=pin) as pinner:
            self.assertEqual(HISTORICAL_READER(), old)
            self.assertEqual(pinner.call_count, 4)
            transitive.assert_called_once_with(None)
            pinner.return_value = dict(pin, sha256="3" * 64)
            with self.assertRaisesRegex(ValueError, "dependency pin"):
                HISTORICAL_READER()
            pinner.return_value = pin
            source = b"changed historical source"
            with self.assertRaisesRegex(ValueError, "source pin"):
                HISTORICAL_READER()

    def run_synthetic(self, output, receipt_path, raw, duration=.6, post_error=False,
                      worker_code=0, worker_stderr=b""):
        receipt = dict(build_directory="/tmp/neutral-unused", worker_argv=["never-executed", "--worker"])
        M.C.write_new(receipt_path, M.C.canonical(receipt))
        raw_bytes = M.C.canonical(raw) + b"\n"
        observed = dict(stdout=raw_bytes, stderr=worker_stderr, failure=None,
                        returncode=worker_code, elapsed_seconds=duration)
        with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current_receipt",
                side_effect=[receipt, ValueError("post-pin") if post_error else receipt]), \
                mock.patch.object(M.C, "capture_worker", return_value=observed) as capture, \
                mock.patch("sys.stdout", new_callable=io.StringIO):
            code = M.run(receipt_path)
            self.assertEqual(capture.call_count, 1)
        summary = M.C.strict_json(M.C.read_regular(output / "summary.json", 65536))
        self.assertEqual(M.C.read_regular(output / "raw.json", M.RAW_CAP), raw_bytes)
        return code, summary, receipt

    def test_immutable_one_shot_publication(self):
        with tempfile.TemporaryDirectory() as temp:
            output, receipt_path = Path(temp) / "out", Path(temp) / "receipt"
            code, summary, receipt = self.run_synthetic(output, receipt_path, self.raw)
            self.assertEqual(code, 0)
            self.assertEqual(summary["outcome"], "PLACEMENT_SUPPORTED")
            manifest = M.C.strict_json(M.C.read_regular(output / "COMPLETE.json", 65536))
            self.assertEqual(set(manifest["files"]), {"CLAIM.json", "raw.json", "stderr.txt", "analysis.json", "summary.json"})
            for name, pin in manifest["files"].items():
                data = M.C.read_regular(output / name, M.TOTAL_CAP)
                self.assertEqual(pin, dict(bytes=len(data), sha256=M.C.sha(data)))
                self.assertEqual((output / name).stat().st_mode & 0o777, 0o400)
                self.assertEqual((output / name).stat().st_nlink, 1)
            with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current_receipt", return_value=receipt), \
                    mock.patch.object(M.C, "capture_worker") as capture, self.assertRaises(FileExistsError):
                M.run(receipt_path)
            capture.assert_not_called()

    def test_invalid_partial_prefix_retained_without_analysis(self):
        raw = dict(self.raw, outcome="INVALID", failure="neutral callback failure", callbacks=1,
                   records=[self.raw["records"][0], [1, 0, 0, 0, 0, 1, 0] + [None] * 12])
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "out"
            code, summary, _ = self.run_synthetic(output, Path(temp) / "receipt", raw, worker_code=1)
            self.assertEqual(code, 1)
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertFalse((output / "analysis.json").exists())

    def test_observed_internal_duration_and_postpin_failures(self):
        for duration, post_error in ((10, False), (.01, False), (.6, True)):
            with self.subTest(duration=duration, post_error=post_error), tempfile.TemporaryDirectory() as temp:
                code, summary, _ = self.run_synthetic(Path(temp) / "out", Path(temp) / "receipt", self.raw,
                                                       duration=duration, post_error=post_error)
                self.assertEqual(code, 1)
                self.assertEqual(summary["outcome"], "INVALID")

    def test_worker_stderr_and_controller_deadline_cannot_be_success(self):
        with tempfile.TemporaryDirectory() as temp:
            code, summary, _ = self.run_synthetic(Path(temp) / "out", Path(temp) / "receipt",
                                                 self.raw, worker_stderr=b"neutral failure")
            self.assertEqual(code, 1)
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertFalse((Path(temp) / "out/analysis.json").exists())
        # Advance only after successful synthetic reduction/post-pin work.
        # Every duration stays in its own clock domain; no real clock call.
        with tempfile.TemporaryDirectory() as temp:
            output, receipt_path = Path(temp) / "out", Path(temp) / "receipt"
            receipt = dict(build_directory="/tmp/neutral-unused", worker_argv=["never-executed", "--worker"])
            M.C.write_new(receipt_path, M.C.canonical(receipt))
            observed = dict(stdout=M.C.canonical(self.raw) + b"\n", stderr=b"", failure=None,
                            returncode=0, elapsed_seconds=.6)
            stage = [100.0]

            def pins(*args, **kwargs):
                if stage[0] == 100.0:
                    stage[0] = 100.1
                else:
                    stage[0] = 119.6
                return receipt

            with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current_receipt", side_effect=pins), \
                    mock.patch.object(M.C, "capture_worker", return_value=observed), \
                    mock.patch.object(M.time, "monotonic", side_effect=lambda: stage[0]), \
                    mock.patch("sys.stdout", new_callable=io.StringIO):
                self.assertEqual(M.run(receipt_path), 1)
            summary = M.C.strict_json(M.C.read_regular(output / "summary.json", 65536))
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertEqual(summary["failure"], "controller deadline")


if __name__ == "__main__":
    unittest.main()
