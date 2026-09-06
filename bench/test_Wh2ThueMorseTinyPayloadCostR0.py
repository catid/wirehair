#!/usr/bin/env python3
"""Neutral JSONL/controller tests; no selected lookup, codec or native timing."""
import copy
import importlib.util
import io
import math
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location(
    "tiny_cost", Path(__file__).with_name("Wh2ThueMorseTinyPayloadCostR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

# Synthetic direct rows only: identity twice. Never construct the selected pair.
_lookup = bytearray(39936)
for _row in range(12):
    _lookup[6 * _row + _row % 6] = 1
NEUTRAL_LOOKUP = bytes(_lookup)
NEUTRAL_SHA = M.sha(NEUTRAL_LOOKUP)
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1)


def neutral_roster():
    """Literal frozen loop law, independent of the controller generator."""
    rows = []
    for rep in range(12):
        for step in range(2):
            for width_step in range(3):
                for metric_step in range(2):
                    for class_step in range(3):
                        order = (rep + step) % 2
                        comparison = (2*rep + step + width_step + metric_step + class_step) % 3
                        for position, side in enumerate(SIDES):
                            pair = rep % 4 if position < 2 else (position - 2)//2
                            rows.append(dict(index=len(rows), rep=rep, order=order,
                                width=(rep + step + width_step) % 3,
                                metric=(rep + step + width_step + metric_step) % 2,
                                **{"class": comparison}, position=position,
                                arm=((0, 0), (1, 1), (0, 1))[comparison][side ^ order],
                                q=((2*(rep + 12*pair) + 1)*1000000)//96))
    return rows


def fixture():
    identity = dict(derived=dict(family=26, model=8, stepping=1, full_apic_id=100,
        initial_apic_id_8=100, core_id=50, package_id=0, thread_id=0,
        threads_per_core=2, ccd_id=6, complex_id=6, logical_processors_per_package=128),
        before_cpu=50, after_cpu=50, requested_cpu=50, raw_capture_complete=True,
        semantic_validation_passed=True, singleton_affinity_verified=True,
        canonical_hex="0000", canonical_bytes=2, canonical_sha256=M.sha(bytes(2)),
        capabilities=dict.fromkeys(("leaf1_ecx", "leaf1_edx", "leaf6_eax", "leaf6_ecx",
            "leaf80000001_ecx", "leaf80000001_edx", "leaf80000008_ebx", "leaf80000021_eax",
            "max_basic_leaf", "max_extended_leaf"), 0))
    runtime = dict(polynomial=0x14d, ssse3=1, avx2=1, gfni=1, avx512=1, address=0x700000)
    fixtures, handles = [], {}
    for width, B in enumerate((2, 64, 1280)):
        source = bytes((37*i + i//11) & 255 for i in range(6*B))
        entries = []
        for rep in range(12):
            for step in range(2):
                arm = (rep + width + step) % 2
                for lane in range(2):
                    address = 0x1000000 + len(handles)*256
                    handles[(width, rep, lane, arm)] = address
                    entries.append(dict(rep=rep, lane=lane, arm=arm, address=address,
                        roundtrip=dict(create_code=0, feed_codes=[1]*5+[0], feed_count=6,
                                       recover_code=0, recover_bytes=6*B)))
        fixtures.append(dict(width_index=width, block_bytes=B, message_bytes=6*B,
            output_stride=B+128, source_hex=source.hex(), source_sha256=M.sha(source),
            source_address=0x100000 + width*0x10000 + 64,
            output_addresses=[0x400000 + width*0x100000 + lane*0x10000 for lane in range(2)],
            rows_hex=NEUTRAL_LOOKUP[:72].hex(), packets_hex=(source*2).hex(),
            packets_sha256=M.sha(source*2), handles=entries))
    start, cpu = 10**12, 5*10**12
    prelude = [start+1000, cpu+1000, start+1200, start+1300, cpu+1400, start+1500, [0]*4, [0]*4]
    header = dict(type="header", protocol=M.PROTOCOL, schema=M.PROTOCOL+".raw.v1",
        claim_sha256="ab"*32, target_cpu=50, worker_start_ns=start, worker_start_cpu_ns=cpu,
        identity_before=identity, runtime_before=runtime, clock_resolutions=[1, 1],
        lookup_sha256=NEUTRAL_SHA, lookup_address=0x200000, fixtures=fixtures,
        prelude=dict(iterations=1 << 20, seed=0x9e3779b97f4a7c15,
                     final_state=0x43935dad1647741b, observation=prelude))
    now = start + 5000
    t0 = (now//2000000 + 2)*2000000
    records, total = [], 0
    for coordinates in neutral_roster():
        row = dict(coordinates, type="record", called=True, checked=True)
        target = t0 + 2000000*row["index"] + row["q"]
        ctarget = cpu + target - start
        duration = (16000 if row["width"] == 0 else 20100) if row["arm"] else 20000
        row.update(target=target, ready=target-200000,
            wait=[target-199900, ctarget-199900, target, ctarget],
            observation=[target+500, ctarget+700, target+1000, target+1000+duration,
                         ctarget+1500+duration, target+2000+duration, [0]*4, [0]*4])
        metric, B = row["metric"], (2, 64, 1280)[row["width"]]
        address = 0x3000000 + row["index"]*256
        addresses = dict(count=64, sha256=M.sha(struct.pack("<Q", address)*64), min=address, max=address)
        if metric == 0:
            addresses = dict(count=0, sha256=M.sha(b""), min=None, max=None)
        row["work"] = dict(complete=True, create_calls=64*metric,
            encode_calls=384 if metric == 0 else 768, free_calls=64*metric,
            create_code=None if metric == 0 else 0, codes=[0]*(6 if metric == 0 else 12),
            required=[B]*(6 if metric == 0 else 12), written=[B]*(6 if metric == 0 else 12),
            handle=handles[(row["width"], row["rep"], row["order"], row["arm"])] if metric == 0 else None,
            addresses=addresses)
        records.append(row)
        total += duration
    footer = dict(type="footer", protocol=M.PROTOCOL, outcome="COMPLETE", failure=None,
        schedule_now_ns=now, t0=t0, records=4320, callbacks=4320, checked_callbacks=4320,
        create_calls=138240, encode_calls=2488320, free_calls=138240, sum_work_ns=total,
        sum_wait_wall_ns=4320*199900, sum_wait_cpu_ns=4320*199900,
        worker_end_ns=records[-1]["observation"][5]+1000,
        worker_end_cpu_ns=records[-1]["observation"][4]+1000,
        identity_after=copy.deepcopy(identity), runtime_after=dict(runtime),
        preflight=dict(encoder_creates=144, encoder_frees=144, encodes=1728,
                       decoder_creates=144, feeds=864, recovers=144, decoder_frees=144))
    return [header] + records + [footer]


def raw(items):
    return b"".join(M.canonical(item) for item in items)


class ReducerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.items = fixture()

    def setUp(self):
        patches = [mock.patch.object(M, "lookup_bytes", return_value=NEUTRAL_LOOKUP),
                   mock.patch.object(M, "LOOKUP_SHA", NEUTRAL_SHA)]
        for name in ("sibling", "current_receipt", "checked", "write_lookup", "run", "controller_run"):
            patches.append(mock.patch.object(M, name, side_effect=AssertionError("scientific work forbidden")))
        patches.append(mock.patch.object(M.A, "capture", side_effect=AssertionError("worker forbidden")))
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def reject(self, mutate):
        items = copy.deepcopy(self.items)
        mutate(items)
        with self.assertRaises(ValueError):
            M.reduce_raw(raw(items))

    def test_full_neutral_reduction_and_exact_work_accounting(self):
        data = raw(self.items)
        self.assertLess(len(data), 4*1024**2)
        result = M.reduce_raw(data)
        self.assertEqual(result["outcome"], "PASS")
        self.assertEqual((result["callbacks"], result["create_calls"], result["encode_calls"], result["free_calls"]),
                         (4320, 138240, 2488320, 138240))
        self.assertEqual(result["raw_sha256"], M.sha(data))
        self.assertEqual(len(result["statistics"]), 36)
        self.assertEqual(sum("control_pass" in row for row in result["statistics"]), 24)
        self.assertEqual(sum("treatment_pass" in row for row in result["statistics"]), 12)
        self.assertEqual(result["failed_controls"], [])
        self.assertEqual(result["failed_treatments"], [])
        for group in result["statistics"]:
            expected = (.8 if group["block_bytes"] == 2 else 1.005) if group["comparison"] == "BC" else 1
            self.assertAlmostEqual(group["estimate"]["ratio"], expected)
            self.assertEqual(len(group["replicate_logs"]), 12)
        for key in ("WH1_compared", "public_WH2_compared", "all_K_claimed", "recovery_rate_claimed",
                    "production_promotion_claimed"):
            self.assertIs(result[key], False)

    def test_exact_roster_and_all48_phases_per_logical_side(self):
        rows = list(M.roster())
        self.assertEqual(rows, neutral_roster())
        groups = {}
        for row in rows:
            if row["position"] >= 2:
                key = (row["width"], row["metric"], row["class"], row["order"],
                       SIDES[row["position"]] ^ row["order"])
                groups.setdefault(key, []).append(row["q"])
        self.assertEqual(len(groups), 72)
        expected = [((2*i+1)*1000000)//96 for i in range(48)]
        for phases in groups.values():
            self.assertEqual(sorted(phases), expected)
        self.assertEqual(sum(row["position"] < 2 for row in rows), 864)
        gaps = [2000000+b["q"]-a["q"] for a, b in zip(rows, rows[1:])]
        self.assertEqual((min(gaps), max(gaps)), (1250000, 2250000))

    def test_control_failure_overrides_passing_treatments(self):
        records = copy.deepcopy(self.items[1:-1])
        for row in records:
            if row["class"] == 0 and row["order"] == 1 and SIDES[row["position"]] ^ row["order"]:
                row["observation"][3] += 3000
        result = M.statistics(records)
        self.assertEqual(result["outcome"], "CONTROL_FAIL")
        self.assertEqual(len(result["failed_controls"]), 6)
        self.assertEqual(result["failed_treatments"], [])

    def test_b2_gain_and_each_wider_bound_are_required(self):
        for width, duration in ((0, 20000), (1, 21000), (2, 21000)):
            records = copy.deepcopy(self.items[1:-1])
            for row in records:
                if row["width"] == width and row["arm"] == 1:
                    row["observation"][3] = row["observation"][2] + duration
            result = M.statistics(records)
            self.assertEqual(result["outcome"], "FAIL")
            self.assertEqual(result["failed_controls"], [])
            self.assertEqual(result["failed_treatments"], [[width, metric, 2, order]
                for metric in range(2) for order in range(2)])

    def test_strict_boundaries_and_t11(self):
        bound = math.log1p(.02)
        with mock.patch.object(M.A, "confidence", return_value=dict(lower95_log=-bound, upper95_log=0)):
            result = M.statistics(self.items[1:-1])
            self.assertEqual(len(result["failed_controls"]), 24)
            self.assertEqual(len(result["failed_treatments"]), 4)
        with mock.patch.object(M.A, "confidence", return_value=dict(lower95_log=0, upper95_log=bound)):
            result = M.statistics(self.items[1:-1])
            self.assertEqual(len(result["failed_controls"]), 24)
            self.assertEqual(len(result["failed_treatments"]), 12)
        for values in ([0]*11, [0]*13, [float("nan")]*12):
            with self.assertRaises(ValueError):
                M.A.confidence(values)
        self.assertAlmostEqual(M.A.confidence([.01, -.01]*6)["upper95_log"],
                               2.200985160082949*math.sqrt(.0001/11))

    def test_missing_duplicate_and_boolean_roster_rejected(self):
        for mutate in (lambda x: x.pop(8), lambda x: x.__setitem__(8, copy.deepcopy(x[7])),
                       lambda x: x[1].update(index=False), lambda x: x[1].update(width=0.0),
                       lambda x: x[7].update(q=x[7]["q"]+1)):
            self.reject(mutate)

    def test_semantic_oracle_rejects_coordinated_wrong_hashes(self):
        for field in ("source", "packets"):
            def mutate(items, field=field):
                entry = items[0]["fixtures"][2]
                data = bytearray.fromhex(entry[field+"_hex"])
                data[-1] ^= 1
                entry.update({field+"_hex": data.hex(), field+"_sha256": M.sha(data)})
            self.reject(mutate)
        self.reject(lambda x: x[0]["fixtures"][0].update(rows_hex="00"*72))

    def test_status_lengths_counters_and_completion_rejected(self):
        for mutate in (
            lambda x: x[1]["work"]["codes"].__setitem__(0, False),
            lambda x: x[1]["work"]["written"].__setitem__(0, 1),
            lambda x: x[1]["work"]["required"].__setitem__(0, 2.0),
            lambda x: x[1]["work"].update(encode_calls=383),
            lambda x: x[1]["work"].update(complete=False),
            lambda x: x[1].update(called=False),
            lambda x: x[1].update(checked=False),
            lambda x: x[-1].update(encode_calls=2488319),
            lambda x: x[-1]["preflight"].update(feeds=863),
            lambda x: x[0]["fixtures"][0]["handles"][0]["roundtrip"].update(feed_count=5),
            lambda x: x[0]["fixtures"][0]["handles"][0]["roundtrip"]["feed_codes"].__setitem__(0, True),
            lambda x: x[-1].update(outcome="INVALID")):
            self.reject(mutate)

    def test_physical_handle_binding_and_spans(self):
        for mutate in (
            lambda x: x[1]["work"].update(handle=x[1]["work"]["handle"]+256),
            lambda x: x[0]["fixtures"][0]["handles"][1].update(
                address=x[0]["fixtures"][0]["handles"][0]["address"]),
            lambda x: x[0]["fixtures"][0]["handles"][0].update(arm=False),
            lambda x: x[0]["fixtures"][0].update(output_addresses=[0x400000]*2),
            lambda x: x[0]["fixtures"][2].update(output_addresses=[0x400000, 0x610000]),
            lambda x: x[0].update(lookup_address=0x100000),
            lambda x: x[0]["fixtures"][0].update(source_address=0x100041)):
            self.reject(mutate)
        items = copy.deepcopy(self.items)
        next(row for row in items[1:-1] if row["metric"] == 1)["work"]["addresses"]["count"] = 63
        with self.assertRaises(ValueError):
            M.reduce_raw(raw(items))

    def test_clocks_margins_and_full_lifetime(self):
        for mutate in (
            lambda x: x[1].update(ready=x[1]["target"]-99999),
            lambda x: x[1]["observation"].__setitem__(2, x[1]["target"]+5001),
            lambda x: x[1].update(observation=[]),
            lambda x: x[1]["observation"].__setitem__(7, [-1, 0, 0, 0]),
            lambda x: x[1]["wait"].__setitem__(1, x[0]["worker_start_cpu_ns"]-1),
            lambda x: x[-1].update(worker_end_ns=x[-2]["observation"][5]-1),
            lambda x: x[-1].update(worker_end_cpu_ns=x[0]["worker_start_cpu_ns"]+14*10**9+1),
            lambda x: x[-1].update(sum_work_ns=x[-1]["sum_work_ns"]+1),
            lambda x: x[-1].update(t0=x[-1]["t0"]+1)):
            self.reject(mutate)

    def test_identity_runtime_types_and_postbinding(self):
        for mutate in (lambda x: x[0].update(target_cpu=50.0),
            lambda x: x[0]["runtime_before"].update(gfni=True),
            lambda x: x[-1]["runtime_after"].update(address=0x710000),
            lambda x: x[0]["identity_before"]["derived"].update(thread_id=False),
            lambda x: x[-1]["identity_after"].update(capabilities={})):
            self.reject(mutate)

    def test_json_framing_and_duplicate_keys(self):
        for data in (b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":Infinity}'):
            with self.assertRaises(ValueError):
                M.decode(data)
        for data in (raw(self.items)[:-1], b"x"*(M.RAW_CAP+1)):
            with self.assertRaises(ValueError):
                M.reduce_raw(data)


class FakeChild:
    """Local tiny pipes only; lifecycle tests never create a subprocess."""
    def __init__(self, stdout=b"prefix\n", stderr=b""):
        self.pid, self.returncode, self.waited = 4242, None, False
        self.stdout, self.stderr = self.pipe(stdout), self.pipe(stderr)

    @staticmethod
    def pipe(data):
        read_fd, write_fd = os.pipe()
        os.write(write_fd, data)
        os.close(write_fd)
        return os.fdopen(read_fd, "rb")

    def poll(self):
        return self.returncode

    def wait(self, timeout=None):
        self.waited, self.returncode = True, 0
        return 0

    def kill(self):
        self.returncode = -9


class ReceiptTests(unittest.TestCase):
    def test_nonmapping_manifest_sources_rejected_before_pinning(self):
        names = ["gf256.cpp", "Wh2ThueMorseNativeCodec.cpp", "Wh2ThueMorseTinyPayloadR0.cpp",
            "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp", "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp",
            "Wh2ThueMorseTinyPayloadCostR0.cpp", "Wh2RdpruTargetIdentityV2.cpp",
            "Wh2PublicBorrowedTargetIdentity.cpp", "Wh2FrozenTrace.cpp",
            "Wh2ThueMorseNativeCodecTest.cpp", "Wh2ThueMorseTinyPayloadR0CodecTest.cpp",
            "Wh2ThueMorseTinyPayloadR0Test.cpp"]
        commands, bridge_arm = [], 0
        for name in names:
            path = M.ROOT / (name if name == "gf256.cpp" else "bench/"+name)
            command = "/usr/bin/g++ -std=c++11 -O3 -fno-lto -fPIC"
            if name == "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp":
                command += " -DWH2_TINY_COST_CANDIDATE="+str(bridge_arm)
                bridge_arm += 1
            commands.append(dict(file=str(path), command=command))
        includes = [M.ROOT/"gf256.h", M.ROOT/"bench/Wh2ThueMorseNativeCodec.h"]
        source_names = sorted({str(Path(item["file"]).relative_to(M.ROOT)) for item in commands} |
                              {str(path.relative_to(M.ROOT)) for path in includes})
        # The matching list passes the old sorted-key comparison but has no
        # .items(): that must be ValueError, not an unsealed AttributeError.
        for invalid in (source_names, None, True):
            with tempfile.TemporaryDirectory(prefix="wh2-tiny-cost-manifest-") as temporary:
                build = Path(temporary)
                files = {build/"CMakeCache.txt": b"WH2_TINY_COST_SANITIZERS:BOOL=OFF\n",
                    build/"lookup.bin": NEUTRAL_LOOKUP,
                    build/"compile_commands.json": M.canonical(commands),
                    build/"native-inputs.json": M.canonical(
                        dict(schema=M.PROTOCOL+".native-inputs.v1", sources=invalid))}
                inspections = [b"ninja: no work to do.\n",
                               ("".join("    "+str(path)+"\n" for path in includes)).encode()]
                with mock.patch.object(M, "checked", side_effect=inspections) as checked, \
                     mock.patch.object(M, "read_regular", side_effect=lambda path, cap: files[Path(path)]), \
                     mock.patch.object(M, "LOOKUP_SHA", NEUTRAL_SHA), \
                     mock.patch.object(M, "lookup_bytes", side_effect=AssertionError("lookup forbidden")), \
                     mock.patch.object(M, "pin", side_effect=AssertionError("real pins forbidden")) as pin:
                    with self.assertRaisesRegex(ValueError, "native input manifest source mapping"):
                        M.current_receipt(build, 123)
                self.assertEqual(checked.call_args_list, [
                    mock.call(["ninja", "-C", build, "-n"], 123),
                    mock.call(["ninja", "-C", build, "-t", "deps"], 123)])
                pin.assert_not_called()


class LifecycleTests(unittest.TestCase):
    def test_outer_selector_failure_does_not_spawn(self):
        with mock.patch.object(M.selectors, "DefaultSelector", side_effect=OSError("no selector")), \
             mock.patch.object(M.subprocess, "Popen") as spawn:
            with self.assertRaises(OSError):
                M.run(Path("/never-read"))
        spawn.assert_not_called()

    def test_outer_uses_new_module_and_kills_group_before_reaping(self):
        child = FakeChild(b"summary\n")
        stdout, stderr = mock.Mock(buffer=io.BytesIO()), mock.Mock(buffer=io.BytesIO())
        def kill_group(pid, sig):
            self.assertFalse(child.waited)
            self.assertEqual(pid, child.pid)
        with mock.patch.object(M.subprocess, "Popen", return_value=child) as spawn, \
             mock.patch.object(M.time, "monotonic", return_value=100), \
             mock.patch.object(M.os, "waitid", return_value=object()) as observed, \
             mock.patch.object(M.os, "killpg", side_effect=kill_group) as killed, \
             mock.patch.object(M.sys, "stdout", stdout), mock.patch.object(M.sys, "stderr", stderr):
            M.run(Path("/never-read"))
        self.assertEqual(spawn.call_args[0][0][1:3], [str(Path(M.__file__).resolve()), "_controller"])
        self.assertTrue(spawn.call_args[1]["start_new_session"])
        self.assertTrue(child.waited and child.stdout.closed and child.stderr.closed)
        self.assertEqual(killed.call_count, 1)
        self.assertTrue(observed.call_args[0][2] & os.WNOWAIT)
        self.assertEqual(stdout.buffer.getvalue(), b"summary\n")

    def test_capture_preserves_prefix_on_cap_and_seal_failures(self):
        for failure in ("cap", "fsync", "fchmod"):
            with tempfile.TemporaryDirectory(prefix="wh2-tiny-cost-capture-") as temporary:
                paths = [Path(temporary)/"raw", Path(temporary)/"err"]
                child = FakeChild(b"0123456789")
                cap = mock.patch.object(M.A, "RAW_CAP", 8) if failure == "cap" else \
                      mock.patch.object(M.os, failure, side_effect=OSError("seal failure"))
                with mock.patch.object(M.subprocess, "Popen", return_value=child), \
                     mock.patch.object(M.time, "monotonic", return_value=100), cap:
                    result = M.A.capture(Path("/never-executed"), "ab"*32, 110, paths)
                self.assertEqual(result[0], b"01234567" if failure == "cap" else b"0123456789")
                self.assertIsNotNone(result[3])
                self.assertTrue(child.waited and child.stdout.closed and child.stderr.closed)
                self.assertEqual(M.read_regular(paths[0], 100), result[0])

    def test_prelaunch_mismatch_does_not_claim_namespace(self):
        with tempfile.TemporaryDirectory(prefix="wh2-tiny-cost-prelaunch-") as temporary:
            root = Path(temporary)
            M.publish(root/"receipt", M.canonical(dict(build=str(root))))
            with mock.patch.object(M, "OUTPUT", root/"bundle"), \
                 mock.patch.object(M, "current_receipt", return_value={}), \
                 mock.patch.object(M.A, "capture") as capture:
                with self.assertRaises(ValueError):
                    M.controller_run(root/"receipt")
            capture.assert_not_called()
            self.assertFalse((root/"bundle").exists())

    def test_invalid_prefix_postpins_and_late_analysis_are_sealed(self):
        for failure in ("capture", "postpins", "late", "before-spool"):
            with tempfile.TemporaryDirectory(prefix="wh2-tiny-cost-controller-") as temporary:
                root = Path(temporary)
                output, receipt = root/"bundle", dict(build=str(root))
                M.publish(root/"receipt", M.canonical(receipt))
                now = [100]
                def capture(executable, claim, deadline, paths):
                    self.assertEqual(executable, root/M.TARGET)
                    self.assertEqual(list(paths), [output/"raw.jsonl", output/"stderr.txt"])
                    if failure == "before-spool":
                        raise OSError("injected before spool")
                    for path, data in zip(paths, (b"prefix\n", b"")):
                        M.publish(path, data)
                    return b"prefix\n", b"", 0, "capture failure" if failure == "capture" else None
                def reduce(data):
                    if failure == "late":
                        now[0] = 121
                    return dict(outcome="PASS", claim_sha256=M.sha(M.canonical(receipt)))
                with mock.patch.object(M, "OUTPUT", output), \
                     mock.patch.object(M.time, "monotonic", side_effect=lambda: now[0]), \
                     mock.patch.object(M, "current_receipt", side_effect=[receipt, {} if failure == "postpins" else receipt]), \
                     mock.patch.object(M.A, "capture", side_effect=capture), \
                     mock.patch.object(M, "reduce_raw", side_effect=reduce), mock.patch("builtins.print"):
                    M.controller_run(root/"receipt")
                summary = M.decode(M.read_regular(output/"summary.json", 10000))
                complete = M.decode(M.read_regular(output/"COMPLETE.json", 10000))
                self.assertEqual(summary["outcome"], "INVALID")
                self.assertEqual(complete["outcome"], "INVALID")
                self.assertIsNotNone(summary["failure"])
                self.assertEqual(M.read_regular(output/"raw.jsonl", 100), b"" if failure == "before-spool" else b"prefix\n")
                self.assertEqual(len(complete["files"]), 5)
                for member in complete["files"]:
                    path = Path(member["path"])
                    self.assertEqual(M.pin(path), member)
                    self.assertEqual(path.stat().st_mode & 0o777, 0o400)
                    self.assertEqual(path.stat().st_nlink, 1)


if __name__ == "__main__":
    unittest.main()
