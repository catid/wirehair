#!/usr/bin/env python3
"""Synthetic-only controller/reducer tests: no codec, target or native timing."""
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

SPEC = importlib.util.spec_from_file_location("aligned_cost", Path(__file__).with_name("Wh2AlignedIntermediateCostR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def fixture():
    source = bytes((37*i + i//11) & 255 for i in range(7680))
    master = source + bytes(38400)
    packets = source * 2
    profile = b"WHV2" + struct.pack("<HHQQI", 1, 32, 0x4b295bbb47f4f9c9, 7680, 1280) + bytes(4)
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
    snapshots, by_key = [], {}
    for rep in range(12):
        for i in range(3):
            arm = (rep + i) % 3
            for lane in range(2):
                index = len(snapshots)
                snapshot = dict(arm=arm, rep=rep, lane=lane, handle=0x800000 + 256*index,
                    source=0x100000, intermediate=0x1000000 + 0x20000*index if arm else None,
                    intermediate_bytes=46080 if arm else None, intermediate_capacity=46080 if arm else None,
                    source_policy=2 if arm else None, profile_hex=profile.hex() if arm else None,
                    source_sha256=M.sha(source), intermediate_sha256=M.sha(master) if arm else None)
                snapshots.append(snapshot)
                by_key[(rep, lane, arm)] = snapshot
    start, cpu = 10**12, 5*10**12
    prelude = [start+1000, cpu+1000, start+1200, start+1300, cpu+1400, start+1500, [0]*4, [0]*4]
    header = dict(type="header", protocol=M.PROTOCOL, schema=M.PROTOCOL+".raw.v1", claim_sha256="ab"*32,
        target_cpu=50, worker_start_ns=start, worker_start_cpu_ns=cpu, identity_before=identity,
        runtime_before=runtime, clock_resolutions=[1, 1], source_hex=source.hex(), source_sha256=M.sha(source),
        outputs=[0x900000, 0xa00000], prelude=dict(iterations=1 << 20, seed=0x9e3779b97f4a7c15,
            final_state=0x43935dad1647741b, observation=prelude),
        preflight=dict(packets_hex=[packets.hex()]*3, packet_sha256=[M.sha(packets)]*3,
            profiles_hex=[None, profile.hex(), profile.hex()], intermediate_hex=[master.hex()]*2,
            intermediate_sha256=[M.sha(master)]*2,
            columns=[[[j % 6, 6+j % 6, 12+j % 6, 18+j % 6] for j in range(12)]]*2,
            decode_statuses=[[1]*5+[0]]*3, first_success=[5]*3, snapshots=snapshots))
    now = start + 5000
    t0 = (now // 2000000 + 2) * 2000000
    records, total, wait_wall, wait_cpu, calls = [], 0, 0, 0, 0
    for coordinates in M.roster():
        row = dict(coordinates, type="record")
        index, arm, metric = row["index"], row["arm"], row["metric"]
        target = t0 + 2000000 * index + row["q"]
        ctarget = cpu + target - start
        duration = (20000, 19000, 16000)[arm]
        row.update(target=target, ready=target-200000,
                   wait=[target-199900, ctarget-199900, target, ctarget],
                   observation=[target+500, ctarget+700, target+1000, target+1000+duration,
                                ctarget+1500+duration, target+2000+duration, [0]*4, [0]*4], checked=True, snapshot=None)
        natural = by_key[(row["rep"], row["order"], arm)]
        work = dict(calls=(1, 384, 14, 9)[metric], create_code=None if metric == 1 else 0,
                    recover_code=0 if metric == 3 else None,
                    codes=([1]*5+[0]) if metric == 3 else [0]*(0, 6, 12)[metric],
                    lengths=[1280]*(0, 6, 12, 6)[metric], first_success=5 if metric == 3 else None,
                    bytes=7680 if metric == 3 else None,
                    handle=natural["handle"] if metric == 1 else 0x90000000+index*256,
                    freed=metric in (2, 3), complete=True)
        row["work"] = work
        if metric == 0:
            snapshot = dict(natural, handle=work["handle"])
            if arm:
                snapshot["intermediate"] = 0x30000000 + 0x20000*index
            row["snapshot"] = snapshot
        records.append(row)
        total += duration
        wait_wall += 199900
        wait_cpu += 199900
        calls += work["calls"]
    footer = dict(type="footer", protocol=M.PROTOCOL, outcome="COMPLETE", failure=None,
        schedule_now_ns=now, t0=t0, records=4800, callbacks=4800, work_calls=calls,
        checked_callbacks=4800, sum_work_ns=total, sum_wait_wall_ns=wait_wall, sum_wait_cpu_ns=wait_cpu,
        worker_end_ns=records[-1]["observation"][5]+1000,
        worker_end_cpu_ns=records[-1]["observation"][4]+1000,
        identity_after=copy.deepcopy(identity), runtime_after=dict(runtime))
    return [header] + records + [footer]


def raw(items):
    return b"".join(M.canonical(item) for item in items)


class CostTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.items = fixture()

    def setUp(self):
        # A reducer/schema test must never accidentally become an experiment.
        for name in ("current_receipt", "capture", "run", "controller_run", "checked"):
            patch = mock.patch.object(M, name, side_effect=AssertionError("native work forbidden"))
            patch.start()
            self.addCleanup(patch.stop)

    def test_full_roster_and_conditional_statistics(self):
        result = M.reduce_raw(raw(self.items))
        self.assertEqual(result["outcome"], "PASS")
        self.assertEqual(len(result["statistics"]), 40)
        self.assertEqual(sum("control_pass" in x for x in result["statistics"]), 24)
        self.assertEqual(result["failed_controls"], [])
        for row in result["statistics"]:
            expected = .8 if row["comparison"] == "WA" else 16/19 if row["comparison"] == "PA" else 1
            self.assertAlmostEqual(row["estimate"]["ratio"], expected)

    def test_every_group_has_all48_phases(self):
        groups = {}
        rows = list(M.roster())
        self.assertEqual(len(rows), 4800)
        for row in rows:
            if row["position"] >= 2:
                key = (row["metric"], row["class"], row["order"], M.SIDES[row["position"]] ^ row["order"])
                groups.setdefault(key, []).append(row["q"])
        expected = [((2*i+1)*1000000)//96 for i in range(48)]
        self.assertEqual(len(groups), 80)
        for phases in groups.values():
            self.assertEqual(sorted(phases), expected)
        gaps = [2000000+b["q"]-a["q"] for a, b in zip(rows, rows[1:])]
        self.assertEqual(min(gaps), 1250000)
        self.assertEqual(max(gaps), 2250000)

    def test_control_failure_overrides_all_treatments(self):
        records = copy.deepcopy(self.items[1:-1])
        for row in records:
            if row["class"] == 0 and row["order"] == 1 and M.SIDES[row["position"]] ^ row["order"]:
                row["observation"][3] += 3000
        result = M.statistics(records)
        self.assertEqual(result["outcome"], "CONTROL_FAIL")
        self.assertEqual(len(result["failed_controls"]), 4)
        self.assertEqual(result["failed_treatments"], [])

    def test_no_pooling_or_repair_only_promotion(self):
        records = copy.deepcopy(self.items[1:-1])
        for row in records:
            if row["arm"] == 2 and row["metric"] == 0:
                row["observation"][3] = row["observation"][2]+25000
        result = M.statistics(records)
        self.assertEqual(result["outcome"], "FAIL")
        self.assertEqual(result["failed_controls"], [])
        self.assertEqual(len(result["failed_treatments"]), 4)

    def test_strict_control_boundary_and_t11(self):
        with mock.patch.object(M, "confidence", return_value=dict(lower95_log=-math.log1p(.02),
                upper95_log=0, ratio=1, lower95=1/1.02, upper95=1)):
            self.assertEqual(len(M.statistics(self.items[1:-1])["failed_controls"]), 24)
        for count in (11, 13):
            with self.assertRaises(ValueError):
                M.confidence([0]*count)
        values = [.01, -.01]*6
        result = M.confidence(values)
        self.assertAlmostEqual(result["upper95_log"], M.T11 * math.sqrt(.0001/11))

    def test_mutation_rejections(self):
        mutations = (
            lambda x: x.pop(8),
            lambda x: x[7].update(index=True),
            lambda x: x[7].update(q=x[7]["q"]+1),
            lambda x: x[7].update(ready=x[7]["target"]-99999),
            lambda x: x[7]["observation"].__setitem__(2, x[7]["target"]+5001),
            lambda x: x[7].update(checked=False),
            lambda x: x[7]["work"].update(complete=False),
            lambda x: x[-1].update(callbacks=4799),
            lambda x: x[-1].update(outcome="INVALID"),
            lambda x: x[0]["preflight"]["snapshots"][4].update(intermediate=17),
            lambda x: x[0]["preflight"].update(first_success=[6, 5, 5]),
            lambda x: x[0]["preflight"]["columns"][0][0].__setitem__(0, 35),
        )
        for mutate in mutations:
            items = copy.deepcopy(self.items)
            mutate(items)
            with self.assertRaises(ValueError):
                M.reduce_raw(raw(items))

    def test_json_rejects_duplicates_nonfinite_and_truncation(self):
        for data in (b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":Infinity}'):
            with self.assertRaises(ValueError):
                M.decode(data)
        with self.assertRaises(ValueError):
            M.reduce_raw(raw(self.items)[:-1])

    def test_numeric_impersonation_and_output_overlap(self):
        mutations = (
            lambda x: x[0]["runtime_before"].update(gfni=True),
            lambda x: x[-1]["runtime_after"].update(gfni=1.0),
            lambda x: x[0]["identity_before"]["derived"].update(thread_id=False),
            lambda x: x[0].update(target_cpu=50.0),
            lambda x: x[0]["preflight"]["snapshots"][0].update(arm=False),
            lambda x: x[0]["preflight"]["snapshots"][2].update(intermediate_capacity=46080.0),
            lambda x: x[0]["preflight"]["decode_statuses"][0].__setitem__(0, True),
            lambda x: x[1]["work"].update(create_code=False),
            lambda x: x[0].update(outputs=[0x900000, 0x901000]),
            lambda x: x[0].update(outputs=[0x1040000, 0x900000]),
        )
        for mutate in mutations:
            items = copy.deepcopy(self.items)
            mutate(items)
            with self.assertRaises(ValueError):
                M.reduce_raw(raw(items))

    def test_second_arm_column_type_and_missing_capabilities(self):
        items = copy.deepcopy(self.items)
        items[0]["preflight"]["columns"][1] = copy.deepcopy(items[0]["preflight"]["columns"][1])
        items[0]["preflight"]["columns"][1][0][0] = 0.0
        with self.assertRaises(ValueError):
            M.reduce_raw(raw(items))
        items = copy.deepcopy(self.items)
        items[0]["identity_before"]["capabilities"] = {}
        with self.assertRaises(ValueError):
            M.reduce_raw(raw(items))

    def test_publication_and_file_guards(self):
        with tempfile.TemporaryDirectory(prefix="wh2-cost-neutral-") as temporary:
            root = Path(temporary)
            path = root / "owned.json"
            M.publish(path, b"abc")
            self.assertEqual(path.stat().st_mode & 0o777, 0o400)
            self.assertEqual(M.read_regular(path, 3), b"abc")
            with self.assertRaises(FileExistsError):
                M.publish(path, b"replacement")
            with self.assertRaises(ValueError):
                M.read_regular(path, 2)
            link = root / "symlink"
            link.symlink_to(path)
            with self.assertRaises(OSError):
                M.read_regular(link, 3)
            hardlink = root / "hardlink"
            os.link(path, hardlink)
            with self.assertRaises(ValueError):
                M.read_regular(path, 3)
            self.assertEqual(M.read_regular(path, 3, installed=True), b"abc")
            fifo = root / "fifo"
            os.mkfifo(fifo)
            with self.assertRaises(ValueError):
                M.read_regular(fifo, 3)


class FakeChild:
    """Local pipes only; no subprocess is created by any lifecycle test."""
    def __init__(self, stdout=b"prefix\n", stderr=b""):
        self.pid = 4242
        self.returncode = None
        self.stdout, self.stderr = self.pipe(stdout), self.pipe(stderr)
        self.waited = False

    @staticmethod
    def pipe(data):
        read_fd, write_fd = os.pipe()
        os.write(write_fd, data)
        os.close(write_fd)
        return os.fdopen(read_fd, "rb")

    def poll(self):
        return self.returncode

    def wait(self, timeout=None):
        self.waited = True
        self.returncode = 0
        return 0

    def kill(self):
        self.returncode = -9


class LifecycleTests(unittest.TestCase):
    def test_outer_selector_failure_never_spawns_controller(self):
        with mock.patch.object(M.selectors, "DefaultSelector", side_effect=OSError("no selector")), \
             mock.patch.object(M.subprocess, "Popen") as spawn:
            with self.assertRaises(OSError):
                M.run(Path("/never-read"))
        spawn.assert_not_called()

    def test_prefix_retained_and_all_fds_closed_on_cleanup_errors(self):
        for failed_action in ("fsync", "fchmod"):
            with tempfile.TemporaryDirectory(prefix="wh2-cost-lifecycle-") as temporary:
                paths = [Path(temporary)/"raw", Path(temporary)/"err"]
                child = FakeChild(b"prefix\n", b"diagnostic\n")
                original_close = M.os.close
                with mock.patch.object(M.subprocess, "Popen", return_value=child), \
                     mock.patch.object(M.time, "monotonic", return_value=100), \
                     mock.patch.object(M.os, failed_action, side_effect=OSError("injected seal failure")), \
                     mock.patch.object(M.os, "close", wraps=original_close) as closed:
                    result = M.capture(Path("/never-executed"), "ab"*32, 110, paths)
                self.assertEqual(result[:3], (b"prefix\n", b"diagnostic\n", 0))
                self.assertIn("injected seal failure", result[3])
                self.assertEqual(closed.call_count, 2)
                self.assertTrue(child.waited and child.stdout.closed and child.stderr.closed)
                self.assertEqual(M.read_regular(paths[0], 100), result[0])

    def test_output_cap_keeps_exact_prefix_and_handles_exit_race(self):
        with tempfile.TemporaryDirectory(prefix="wh2-cost-cap-") as temporary:
            paths = [Path(temporary)/"raw", Path(temporary)/"err"]
            child = FakeChild(b"0123456789")
            with mock.patch.object(M.subprocess, "Popen", return_value=child), \
                 mock.patch.object(M.time, "monotonic", return_value=100), \
                 mock.patch.object(M, "RAW_CAP", 8), \
                 mock.patch.object(child, "kill", side_effect=ProcessLookupError):
                result = M.capture(Path("/never-executed"), "ab"*32, 110, paths)
            self.assertEqual(result[0], b"01234567")
            self.assertIn("output cap", result[3])
            self.assertTrue(child.waited)
            self.assertEqual(M.read_regular(paths[0], 8), result[0])

    def test_outer_observer_cleans_group_before_reaping_dead_leader(self):
        child = FakeChild(b"summary\n")
        stdout, stderr = mock.Mock(buffer=io.BytesIO()), mock.Mock(buffer=io.BytesIO())
        def kill_group(pid, sig):
            self.assertFalse(child.waited)
            self.assertEqual(pid, child.pid)
        with mock.patch.object(M.subprocess, "Popen", return_value=child), \
             mock.patch.object(M.time, "monotonic", return_value=100), \
             mock.patch.object(M.os, "waitid", return_value=object()) as observed, \
             mock.patch.object(M.os, "killpg", side_effect=kill_group) as killed, \
             mock.patch.object(M.sys, "stdout", stdout), mock.patch.object(M.sys, "stderr", stderr):
            M.run(Path("/never-read"))
        self.assertTrue(child.waited)
        self.assertEqual(killed.call_count, 1)
        self.assertTrue(observed.call_args[0][2] & os.WNOWAIT)
        self.assertEqual(stdout.buffer.getvalue(), b"summary\n")

    def test_late_reducer_cannot_publish_pass(self):
        with tempfile.TemporaryDirectory(prefix="wh2-cost-deadline-") as temporary:
            root = Path(temporary)
            output = root/"bundle"
            receipt = {"build": str(root)}
            receipt_path = root/"receipt.json"
            M.publish(receipt_path, M.canonical(receipt))
            now = [100]
            def capture(*args):
                for path, data in zip(args[3], (b"prefix\n", b"")):
                    M.publish(path, data)
                return b"prefix\n", b"", 0, None
            def late_reducer(data):
                now[0] = 121
                return dict(outcome="PASS", claim_sha256=M.sha(M.canonical(receipt)))
            with mock.patch.object(M, "OUTPUT", output), \
                 mock.patch.object(M.time, "monotonic", side_effect=lambda: now[0]), \
                 mock.patch.object(M, "current_receipt", return_value=receipt), \
                 mock.patch.object(M, "capture", side_effect=capture), \
                 mock.patch.object(M, "reduce_raw", side_effect=late_reducer), \
                 mock.patch("builtins.print"):
                M.controller_run(receipt_path)
            summary = M.decode(M.read_regular(output/"summary.json", 10000))
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertIn("deadline", summary["failure"])


if __name__ == "__main__":
    unittest.main()
