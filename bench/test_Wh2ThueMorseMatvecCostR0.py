#!/usr/bin/env python3
"""Neutral .73 schedule/ledger/controller tests; never selected codec work."""
import copy
import hashlib
import importlib.util
import io
import math
import os
from pathlib import Path
import struct
import tempfile
import time
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location(
    "matvec_cost", Path(__file__).with_name("Wh2ThueMorseMatvecCostR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

WIDTHS = (2, 64, 1280)
PAIRS = ((0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 1), (3, 1))
IDS = tuple(range(12)) + tuple(0xffffffff - 2*j for j in range(6))
# Literal two-pass order, independent of M.side().
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0)
ROWS = b"".join(bytes(int(c == (packet & 1023) % 6) for c in range(6)) for packet in IDS)
_lookup = bytearray(39936)
for _phase in range(2):
    for _id in range(1024):
        _lookup[_phase * 6144 + _id * 6 + _id % 6] = 1
for _offset in range(12288, 39936, 36):
    for _column in range(6):
        _lookup[_offset + 7 * _column] = 1
LOOKUP = bytes(_lookup)
LOOKUP_SHA = hashlib.sha256(LOOKUP).hexdigest()
COUNTERS = ("create_calls", "encode_calls", "free_calls", "decoder_create_calls",
            "feed_calls", "recover_calls", "decoder_free_calls")


def answer(status=0, required=0, written=0, kind=0):
    return dict(status=status, code=status, required=required, written=written, length_kind=kind)


def neutral_roster():
    index = 0
    for rep in range(12):
        for step in range(2):
            for ws in range(3):
                for ms in range(6):
                    for cs in range(7):
                        order = (rep + step) % 2
                        comparison = (2*rep + step + ws + ms + cs) % 7
                        for position, logical in enumerate(SIDES):
                            phase_pass = int(position >= 10)
                            j = rep % 4 if position < 2 else ((position - 2) % 8) // 2
                            phase = (rep + 6*phase_pass) % 12 + 12*j
                            yield dict(index=index, rep=rep, order=order,
                                       width=(rep + step + ws) % 3,
                                       metric=(rep + step + ws + ms) % 6,
                                       **{"class": comparison}, position=position,
                                       arm=PAIRS[comparison][logical ^ order],
                                       q=((2*phase + 1)*1000000)//96)
                            index += 1


def identity():
    return dict(derived=dict(family=26, model=8, stepping=1, full_apic_id=100,
        initial_apic_id_8=100, core_id=50, package_id=0, thread_id=0,
        threads_per_core=2, ccd_id=6, complex_id=6, logical_processors_per_package=128),
        before_cpu=50, after_cpu=50, requested_cpu=50, raw_capture_complete=True,
        semantic_validation_passed=True, singleton_affinity_verified=True,
        canonical_hex="0000", canonical_bytes=2, canonical_sha256=M.sha(bytes(2)),
        capabilities=dict.fromkeys(("leaf1_ecx", "leaf1_edx", "leaf6_eax", "leaf6_ecx",
            "leaf80000001_ecx", "leaf80000001_edx", "leaf80000008_ebx", "leaf80000021_eax",
            "max_basic_leaf", "max_extended_leaf"), 0))


def fixture():
    fixtures, handles, references = [], {}, {}
    feeds_total = 0
    for width, block in enumerate(WIDTHS):
        source = bytes((37*i + i//11) & 255 for i in range(6*block))
        corpus = b"".join(source[((packet & 1023) % 6)*block:
                                 (((packet & 1023) % 6)+1)*block] for packet in IDS)
        entries = []
        for rep in range(12):
            for step in range(4):
                arm = (rep + width + step) % 4
                for lane in range(2):
                    address = 0x1000000 + len(handles)*256
                    handles[(width, rep, lane, arm)] = address
                    traces = []
                    for pattern in range(2):
                        # Native cyclic far rows span three columns, then the
                        # fifth systematic reaches six. Controls are synthetic
                        # distinct profiles, with deliberately unequal lengths.
                        count = (6, 11)[pattern] if arm < 2 else (7, 9)[pattern] if arm == 2 else (8, 10)[pattern]
                        trace = dict(create=answer(), feeds=[answer(int(j+1 < count), block, 0,
                                     1 if arm < 2 else 2) for j in range(count)],
                                     recover=answer(0, 6*block, 6*block, 2 if arm == 2 else 1))
                        traces.append(trace)
                        references[(width, arm, pattern)] = trace
                        feeds_total += count
                    profile = b""
                    if arm == 3:
                        profile = b"WHV2" + struct.pack("<HHQQI", 1, 32, 0x4b295bbb47f4f9c9,
                                                       6*block, block) + bytes(4)
                    entries.append(dict(rep=rep, lane=lane, arm=arm, address=address,
                                        create=answer(), profile_hex=profile.hex(), decode=traces))
        fixtures.append(dict(width_index=width, block_bytes=block, message_bytes=6*block,
            output_stride=block+128, source_hex=source.hex(), source_sha256=M.sha(source),
            source_address=0x100000 + width*0x10000 + 64,
            output_addresses=[0x400000 + width*0x100000 + lane*0x10000 for lane in range(2)],
            rows_hex=ROWS.hex(), packets_hex=[corpus.hex()]*4,
            packets_sha256=[M.sha(corpus)]*4, handles=entries))
    runtime = dict(polynomial=0x14d, ssse3=1, avx2=1, gfni=1, avx512=1, address=0x700000)
    start, cpu = 10**12, 5*10**12
    prelude = [start+1000, cpu+1000, start+1200, start+1300, cpu+1400, start+1500, [0]*4, [0]*4]
    header = dict(type="header", protocol=M.PROTOCOL, schema=M.PROTOCOL+".raw.v1",
        claim_sha256="ab"*32, target_cpu=50, worker_start_ns=start, worker_start_cpu_ns=cpu,
        identity_before=identity(), runtime_before=runtime, clock_resolutions=[1, 1],
        lookup_sha256=LOOKUP_SHA, lookup_address=0x200000, fixtures=fixtures,
        prelude=dict(iterations=1 << 20, seed=0x9e3779b97f4a7c15,
                     final_state=0x43935dad1647741b, observation=prelude))
    now = start + 5000
    t0 = (now//2000000 + 2)*2000000
    records, total, ledger = [], 0, dict.fromkeys(COUNTERS, 0)
    for coordinates in neutral_roster():
        row = dict(coordinates, type="record", called=True, checked=True)
        target = t0 + 2000000*row["index"] + row["q"]
        ctarget = cpu + target - start
        duration = (20000, 16000, 24000, 22000)[row["arm"]]
        row.update(target=target, ready=target-200000,
            wait=[target-199900, ctarget-199900, target, ctarget],
            observation=[target+500, ctarget+700, target+1000, target+1000+duration,
                         ctarget+1500+duration, target+2000+duration, [0]*4, [0]*4])
        metric, arm, block = row["metric"], row["arm"], WIDTHS[row["width"]]
        create, decode = metric in (0, 3), metric >= 4
        count = 6 if metric in (1, 2) else 18 if metric == 3 else 0
        trace = references[(row["width"], arm, metric-4)] if decode else None
        calls = dict(create_calls=128 if create else 0, encode_calls=128*count,
            free_calls=128 if create else 0, decoder_create_calls=128 if decode else 0,
            feed_calls=128*len(trace["feeds"]) if decode else 0,
            recover_calls=128 if decode else 0, decoder_free_calls=128 if decode else 0)
        fresh = 0x3000000 + row["index"]*256
        addresses = dict(count=128, sha256=M.sha(struct.pack("<Q", fresh)*128), min=fresh, max=fresh) if create or decode else \
            dict(count=0, sha256=M.sha(b""), min=None, max=None)
        row["work"] = dict(calls, complete=True, create=answer() if create or decode else None,
            encode=[answer(0, block, block, 1) for _ in range(count)],
            feeds=trace["feeds"] if decode else [], recover=trace["recover"] if decode else None,
            handle=None if create else handles[(row["width"], row["rep"], row["order"], arm)],
            addresses=addresses)
        for key in COUNTERS:
            ledger[key] += calls[key]
        records.append(row)
        total += duration
    preflight = dict(encoder_creates=288, encoder_frees=288, encodes=5184,
                     decoder_creates=576, feeds=feeds_total, recovers=576, decoder_frees=576)
    footer = dict(ledger, type="footer", protocol=M.PROTOCOL, outcome="COMPLETE", failure=None,
        schedule_now_ns=now, t0=t0, records=54432, callbacks=54432, checked_callbacks=54432,
        sum_work_ns=total, sum_wait_wall_ns=54432*199900, sum_wait_cpu_ns=54432*199900,
        worker_end_ns=records[-1]["observation"][5]+1000,
        worker_end_cpu_ns=records[-1]["observation"][4]+1000,
        identity_after=identity(), runtime_after=dict(runtime), preflight=preflight)
    return header, records, footer, (handles, references, preflight)


class NeutralTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.header, cls.records, cls.footer, cls.preflight = fixture()

    def setUp(self):
        patches = [mock.patch.object(M, "lookup_bytes", return_value=LOOKUP),
                   mock.patch.object(M, "selected_rows", return_value=ROWS),
                   mock.patch.object(M, "LOOKUP_SHA", LOOKUP_SHA)]
        for name in ("write_lookup", "run", "controller_run", "capture"):
            patches.append(mock.patch.object(M, name, create=True,
                           side_effect=AssertionError("selected/native work forbidden")))
        patches.append(mock.patch.object(M.T, "lookup_bytes", side_effect=AssertionError("selected lookup forbidden")))
        patches.append(mock.patch.object(M.subprocess, "Popen", side_effect=AssertionError("child forbidden")))
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def reduce(self, header=None, records=None, footer=None):
        with mock.patch.object(M, "verify_preflight", return_value=self.preflight):
            return M.reduce_records(self.header if header is None else header,
                                    self.records if records is None else records,
                                    self.footer if footer is None else footer, "de"*32)

    def test_full_serialization_replay_counts_and_budget(self):
        started = time.monotonic()
        raw = b"".join(M.canonical(item) for item in [self.header] + self.records + [self.footer])
        self.assertLess(len(raw), 96 * 1024**2)
        analysis = M.reduce_raw(raw)
        self.assertLess(time.monotonic() - started, 135)
        self.assertEqual(analysis["outcome"], "PASS")
        self.assertEqual(analysis["callbacks"], 54432)
        self.assertEqual(analysis["work_calls"], {key: self.footer[key] for key in COUNTERS})
        self.assertEqual(analysis["raw_sha256"], M.sha(raw))
        self.assertEqual(len(analysis["statistics"]), 252)
        self.assertEqual(sum("control_pass" in item for item in analysis["statistics"]), 144)
        self.assertEqual(sum("treatment_pass" in item for item in analysis["statistics"]), 108)
        self.assertEqual(len(analysis["first_success"]), 24)
        self.assertEqual({row["packets"] for row in analysis["first_success"]}, {6, 7, 8, 9, 10, 11})
        for claim in ("all_K_claimed", "recovery_rate_claimed", "production_promotion_claimed"):
            self.assertIs(analysis[claim], False)

    def test_exact_roster_two_passes_and_chronology(self):
        self.assertEqual(list(M.roster()), list(neutral_roster()))
        self.assertEqual(len(self.records), 54432)
        visits = {}
        for start in range(0, len(self.records), 18):
            panel = self.records[start:start+18]
            first = panel[0]
            for p in range(2, 18, 2):
                left, right = panel[p:p+2]
                self.assertEqual(left["q"], right["q"])
                self.assertNotEqual(SIDES[p] ^ first["order"], SIDES[p+1] ^ first["order"])
                key = (first["width"], first["metric"], first["class"], first["order"], left["q"])
                visits.setdefault(key, []).append((first["rep"], SIDES[p] ^ first["order"]))
        self.assertEqual(len(visits), 3*6*7*2*48)
        for pair in visits.values():
            self.assertEqual(len(pair), 2)
            self.assertEqual((pair[0][0] - pair[1][0]) % 12, 6)
            self.assertNotEqual(pair[0][1], pair[1][1])
        self.assertEqual(sum(row["position"] < 2 for row in self.records), 6048)

    def test_preflight_real_verifier_neutral_payloads(self):
        handles, references, counts = M.verify_preflight(self.header)
        self.assertEqual((handles, references, counts), self.preflight)
        self.assertEqual(len(handles), 288)
        self.assertEqual(counts["feeds"], 4896)

    def test_independent_neutral_coefficient_oracles(self):
        cyclic = (1, 0, 0, 0, 0, 0)
        ids = IDS + tuple((1 << bit) + delta for bit in range(32) for delta in (-1, 0, 1))
        expected = b"".join(bytes(int(column == packet % 6) for column in range(6)) for packet in ids)
        self.assertEqual(M.coefficient_rows((cyclic, cyclic), ids), expected)

        def multiply(a, b):
            value = 0
            for bit in range(8):
                if b & (1 << bit):
                    value ^= a << bit
            for bit in range(14, 7, -1):
                if value & (1 << bit):
                    value ^= 0x14d << (bit-8)
            return value

        feedbacks = ((1, 1, 0, 1, 0, 1), (3, 1, 0, 1, 0, 1))
        matrices = []
        for feedback in feedbacks:
            matrix = [[0]*6 for _ in range(6)]
            for index in range(5):
                matrix[index+1][index] = 1
            for index in range(6):
                matrix[index][5] = feedback[index]
            matrices.append(matrix)
        product = [[int(r == c) for c in range(6)] for r in range(6)]
        chronological = bytearray()
        for packet in range(65):
            chronological.extend(product[row][0] for row in range(6))
            matrix = matrices[bin(packet).count("1") % 2]
            next_product = [[0]*6 for _ in range(6)]
            for r in range(6):
                for c in range(6):
                    for k in range(6):
                        next_product[r][c] ^= multiply(product[r][k], matrix[k][c])
            product = next_product
        self.assertEqual(M.coefficient_rows(feedbacks, range(65)), bytes(chronological))

    def test_preflight_malformed_payload_profile_status_and_first_count(self):
        header = copy.deepcopy(self.header)
        fixture = header["fixtures"][0]
        original = fixture["packets_hex"][0]
        fixture["packets_hex"][0] = original[:-2] + ("01" if original[-2:] != "01" else "02")
        fixture["packets_sha256"][0] = M.sha(bytes.fromhex(fixture["packets_hex"][0]))
        with self.assertRaises(ValueError):
            M.verify_preflight(header)
        fixture["packets_hex"][0] = original
        fixture["packets_sha256"][0] = M.sha(bytes.fromhex(original))
        entry = fixture["handles"][0]
        entry["rep"] = False
        with self.assertRaises(ValueError):
            M.verify_preflight(header)
        entry["rep"] = 0
        for field, value in (("code", True), ("status", 4), ("length_kind", 2), ("written", 1)):
            feeds = entry["decode"][0]["feeds"]
            original = feeds[0][field]
            feeds[0][field] = value
            with self.assertRaises(ValueError):
                M.verify_preflight(header)
            feeds[0][field] = original
        for feeds in ([], [answer(1, 2, 0, 1)]*13, [answer(0, 2, 0, 1)]*6):
            old = entry["decode"][0]["feeds"]
            entry["decode"][0]["feeds"] = feeds
            with self.assertRaises(ValueError):
                M.verify_preflight(header)
            entry["decode"][0]["feeds"] = old
        public = next(item for item in fixture["handles"] if item["arm"] == 3)
        public["profile_hex"] = "00" + public["profile_hex"][2:]
        with self.assertRaises(ValueError):
            M.verify_preflight(header)

    def test_work_all_metrics_and_raw_result_errors(self):
        handles, references, _ = self.preflight
        for metric in range(6):
            row = next(item for item in self.records if item["metric"] == metric)
            original = row["work"]
            self.assertEqual(M.verify_work(original, row, handles, references),
                             {key: original[key] for key in COUNTERS})
            cases = [("complete", 1), ("handle", False)] + [(key, True) for key in COUNTERS]
            for field, value in cases:
                altered = dict(original, **{field: value})
                with self.assertRaises(ValueError):
                    M.verify_work(altered, row, handles, references)
            for key in ("encode", "feeds"):
                if original[key]:
                    for field, value in (("status", False), ("code", 7), ("required", 999),
                                         ("written", 999), ("length_kind", 0)):
                        altered = copy.deepcopy(original)
                        altered[key][0][field] = value
                        with self.assertRaises(ValueError):
                            M.verify_work(altered, row, handles, references)
            if original["recover"] is not None:
                altered = copy.deepcopy(original)
                altered["recover"]["length_kind"] = 2 if row["arm"] != 2 else 1
                with self.assertRaises(ValueError):
                    M.verify_work(altered, row, handles, references)

    def test_roster_boolean_duplicate_missing_and_missed_slots(self):
        records = list(self.records)
        for field, value in (("index", False), ("order", False), ("arm", False),
                             ("called", 1), ("checked", 1), ("q", True)):
            records[0] = dict(self.records[0], **{field: value})
            with self.assertRaises(ValueError):
                self.reduce(records=records)
        records[0] = self.records[0]
        records[1] = records[0]
        with self.assertRaises(ValueError):
            self.reduce(records=records)
        with self.assertRaises(ValueError):
            self.reduce(records=self.records[:-1])
        for field, value in (("target", self.records[0]["target"]+2000000),
                             ("ready", self.records[0]["target"]-99999),
                             ("observation", []), ("wait", [1, 2])):
            records[0] = dict(self.records[0], **{field: value})
            records[1] = self.records[1]
            with self.assertRaises(ValueError):
                self.reduce(records=records)
        observation = list(self.records[0]["observation"])
        observation[2] = self.records[0]["target"]+5001
        records[0] = dict(self.records[0], observation=observation)
        with self.assertRaises(ValueError):
            self.reduce(records=records)

    def test_known_ratios_and_exact_treatment_gates(self):
        analysis = M.statistics(self.records)
        ratios = {"BB": 1, "CC": 1, "WW": 1, "PP": 1,
                  "BC": .8, "WC": 2/3, "PC": 8/11}
        for row in analysis["statistics"]:
            self.assertAlmostEqual(row["estimate"]["ratio"], ratios[row["comparison"]], places=12)
            if "treatment_pass" in row:
                strict = row["comparison"] == "WC" or (row["comparison"] == "BC" and
                    row["metric"] in ("distant_repair", "encoder_lifecycle", "distant_decoder_lifecycle"))
                self.assertEqual(row["upper_ratio_limit"], 1.0 if strict else 1.02)

    def altered_statistics(self, predicate, multiplier):
        changed = list(self.records)
        for index, row in enumerate(self.records):
            if row["position"] >= 2 and predicate(row) and (SIDES[row["position"]] ^ row["order"]):
                observation = list(row["observation"])
                observation[3] = observation[2] + int((observation[3]-observation[2])*multiplier)
                changed[index] = dict(row, observation=observation)
        return M.statistics(changed)

    def test_control_fail_priority_and_outlier_retention(self):
        analysis = self.altered_statistics(lambda row: row["width"] == 0 and row["order"] == 0 and
            row["metric"] == 0 and row["class"] in (0, 4) and row["rep"] == 11, 4)
        self.assertEqual(analysis["outcome"], "CONTROL_FAIL")
        self.assertIn([0, 0, 0, 0], analysis["failed_controls"])
        self.assertIn([0, 0, 4, 0], analysis["failed_treatments"])
        control = next(row for row in analysis["statistics"] if row["block_bytes"] == 2 and
            row["metric"] == "encoder_create_free" and row["comparison"] == "BB" and row["order"] == 0)
        logs = control["replicate_logs"]
        self.assertEqual(len(logs), 12)
        self.assertEqual(logs[:11], [0]*11)
        self.assertAlmostEqual(logs[11], math.log(4), places=12)
        mean = math.log(4)/12
        radius = 2.200985160082949*math.sqrt(sum((x-mean)**2 for x in logs)/11/12)
        self.assertAlmostEqual(control["estimate"]["upper95_log"], mean+radius, places=12)

    def test_gain_required_and_unaffected_two_percent_bound(self):
        equal = self.altered_statistics(lambda row: row["class"] == 4, 1.25)
        self.assertEqual(equal["outcome"], "FAIL")
        self.assertEqual(len(equal["failed_treatments"]), 3*3*2)
        self.assertTrue(all(item[1] in (2, 3, 5) for item in equal["failed_treatments"]))
        over = self.altered_statistics(lambda row: row["class"] == 4 and row["metric"] == 1, 1.2875)
        self.assertEqual(len(over["failed_treatments"]), 3*2)
        self.assertEqual(over["failed_controls"], [])

    def test_footer_budgets_and_count_types(self):
        wall_span = self.footer["worker_end_ns"]-self.header["worker_start_ns"]
        for field, value in (("callbacks", True), ("sum_work_ns", self.footer["sum_work_ns"]+1),
                             ("worker_end_ns", self.header["worker_start_ns"]+125*10**9),
                             ("worker_end_cpu_ns", self.header["worker_start_cpu_ns"]+121*10**9),
                             ("worker_end_cpu_ns", self.header["worker_start_cpu_ns"]+wall_span+1)):
            with self.assertRaises(ValueError):
                self.reduce(footer=dict(self.footer, **{field: value}))
        self.assertLess(self.footer["sum_work_ns"], 8*10**9)
        self.assertLess(self.footer["worker_end_ns"]-self.header["worker_start_ns"], 125*10**9)
        with mock.patch.object(M, "WORK_CAP", self.footer["sum_work_ns"]-1):
            with self.assertRaises(ValueError):
                self.reduce()
        with self.assertRaises(ValueError):
            M.time_left(time.monotonic()-1)
        observation = list(self.records[0]["observation"])
        observation[4] = observation[1]+observation[5]-observation[0]+1000
        with self.assertRaises(ValueError):
            M.clocks(observation)

    def test_raw_framing_and_duplicate_json_keys(self):
        for raw in (b"", b"{}", b"{}\n", b"{}\n"*54433):
            with self.assertRaises(ValueError):
                M.reduce_raw(raw)
        with self.assertRaises(ValueError):
            M.decode(b'{"index":0,"index":1}')

    def test_extra_keys_fail_closed_at_each_container(self):
        for name in ("header", "footer", "prelude", "record", "identity_before", "identity_after",
                     "runtime_before", "runtime_after"):
            header, footer, records = dict(self.header), dict(self.footer), list(self.records)
            if name == "header":
                header["unexpected"] = 1
            elif name == "footer":
                footer["unexpected"] = 1
            elif name == "prelude":
                header["prelude"] = dict(header["prelude"], unexpected=1)
            elif name.endswith("before"):
                header[name] = dict(header[name], unexpected=1)
            elif name.endswith("after"):
                footer[name] = dict(footer[name], unexpected=1)
            else:
                records[0] = dict(records[0], unexpected=1)
            with self.assertRaises(ValueError):
                self.reduce(header=header, records=records, footer=footer)
        for name in ("work", "addresses"):
            row = self.records[0]
            work = copy.deepcopy(row["work"])
            (work if name == "work" else work["addresses"])["unexpected"] = 1
            with self.assertRaises(ValueError):
                M.verify_work(work, row, self.preflight[0], self.preflight[1])
        for name in ("fixture", "handle", "decode"):
            header = copy.deepcopy(self.header)
            value = header["fixtures"][0]
            if name != "fixture":
                value = value["handles"][0]
            if name == "decode":
                value = value["decode"][0]
            value["unexpected"] = 1
            with self.assertRaises(ValueError):
                M.verify_preflight(header)


class FakeChild:
    """Tiny local pipes, never an actual subprocess or native executable."""
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


class LifecycleTests(unittest.TestCase):
    def test_selector_failure_never_spawns(self):
        with mock.patch.object(M.selectors, "DefaultSelector", side_effect=OSError("selector failure")), \
             mock.patch.object(M.subprocess, "Popen") as spawn:
            with self.assertRaises(OSError):
                M.run(Path("/never-read"))
            with self.assertRaises(OSError):
                M.capture(Path("/never-run"), "ab"*32, 110, ())
        spawn.assert_not_called()

    def test_outer_owns_group_and_kills_before_reaping(self):
        child = FakeChild(b"summary\n")
        stdout, stderr = mock.Mock(buffer=io.BytesIO()), mock.Mock(buffer=io.BytesIO())

        def kill_group(pid, sig):
            self.assertFalse(child.waited)
            self.assertEqual(pid, child.pid)
            self.assertEqual(sig, M.signal.SIGKILL)

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

    def test_capture_caps_sealing_and_deadline_preserve_prefix(self):
        for failure in ("stdout", "stderr", "fsync", "fchmod", "deadline"):
            with tempfile.TemporaryDirectory(prefix="wh2-matvec-capture-") as temporary:
                paths = [Path(temporary)/"raw", Path(temporary)/"err"]
                child = FakeChild(b"0123456789", b"abcdefghij" if failure == "stderr" else b"")
                patch = mock.patch.object(M, "RAW_CAP" if failure == "stdout" else "ERR_CAP", 8) \
                    if failure in ("stdout", "stderr") else \
                    mock.patch.object(M.os, failure, side_effect=OSError("seal failure")) \
                    if failure in ("fsync", "fchmod") else mock.patch.object(M, "WORK_CAP", M.WORK_CAP)
                with mock.patch.object(M.subprocess, "Popen", return_value=child), \
                     mock.patch.object(M.time, "monotonic", return_value=111 if failure == "deadline" else 100), patch:
                    result = M.capture(Path("/never-executed"), "ab"*32, 110, paths)
                self.assertIsNotNone(result[3])
                self.assertTrue(child.waited and child.stdout.closed and child.stderr.closed)
                for path, data in zip(paths, result[:2]):
                    self.assertEqual(M.read_regular(path, 100), data)
                if failure == "stdout":
                    self.assertEqual(result[0], b"01234567")
                elif failure == "stderr":
                    self.assertEqual(result[1], b"abcdefgh")
                elif failure == "deadline":
                    self.assertEqual(result[:2], (b"", b""))
                else:
                    self.assertEqual(result[0], b"0123456789")

    def test_capture_never_overwrites_prior_spool(self):
        with tempfile.TemporaryDirectory(prefix="wh2-matvec-existing-") as temporary:
            paths = [Path(temporary)/"raw", Path(temporary)/"err"]
            M.publish(paths[0], b"original\n")
            with mock.patch.object(M.subprocess, "Popen") as spawn:
                result = M.capture(Path("/never-executed"), "ab"*32, 110, paths)
            spawn.assert_not_called()
            self.assertIsNotNone(result[3])
            self.assertEqual(M.read_regular(paths[0], 100), b"original\n")

    def test_prelaunch_mismatch_never_claims(self):
        with tempfile.TemporaryDirectory(prefix="wh2-matvec-prelaunch-") as temporary:
            root = Path(temporary)
            M.publish(root/"receipt", M.canonical(dict(build=str(root))))
            with mock.patch.object(M, "OUTPUT", root/"bundle"), \
                 mock.patch.object(M, "current_receipt", return_value={}), \
                 mock.patch.object(M, "capture") as capture:
                with self.assertRaises(ValueError):
                    M.controller_run(root/"receipt")
            capture.assert_not_called()
            self.assertFalse((root/"bundle").exists())

    def test_invalid_prefix_postpins_and_late_analysis_sealed(self):
        for failure in ("capture", "postpins", "late", "before-spool"):
            with tempfile.TemporaryDirectory(prefix="wh2-matvec-controller-") as temporary:
                root = Path(temporary)
                output, receipt = root/"bundle", dict(build=str(root))
                M.publish(root/"receipt", M.canonical(receipt))
                now = [100]

                def capture(executable, claim, deadline, paths):
                    self.assertEqual(executable, root/M.TARGET)
                    self.assertEqual(claim, M.sha(M.canonical(receipt)))
                    self.assertEqual(deadline, 235)
                    self.assertEqual(list(paths), [output/"raw.jsonl", output/"stderr.txt"])
                    if failure == "before-spool":
                        raise OSError("injected before spool")
                    for path, data in zip(paths, (b"prefix\n", b"")):
                        M.publish(path, data)
                    return b"prefix\n", b"", 0, "capture failure" if failure == "capture" else None

                def reduce(data):
                    if failure == "late":
                        now[0] = 236
                    return dict(outcome="PASS", claim_sha256=M.sha(M.canonical(receipt)))

                with mock.patch.object(M, "OUTPUT", output), \
                     mock.patch.object(M.time, "monotonic", side_effect=lambda: now[0]), \
                     mock.patch.object(M, "current_receipt", side_effect=[receipt, {} if failure == "postpins" else receipt]), \
                     mock.patch.object(M, "capture", side_effect=capture), \
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
