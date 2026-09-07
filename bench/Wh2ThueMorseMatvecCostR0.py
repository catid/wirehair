#!/usr/bin/env python3
"""One-shot four-arm .73 full-cost screen; replay never executes a codec."""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import selectors
import shlex
import signal
import struct
import subprocess
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Pure utility reuse only: never mutate an older module's namespace or caps.
A = sibling("_matvec_cost_common", "Wh2AlignedIntermediateCostR0.py")
T = sibling("_matvec_cost_lookup", "Wh2ThueMorseTinyPayloadCostR0.py")
require, exact, integer = A.require, A.exact, A.integer
canonical, decode, sha, hex_bytes = A.canonical, A.decode, A.sha, A.hex_bytes
read_regular, pin, checked, publish = A.read_regular, A.pin, A.checked, A.publish
time_left = A.time_left
ROOT = A.ROOT
PROTOCOL = "wirehair.wh2.thue-matvec-cost-r0"
OUTPUT = Path("/var/tmp/wh2-thue-matvec-cost-r0")
TARGET = "wh2_thue_matvec_cost_r0"
RAW_CAP, ERR_CAP, TOTAL_CAP = 96 * 1024**2, 65536, 256 * 1024**2
WORKER_WALL, WORKER_CPU, CONTROLLER_WALL, OBSERVER_WALL = 125, 120, 135, 145
WORK_CAP = 8 * 10**9
BATCH = 128
LOOKUP_SHA = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d"
WIDTHS = (2, 64, 1280)
METRICS = ("encoder_create_free", "low_repair", "distant_repair", "encoder_lifecycle",
           "low_decoder_lifecycle", "distant_decoder_lifecycle")
PAIRS = ((0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 1), (3, 1))
COMPARISONS = ("BB", "CC", "WW", "PP", "BC", "WC", "PC")
PACKET_IDS = tuple(range(12)) + tuple(0xffffffff - 2*j for j in range(6))
CALLBACKS = 12 * 2 * 3 * 6 * 7 * 18
GENERATED_SHA = "2bc9b40321a2a95af982d66c7e9b3623f57674ba80ca3faf3494d9a54e9f84a8"
NEW_SOURCES = tuple("bench/" + name for name in (
    "Wh2ThueMorseMatvecCostR0.py", "test_Wh2ThueMorseMatvecCostR0.py",
    "Wh2ThueMorseMatvecCostR0.cpp", "Wh2ThueMorseMatvecCostR0Bridge.h",
    "Wh2ThueMorseMatvecCostR0Bridge.cpp", "Wh2ThueMorseMatvecCostR0BridgeTest.cpp",
    "Wh2ThueMorseMatvecCostR0/CMakeLists.txt"))


def lookup_bytes():
    return T.lookup_bytes()


def write_lookup(path):
    packed = lookup_bytes()
    path = Path(path)
    if path.exists() or path.is_symlink():
        require(read_regular(path, 39936) == packed, "existing lookup differs; use a fresh build")
    else:
        publish(path, packed)


def side(position, order):
    if position < 2:
        return position ^ order
    phase_pass, within = divmod(position - 2, 8)
    return (0, 1, 1, 0)[within // 2] ^ phase_pass ^ (within % 2) ^ order


def roster():
    index = 0
    for rep in range(12):
        for step in range(2):
            order = (rep + step) % 2
            for ws in range(3):
                width = (rep + step + ws) % 3
                for ms in range(6):
                    metric = (rep + step + ws + ms) % 6
                    for cs in range(7):
                        comparison = (2*rep + step + ws + ms + cs) % 7
                        for position in range(18):
                            phase_pass = 0 if position < 2 else (position - 2) // 8
                            j = rep % 4 if position < 2 else ((position - 2) % 8) // 2
                            phase = (rep + 6*phase_pass) % 12 + 12*j
                            yield dict(index=index, rep=rep, order=order, width=width,
                                       metric=metric, **{"class": comparison}, position=position,
                                       arm=PAIRS[comparison][side(position, order)],
                                       q=((2*phase + 1)*1000000)//96)
                            index += 1


def statistics(records):
    require(len(records) == CALLBACKS, "all statistical records")
    groups = {}
    for start in range(0, len(records), 18):
        panel = records[start:start + 18]
        first, contrasts = panel[0], []
        for j in range(8):
            by_side = {}
            for row in panel[2 + 2*j:4 + 2*j]:
                key = side(row["position"], row["order"])
                require(key not in by_side, "duplicate logical side")
                by_side[key] = row["observation"][3] - row["observation"][2]
            require(set(by_side) == {0, 1} and min(by_side.values()) > 0, "paired positive times")
            contrasts.append(math.log(by_side[1]) - math.log(by_side[0]))
        key = (first["width"], first["metric"], first["class"], first["order"])
        groups.setdefault(key, []).append(math.fsum(contrasts) / 8)
    require(len(groups) == 252, "all separate width/metric/class/order groups")
    results, controls, treatments = [], [], []
    bound = math.log1p(.02)
    for (width, metric, comparison, order), values in sorted(groups.items()):
        estimate = A.confidence(values)
        result = dict(block_bytes=WIDTHS[width], metric=METRICS[metric],
                      comparison=COMPARISONS[comparison], order=order,
                      estimate=estimate, replicate_logs=values)
        if comparison < 4:
            result["control_pass"] = estimate["lower95_log"] > -bound and estimate["upper95_log"] < bound
            if not result["control_pass"]:
                controls.append([width, metric, comparison, order])
        else:
            # The Map-negative paths need equivalence, not an invented gain.
            # Every direct WH1 comparison must be a resolved win; no 5% floor.
            limit = 0.0 if comparison == 5 or (comparison == 4 and metric in (2, 3, 5)) else bound
            result["treatment_pass"] = estimate["upper95_log"] < limit
            result["upper_ratio_limit"] = math.exp(limit)
            if not result["treatment_pass"]:
                treatments.append([width, metric, comparison, order])
        results.append(result)
    return dict(outcome="CONTROL_FAIL" if controls else "FAIL" if treatments else "PASS",
                statistics=results, failed_controls=controls, failed_treatments=treatments,
                WH1_compared=True, public_WH2_compared=True, all_K_claimed=False,
                recovery_rate_claimed=False, production_promotion_claimed=False)


def coefficient_rows(feedbacks, packet_ids):
    """Independent polynomial arithmetic and binary prefix decomposition."""
    multiply = T.multiply

    def product(a, b):
        out = [0] * 36
        for r in range(6):
            for c in range(6):
                for k in range(6):
                    out[6*r+c] ^= multiply(a[6*r+k], b[6*k+c])
        return out

    pair = []
    for feedback in feedbacks:
        matrix = [0] * 36
        for i in range(5):
            matrix[6*(i+1)+i] = 1
        for i, value in enumerate(feedback):
            matrix[6*i+5] = value
        pair.append(matrix)
    levels = [pair]
    for _ in range(31):
        a, b = levels[-1]
        levels.append([product(a, b), product(b, a)])
    rows = bytearray()
    for packet_id in packet_ids:
        vector = [1, 0, 0, 0, 0, 0]
        for bit in range(32):
            if packet_id & (1 << bit):
                phase = bin(packet_id >> (bit + 1)).count("1") % 2
                matrix, out = levels[bit][phase], [0] * 6
                for r in range(6):
                    for c in range(6):
                        out[r] ^= multiply(matrix[6*r+c], vector[c])
                vector = out
        rows.extend(vector)
    return bytes(rows)


def selected_rows():
    return coefficient_rows(((124, 127, 152, 84, 241, 63),
                             (125, 127, 152, 84, 241, 63)), PACKET_IDS)


def result(status=0, required=0, written=0, kind=0):
    return dict(status=status, code=status, required=required, written=written, length_kind=kind)


def fields(value, names):
    require(type(value) is dict and set(value) == set(names.split()), "exact object fields: " + names)


def clocks(observation, previous=None):
    duration = A.clocks(observation, previous)
    require(observation[4]-observation[1] <= observation[5]-observation[0],
            "thread CPU cannot exceed its outer monotonic envelope")
    return duration


def verify_preflight(header):
    exact(header["lookup_sha256"], LOOKUP_SHA, "lookup identity")
    lookup_address = integer(header["lookup_address"], 1, (1 << 64) - 39936)
    lookup = lookup_bytes()
    rows = selected_rows()
    require(len(lookup) == 39936 and sha(lookup) == LOOKUP_SHA, "selected lookup exact bytes")
    require(len(rows) == 108 and rows[:72] == lookup[:72], "independent systematic/low rows")
    fixtures = header["fixtures"]
    require(type(fixtures) is list and len(fixtures) == 3, "three frozen widths")
    handles, references = {}, {}
    spans, seen_handles = [(lookup_address, lookup_address + 39936)], set()
    feed_total = 0
    for width, (B, fixture) in enumerate(zip(WIDTHS, fixtures)):
        fields(fixture, "width_index block_bytes message_bytes output_stride source_hex source_sha256 "
               "source_address output_addresses rows_hex packets_hex packets_sha256 handles")
        exact([fixture["width_index"], fixture["block_bytes"], fixture["message_bytes"], fixture["output_stride"]],
              [width, B, 6*B, B+128], "fixture dimensions")
        source = bytes((37*i + i//11) & 255 for i in range(6*B))
        exact(fixture["source_hex"], source.hex(), "source formula")
        exact(fixture["source_sha256"], sha(source), "source hash")
        address = integer(fixture["source_address"], 64, (1 << 64) - (6*B+64))
        require(address % 4096 == 64, "guarded source alignment")
        spans.append((address-64, address+6*B+64))
        addresses = fixture["output_addresses"]
        require(type(addresses) is list and len(addresses) == 2, "two output lanes")
        for address in addresses:
            integer(address, 1, (1 << 64) - 18*(B+128))
            require(address % 4096 == 0, "output base alignment")
            spans.append((address, address + 18*(B+128)))
        exact(fixture["rows_hex"], rows.hex(), "independent distant/low coefficient oracle")
        require(type(fixture["packets_hex"]) is list and len(fixture["packets_hex"]) == 4,
                "four distinct codec packet corpora")
        packets = [hex_bytes(value, 18*B) for value in fixture["packets_hex"]]
        exact(fixture["packets_sha256"], [sha(value) for value in packets], "packet hashes")
        require(all(value[:6*B] == source for value in packets), "all-arm systematic oracle")
        oracle = bytearray()
        for packet in range(18):
            for byte in range(B):
                value = 0
                for column in range(6):
                    value ^= T.multiply(rows[6*packet+column], source[column*B+byte])
                oracle.append(value)
        require(packets[0] == packets[1] == bytes(oracle), "independent GF256 payload oracle")
        entries = fixture["handles"]
        require(type(entries) is list and len(entries) == 96, "retained handle count")
        profiles = {}
        expected_order = [(rep, lane, (rep + width + step) % 4)
                          for rep in range(12) for step in range(4) for lane in range(2)]
        for entry, (rep, lane, arm) in zip(entries, expected_order):
            fields(entry, "rep lane arm address create profile_hex decode")
            exact([entry["rep"], entry["lane"], entry["arm"]], [rep, lane, arm], "allocation order")
            # Exact native sizes are known; public pointers are opaque. Do not
            # fabricate a public object span from the bridge's POD size.
            extent = 96 if arm < 2 else 1
            address = integer(entry["address"], 1, (1 << 64) - extent)
            require(address not in seen_handles, "retained handles cannot alias")
            seen_handles.add(address)
            spans.append((address, address + extent))
            handles[(width, rep, lane, arm)] = address
            exact(entry["create"], result(), "successful real encoder create")
            profile = hex_bytes(entry["profile_hex"], 32 if arm == 3 else 0)
            if arm == 3:
                require(profile[:28] == b"WHV2" + struct.pack("<HHQQI", 1, 32,
                        0x4b295bbb47f4f9c9, 6*B, B) and profile[29:] == bytes(3),
                        "current certified profile, exact M/B and reserved bytes")
            if arm in profiles:
                exact(profile, profiles[arm], "same emitted profile across physical handles")
            else:
                profiles[arm] = profile
            traces = entry["decode"]
            require(type(traces) is list and len(traces) == 2, "low and distant fresh decode")
            for pattern, trace in enumerate(traces):
                fields(trace, "create feeds recover")
                exact(trace["create"], result(), "fresh decoder creation")
                feeds = trace["feeds"]
                require(type(feeds) is list and 6 <= len(feeds) <= 12, "fixed trace first success cap")
                expected = [result(1 if j+1 < len(feeds) else 0, B, 0, 1 if arm < 2 else 2)
                            for j in range(len(feeds))]
                exact(feeds, expected, "honest first-success statuses and length provenance")
                exact(trace["recover"], result(0, 6*B, 6*B, 2 if arm == 2 else 1),
                      "full recovered message result")
                key = (width, arm, pattern)
                if key in references:
                    exact(trace, references[key], "same actual decode prefix across handles")
                else:
                    references[key] = trace
                feed_total += len(feeds)
        for pattern in range(2):
            exact(references[(width, 0, pattern)], references[(width, 1, pattern)],
                  "GFNI changes arithmetic only, not decoder behavior")
    for i, (begin, end) in enumerate(spans):
        require(all(end <= other_begin or other_end <= begin for other_begin, other_end in spans[:i]),
                "input/output/known handle spans must be disjoint")
    preflight = dict(encoder_creates=288, encoder_frees=288, encodes=5184,
                     decoder_creates=576, feeds=feed_total, recovers=576, decoder_frees=576)
    return handles, references, preflight


def verify_work(work, coordinates, handles, references):
    fields(work, "create_calls encode_calls free_calls decoder_create_calls feed_calls recover_calls "
           "decoder_free_calls create encode feeds recover complete handle addresses")
    c = coordinates
    metric, arm, B = c["metric"], c["arm"], WIDTHS[c["width"]]
    decoding = metric >= 4
    creating = metric in (0, 3)
    count = 6 if metric in (1, 2) else 18 if metric == 3 else 0
    trace = references[(c["width"], arm, metric-4)] if decoding else None
    feeds = len(trace["feeds"]) if decoding else 0
    expected_calls = dict(create_calls=BATCH if creating else 0, encode_calls=BATCH*count,
                          free_calls=BATCH if creating else 0,
                          decoder_create_calls=BATCH if decoding else 0, feed_calls=BATCH*feeds,
                          recover_calls=BATCH if decoding else 0, decoder_free_calls=BATCH if decoding else 0)
    for key, value in expected_calls.items():
        exact(work[key], value, "exact work accounting: " + key)
    exact(work["create"], result() if creating or decoding else None, "fresh create result")
    exact(work["encode"], [result(0, B, B, 1)] * count, "all encoded result checks")
    exact(work["feeds"], trace["feeds"] if decoding else [], "all actual feed result checks")
    exact(work["recover"], trace["recover"] if decoding else None, "actual recovery result")
    exact(work["complete"], True, "complete workload")
    exact(work["handle"], None if creating else handles[(c["width"], c["rep"], c["order"], arm)],
          "real persistent encoder descriptor/handle")
    addresses = work["addresses"]
    fields(addresses, "count sha256 min max")
    if creating or decoding:
        exact(addresses["count"], BATCH, "every fresh allocation address recorded")
        hex_bytes(addresses["sha256"], 32)
        require(integer(addresses["min"], 1) <= integer(addresses["max"], 1), "fresh address range")
    else:
        exact(addresses, dict(count=0, sha256=sha(b""), min=None, max=None), "no new repair handles")
    return expected_calls


def reduce_raw(raw):
    require(type(raw) is bytes and 0 < len(raw) <= RAW_CAP and raw.endswith(b"\n"), "bounded complete JSONL")
    lines = raw.splitlines()
    require(len(lines) == CALLBACKS+2, "exact header/record/footer count")
    header, *records, footer = [decode(line) for line in lines]
    return reduce_records(header, records, footer, sha(raw))


def reduce_records(header, records, footer, raw_sha):
    fields(header, "type protocol schema claim_sha256 target_cpu worker_start_ns worker_start_cpu_ns "
           "identity_before runtime_before clock_resolutions lookup_address lookup_sha256 fixtures prelude")
    fields(footer, "type protocol outcome failure schedule_now_ns t0 records callbacks checked_callbacks "
           "create_calls encode_calls free_calls decoder_create_calls feed_calls recover_calls decoder_free_calls "
           "sum_work_ns sum_wait_wall_ns sum_wait_cpu_ns worker_end_ns worker_end_cpu_ns identity_after runtime_after preflight")
    exact([header["type"], footer["type"], header["protocol"], footer["protocol"], header["schema"]],
          ["header", "footer", PROTOCOL, PROTOCOL, PROTOCOL+".raw.v1"], "framing")
    exact([footer["outcome"], footer["failure"]], ["COMPLETE", None], "complete correctness")
    exact(header["target_cpu"], 50, "target CPU")
    hex_bytes(header["claim_sha256"], 32)
    for value in (header["identity_before"], footer["identity_after"]):
        fields(value, "after_cpu before_cpu canonical_bytes canonical_hex canonical_sha256 capabilities "
               "derived raw_capture_complete requested_cpu semantic_validation_passed singleton_affinity_verified")
    for value in (header["runtime_before"], footer["runtime_after"]):
        fields(value, "polynomial address ssse3 avx2 gfni avx512")
    A.verify_identity(header["identity_before"])
    A.verify_identity(footer["identity_after"])
    exact(header["identity_before"], footer["identity_after"], "stable target identity")
    A.verify_runtime(header["runtime_before"])
    A.verify_runtime(footer["runtime_after"])
    exact(header["runtime_before"], footer["runtime_after"], "stable GF runtime")
    exact(header["clock_resolutions"], [1, 1], "clock resolutions")
    handles, references, preflight_counts = verify_preflight(header)
    prelude = header["prelude"]
    fields(prelude, "iterations seed final_state observation")
    exact([prelude["iterations"], prelude["seed"], prelude["final_state"]],
          [1 << 20, 0x9e3779b97f4a7c15, 0x43935dad1647741b], "single fixed prelude")
    clocks(prelude["observation"])
    start, cpu_start = integer(header["worker_start_ns"]), integer(header["worker_start_cpu_ns"])
    end, cpu_end = integer(footer["worker_end_ns"]), integer(footer["worker_end_cpu_ns"])
    require(start <= prelude["observation"][0] and cpu_start <= prelude["observation"][1], "prelude lifetime")
    require(0 < end-start < WORKER_WALL*10**9 and 0 <= cpu_end-cpu_start <= WORKER_CPU*10**9, "worker budgets")
    require(cpu_end-cpu_start <= end-start, "thread CPU cannot exceed worker wall time")
    now, t0 = integer(footer["schedule_now_ns"]), integer(footer["t0"])
    require(now >= prelude["observation"][5] and t0 == (now//2000000+2)*2000000, "schedule epoch")
    previous = prelude["observation"]
    work_sum = wait_wall = wait_cpu = 0
    ledger = dict.fromkeys(("create_calls", "encode_calls", "free_calls", "decoder_create_calls",
                           "feed_calls", "recover_calls", "decoder_free_calls"), 0)
    require(type(records) is list and len(records) == CALLBACKS, "complete record roster")
    for expected, row in zip(roster(), records):
        fields(row, "type index rep order width metric class position arm q target ready wait observation called checked work")
        exact({key: row[key] for key in expected}, expected, "exact callback roster")
        exact([row["type"], row["called"], row["checked"]], ["record", True, True], "checked callback")
        target = integer(row["target"])
        require(target == t0 + 2000000*row["index"] + row["q"], "phase target")
        ready = integer(row["ready"])
        require(previous[5] <= ready <= target-100000, "preparation margin")
        duration = clocks(row["observation"], previous)
        wait = row["wait"]
        require(type(wait) is list and len(wait) == 4, "wait shape")
        wm0, wc0, wm1, wc1 = [integer(x) for x in wait]
        require(ready <= wm0 <= target <= wm1 <= row["observation"][0] and
                previous[4] <= wc0 <= wc1 <= row["observation"][1], "wait ordering")
        require(target <= row["observation"][2] <= target+5000, "start lateness")
        work_sum += duration
        wait_wall += wm1-wm0
        wait_cpu += wc1-wc0
        previous = row["observation"]
        counts = verify_work(row["work"], row, handles, references)
        for key, value in counts.items():
            ledger[key] += value
    require(previous[5] <= end and previous[4] <= cpu_end and work_sum <= WORK_CAP, "final lifetime/work cap")
    totals = dict(records=CALLBACKS, callbacks=CALLBACKS, checked_callbacks=CALLBACKS,
                  sum_work_ns=work_sum, sum_wait_wall_ns=wait_wall, sum_wait_cpu_ns=wait_cpu, **ledger)
    for key, value in totals.items():
        exact(footer[key], value, "footer " + key)
    exact(footer["preflight"], preflight_counts, "complete untimed API ledger")
    require(wait_wall <= end-start and wait_cpu <= cpu_end-cpu_start, "conditioning included in lifetime")
    analysis = statistics(records)
    analysis.update(protocol=PROTOCOL, claim_sha256=header["claim_sha256"], raw_sha256=raw_sha,
                    callbacks=CALLBACKS, batch=BATCH, work_calls=ledger, sum_work_ns=work_sum,
                    conditioning=dict(prelude_wall_ns=prelude["observation"][3]-prelude["observation"][2],
                                      wait_wall_ns=wait_wall, wait_cpu_ns=wait_cpu,
                                      worker_wall_ns=end-start, worker_cpu_ns=cpu_end-cpu_start),
                    first_success=[dict(block_bytes=WIDTHS[w], arm=a, pattern=p, packets=len(v["feeds"]))
                                   for (w, a, p), v in sorted(references.items())])
    return analysis


def current_receipt(build, deadline):
    """Pin the fresh build and its actual dependencies before and after WORK."""
    build = Path(build).resolve(strict=True)
    require(build.is_dir() and ROOT not in build.parents and build != ROOT, "external build")
    require(checked(["ninja", "-C", build, "-n"], deadline).decode().strip().endswith(
        "ninja: no work to do."), "current build required")
    cache = read_regular(build / "CMakeCache.txt", 1024**2).decode()
    for name in ("WH2_MATVEC_COST_SANITIZERS", "BUILD_SHARED_LIBS", "BUILD_TESTS", "BUILD_CODEC_V2",
                 "MARCH_NATIVE", "WIREHAIR_BUILD_BOTH", "WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0"):
        require(name + ":BOOL=OFF" in cache, "native isolated build option: " + name)
    require(sha(read_regular(build / "lookup.bin", 39936)) == LOOKUP_SHA, "build lookup")
    generated = build / "codec-neutral/Wh2ThueMorseMatvecCodecR0.generated.cpp"
    require(sha(read_regular(generated, 65536)) == GENERATED_SHA, "exact Apply-only generated codec")
    commands = decode(read_regular(build / "compile_commands.json", 2 * 1024**2))
    require(type(commands) is list, "compile commands")
    public_files = ("wirehair.cpp", "gf256.cpp", "WirehairCodec.cpp", "WirehairTools.cpp") + tuple(
        "codec/WirehairV2" + name + ".cpp" for name in
        ("Codec", "Peel", "Plan", "Policy", "Precode", "PrecodeDecode", "PrecodeEncode", "Profile", "Seeds", "Solve"))
    support = ("Wh2ThueMorseMatvecGfniR0.cpp", "Wh2ThueMorseNativeCodec.cpp",
               "Wh2ThueMorseMatvecCodecR0Test.cpp", "Wh2ThueMorseMatvecCodecR0Test.cpp",
               "Wh2ThueMorseMatvecCostR0Bridge.cpp", "Wh2ThueMorseMatvecCostR0Bridge.cpp",
               "Wh2ThueMorseMatvecCostR0Bridge.cpp", "Wh2ThueMorseMatvecCostR0Bridge.cpp",
               "Wh2ThueMorseMatvecCostR0.cpp", "Wh2RdpruTargetIdentityV2.cpp",
               "Wh2PublicBorrowedTargetIdentity.cpp", "Wh2FrozenTrace.cpp", "Wh2ThueMorseMatvecCostR0BridgeTest.cpp")
    expected_paths = [str(ROOT / name) for name in public_files] + [str(ROOT / "bench" / name) for name in support] + [str(generated)]
    exact(sorted(item["file"] for item in commands), sorted(expected_paths), "exact translation units")
    arms = []
    for item in commands:
        words = shlex.split(item["command"])
        require(Path(words[0]).resolve(strict=True) == Path("/usr/bin/g++").resolve(strict=True), "compiler")
        require(all(flag in words for flag in ("-std=c++11", "-O3", "-fno-lto", "-Werror")), "strict compiler flags")
        require("-fPIC" in words or "-fPIE" in words, "position independent compile")
        require(not any(flag in item["command"] for flag in ("WH_COUNT", "WIREHAIR_TESTING", "WIREHAIR_V2_ENABLE_TEST_HOOKS",
                    "-flto", "-fprofile", "-fsanitize", "-march=", "-mavx", "-mgfni", "WH_ALIGN64", "WH_HUGEPAGE")),
                "forbidden instrumentation/ISA")
        definitions = [word for word in words if word.startswith("-DWH2_MATVEC_COST_ARM")]
        if Path(item["file"]).name == "Wh2ThueMorseMatvecCostR0Bridge.cpp":
            require(len(definitions) == 1, "one bridge arm")
            arms.extend(definitions)
        else:
            require(not definitions, "bridge definition cannot leak")
    exact(sorted(arms), ["-DWH2_MATVEC_COST_ARM="+str(i) for i in range(4)], "four exact bridges")
    dependencies = {Path(item["file"]) for item in commands}
    for line in checked(["ninja", "-C", build, "-t", "deps"], deadline).decode().splitlines():
        if line.startswith("    "):
            path = Path(line.strip())
            dependencies.add(path if path.is_absolute() else build / path)
    require(any(path.name == "WirehairV2Solve.h" for path in dependencies) and
            any(path.name == "gf256.h" for path in dependencies), "native include dependencies")
    manifest = decode(read_regular(build / "native-inputs.json", 65536))
    exact(manifest["schema"], PROTOCOL + ".native-inputs.v1", "prebuild manifest schema")
    require(type(manifest["sources"]) is dict, "prebuild source mapping")
    native_names = {str(path.relative_to(ROOT)) for path in dependencies if ROOT in path.parents}
    require(native_names <= set(manifest["sources"]), "entire native dependency closure was pinned before build")
    for name, digest in manifest["sources"].items():
        require(not Path(name).is_absolute() and ".." not in Path(name).parts,
                "relative repository manifest path")
        exact(pin(ROOT / name)["sha256"], digest, "unchanged prebuild source: " + name)
    source_names = set(manifest["sources"]) | set(NEW_SOURCES) | {"AGENTS.md"}
    source_names.update("bench/" + name for name in (
        "Wh2AlignedIntermediateCostR0.py", "Wh2ThueMorseTinyPayloadCostR0.py",
        "Wh2ThueMorseNativeDataR0.py", "Wh2ThueMorseRecoveryHistoryR0.py", "Wh2NoncommutingRadixRunR0.py",
        "Wh2ThueMorseR0.py", "Wh2NoncommutingRadixR0.py"))
    head = checked(["git", "rev-parse", "HEAD"], deadline).decode().strip()
    sources = []
    for name in sorted(source_names):
        time_left(deadline)
        value = pin(ROOT / name)
        require(value["sha256"] == sha(checked(["git", "cat-file", "blob", head+":"+name], deadline)),
                "uncommitted source: " + name)
        sources.append(value)
    link_map = read_regular(build / "matvec-cost.map", 8 * 1024**2).decode()
    require("libwirehair.a(gf256.cpp.o)" not in link_map and
            "libmatvec_codec_gf.a(Wh2ThueMorseMatvecGfniR0.cpp.o)" in link_map, "single extracted GF object")
    symbols = checked(["/usr/bin/nm", "--defined-only", build / TARGET], deadline).decode().splitlines()
    for name in ("GF256Ctx", "gf256_init_", "gf256_mul_mem", "gf256_muladd_mem"):
        require(sum(line.split()[-1] == name for line in symbols if line.split()) == 1, "single GF symbol: " + name)
    outputs = {build / name for name in (TARGET, "wh2_matvec_cost_bridge_test", "lookup.bin", "matvec-cost.map",
        "matvec-cost-bridge-test.map", "compile_commands.json", "CMakeCache.txt", "build.ninja", "CMakeFiles/rules.ninja",
        "native-inputs.json", "codec-neutral/matvec-codec-source-manifest.json",
        "codec-neutral/wh2_matvec_codec_baseline_test", "codec-neutral/wh2_matvec_codec_candidate_test")}
    outputs.add(generated)
    outputs.update(build.rglob("*.o"))
    outputs.update(build.rglob("*.a"))
    installed = {Path(program).resolve(strict=True) for program in (
        "/usr/bin/cmake", "/usr/bin/ninja", "/usr/bin/g++", "/usr/bin/as", "/usr/bin/ld", "/usr/bin/git",
        "/usr/bin/ldd", "/usr/bin/nm", "/bin/bash", "/usr/bin/python3")}
    for name in ("cc1plus", "collect2"):
        installed.add(Path(checked(["g++", "-print-prog-name="+name], deadline).decode().strip()).resolve(strict=True))
    linked = A.runtime_libraries(checked(["ldd", build / TARGET], deadline).decode())
    installed.update(Path(entry["path"]) for entry in linked if entry["path"] is not None)
    dependency_pins = []
    for path in sorted(dependencies):
        time_left(deadline)
        dependency_pins.append(pin(path, installed=build not in path.parents and ROOT not in path.parents))
    return dict(protocol=PROTOCOL, schema=PROTOCOL+".claim.v1", build=str(build), source_head=head,
                sources=sources, dependencies=dependency_pins, outputs=[pin(path) for path in sorted(outputs)],
                installed=[pin(path, installed=True) for path in sorted(installed)],
                controller=pin(Path(__file__).resolve()), interpreter=pin(Path(sys.executable), installed=True),
                compiler_version=checked(["g++", "--version"], deadline).decode(), runtime_ldd=linked,
                target_cpu=50, raw_cap=RAW_CAP, stderr_cap=ERR_CAP, total_cap=TOTAL_CAP,
                worker_wall_seconds=WORKER_WALL, worker_cpu_seconds=WORKER_CPU,
                controller_wall_seconds=CONTROLLER_WALL, callbacks=CALLBACKS, batch=BATCH)


def controller_run(receipt_path):
    start = time.monotonic()
    deadline = start + CONTROLLER_WALL
    receipt_bytes = read_regular(receipt_path, 1024**2)
    receipt = decode(receipt_bytes)
    require(receipt_bytes == canonical(receipt), "canonical launch receipt")
    exact(receipt, current_receipt(receipt["build"], deadline), "prelaunch closure")
    os.mkdir(str(OUTPUT), 0o700)
    publish(OUTPUT / "CLAIM.json", receipt_bytes)
    raw, stderr, returncode, failure = b"", b"", None, None
    analysis = dict(outcome="INVALID", failure="worker not started")
    try:
        raw, stderr, returncode, failure = capture(Path(receipt["build"]) / TARGET, sha(receipt_bytes), deadline,
                                                 (OUTPUT / "raw.jsonl", OUTPUT / "stderr.txt"))
        exact(receipt, current_receipt(receipt["build"], deadline), "postlaunch closure")
        require(failure is None and returncode == 0 and not stderr, failure or "worker exit/stderr")
        analysis = reduce_raw(raw)
        time_left(deadline)
        require(analysis["claim_sha256"] == sha(receipt_bytes), "worker claim binding")
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        failure = str(error)
        analysis = dict(outcome="INVALID", failure=failure)
    for name in ("raw.jsonl", "stderr.txt"):
        if not (OUTPUT / name).exists():
            publish(OUTPUT / name, b"")
    if time.monotonic() >= deadline:
        failure = "controller deadline before sealing"
        analysis = dict(outcome="INVALID", failure=failure)
    publish(OUTPUT / "analysis.json", canonical(analysis))
    summary = dict(protocol=PROTOCOL, outcome=analysis["outcome"], failure=failure,
                   returncode=returncode, elapsed_seconds=time.monotonic()-start, raw_sha256=sha(raw),
                   production_promotion_claimed=False, all_K_claimed=False)
    publish(OUTPUT / "summary.json", canonical(summary))
    members = [pin(OUTPUT / name) for name in ("CLAIM.json", "raw.jsonl", "stderr.txt", "analysis.json", "summary.json")]
    require(sum(item["bytes"] for item in members) < TOTAL_CAP-65536, "bundle cap")
    publish(OUTPUT / "COMPLETE.json", canonical(dict(protocol=PROTOCOL, files=members, outcome=analysis["outcome"])))
    print(canonical(summary).decode(), end="")


def capture(executable, claim_hash, deadline, spool_paths):
    """Retain bounded raw prefixes and always reap this exact child process."""
    buffers = [bytearray(), bytearray()]
    failure = None
    # Durable prefixes survive even an emergency kill of the controller. All
    # descendants share the outer observer's dedicated process group.
    spools = []
    child = None
    selector = selectors.DefaultSelector()
    try:
        for path in spool_paths:
            spools.append(os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600))
        child = subprocess.Popen([str(executable), "--worker", claim_hash],
                                 stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, close_fds=True)
        for index, stream in enumerate((child.stdout, child.stderr)):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, index)
        while selector.get_map():
            remaining = time_left(deadline)
            for key, _ in selector.select(min(remaining, .05)):
                block = os.read(key.fileobj.fileno(), 65536)
                if not block:
                    selector.unregister(key.fileobj)
                    continue
                index = key.data
                cap = (RAW_CAP, ERR_CAP)[index]
                allowed = cap - len(buffers[index])
                retained = block[:allowed]
                pending = memoryview(retained)
                while pending:
                    written = os.write(spools[index], pending)
                    require(written > 0, "spool write made no progress")
                    buffers[index].extend(pending[:written])
                    pending = pending[written:]
                require(len(block) <= allowed, "worker output cap")
        child.wait(timeout=time_left(deadline))
    except (ValueError, OSError, subprocess.TimeoutExpired) as error:
        failure = str(error)
    finally:
        if child is not None:
            if child.poll() is None:
                try:
                    child.kill()
                except ProcessLookupError:
                    pass
            child.wait()
            for stream in (child.stdout, child.stderr):
                try:
                    stream.close()
                except OSError as error:
                    failure = failure or "pipe cleanup: " + str(error)
        for fd in spools:
            for action in (lambda: os.fsync(fd), lambda: os.fchmod(fd, 0o400), lambda: os.close(fd)):
                try:
                    action()
                except OSError as error:
                    failure = failure or "spool cleanup: " + str(error)
        try:
            selector.close()
        except OSError as error:
            failure = failure or "selector cleanup: " + str(error)
    return bytes(buffers[0]), bytes(buffers[1]), None if child is None else child.returncode, failure


def run(receipt_path):
    """Owned process-group observer for this bounded full-cost screen."""
    buffers = [bytearray(), bytearray()]
    deadline = time.monotonic() + OBSERVER_WALL
    selector = selectors.DefaultSelector()
    child = None
    try:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "_controller", str(receipt_path)],
                                 stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 close_fds=True, start_new_session=True)
        for index, stream in enumerate((child.stdout, child.stderr)):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, index)
        while selector.get_map():
            for key, _ in selector.select(min(time_left(deadline), .05)):
                block = os.read(key.fileobj.fileno(), 65536)
                if not block:
                    selector.unregister(key.fileobj)
                else:
                    buffers[key.data].extend(block)
                    require(len(buffers[key.data]) <= ERR_CAP, "outer observer output cap")
        while os.waitid(os.P_PID, child.pid, os.WEXITED | os.WNOWAIT | os.WNOHANG) is None:
            time_left(deadline)
            time.sleep(.001)
    finally:
        try:
            if child is not None:
                os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        finally:
            if child is not None:
                child.wait()
                for stream in (child.stdout, child.stderr):
                    stream.close()
            selector.close()
    sys.stdout.buffer.write(buffers[0])
    sys.stderr.buffer.write(buffers[1])
    require(child.returncode == 0, "controller failed; preserve namespace and prefixes")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    lookup = sub.add_parser("lookup")
    lookup.add_argument("output", type=Path)
    receipt = sub.add_parser("receipt")
    receipt.add_argument("build", type=Path)
    receipt.add_argument("output", type=Path)
    launch = sub.add_parser("run")
    launch.add_argument("receipt", type=Path)
    controller = sub.add_parser("_controller", help=argparse.SUPPRESS)
    controller.add_argument("receipt", type=Path)
    replay = sub.add_parser("replay")
    replay.add_argument("raw", type=Path)
    args = parser.parse_args()
    if args.command == "lookup":
        write_lookup(args.output)
    elif args.command == "receipt":
        publish(args.output, canonical(current_receipt(args.build, time.monotonic() + 60)))
    elif args.command == "run":
        run(args.receipt)
    elif args.command == "_controller":
        controller_run(args.receipt)
    else:
        print(canonical(reduce_raw(read_regular(args.raw, RAW_CAP))).decode(), end="")


if __name__ == "__main__":
    main()
