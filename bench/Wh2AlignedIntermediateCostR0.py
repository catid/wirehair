#!/usr/bin/env python3
"""One-shot full-cost alignment screen (.68); replay never executes a codec.

The scientific contract is frozen in Beads. This additive controller does not
reuse or re-score any historical experiment. Build and neutral tests are
separate commands; only `run` can consume the new scientific namespace.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import selectors
import shlex
import signal
import stat
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = "wirehair.wh2.aligned-intermediate-cost-r0"
OUTPUT = Path("/var/tmp/wh2-aligned-intermediate-cost-r0")
TARGET = "wh2_aligned_intermediate_cost_r0"
RAW_CAP, ERR_CAP, TOTAL_CAP = 4 * 1024**2, 65536, 12 * 1024**2
T11 = 2.200985160082949
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1)
PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2))
METRICS = ("create", "repair", "encoder_lifecycle", "decoder_lifecycle")
NEW_SOURCES = (
    "bench/Wh2AlignedIntermediateCostR0.py",
    "bench/test_Wh2AlignedIntermediateCostR0.py",
    "bench/Wh2AlignedIntermediateCostR0.cpp",
    "bench/Wh2AlignedIntermediateCostR0Bridge.h",
    "bench/Wh2AlignedIntermediateCostR0Inspect.inc",
    "bench/Wh2AlignedIntermediateCostR0/CMakeLists.txt",
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def integer(value, low=0, high=(1 << 63) - 1):
    require(type(value) is int and low <= value <= high, "integer bounds/type")
    return value


def exact(actual, expected, message):
    """JSON booleans/floats must never impersonate integral evidence."""
    require(type(actual) is type(expected), message + " type")
    if isinstance(expected, dict):
        require(set(actual) == set(expected), message + " keys")
        for key in expected:
            exact(actual[key], expected[key], message + "." + key)
    elif isinstance(expected, (list, tuple)):
        require(len(actual) == len(expected), message + " length")
        for a, b in zip(actual, expected):
            exact(a, b, message)
    else:
        require(actual == expected, message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def decode(data):
    return json.loads(data, object_pairs_hook=unique_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(
                          ValueError("nonfinite JSON constant: " + value)))


def time_left(deadline):
    left = deadline - time.monotonic()
    require(left > 0, "controller deadline")
    return left


def read_regular(path, cap, installed=False):
    path = Path(path)
    if installed:
        path = path.resolve(strict=True)
    fd = os.open(str(path), os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_size <= cap and
                (installed or before.st_nlink == 1), "regular bounded file: " + str(path))
        chunks, size = [], 0
        while True:
            part = os.read(fd, min(65536, cap + 1 - size))
            if not part:
                break
            chunks.append(part)
            size += len(part)
            require(size <= cap, "file grew past cap")
        after = os.fstat(fd)
        require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
                 before.st_ctime_ns, before.st_nlink) ==
                (after.st_dev, after.st_ino, size, after.st_mtime_ns,
                 after.st_ctime_ns, after.st_nlink), "file changed during read")
        return b"".join(chunks)
    finally:
        os.close(fd)


def pin(path, installed=False):
    path = Path(path).resolve(strict=True) if installed else Path(path)
    data = read_regular(path, 256 * 1024**2, installed=installed)
    return {"path": str(path), "bytes": len(data), "sha256": sha(data)}


def checked(args, deadline, cwd=ROOT):
    args = list(args)
    if args[0] in ("git", "ninja", "g++", "ldd"):
        args[0] = "/usr/bin/" + args[0]
    result = subprocess.run([str(x) for x in args], cwd=str(cwd),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=time_left(deadline), check=False)
    require(result.returncode == 0, "command failed: " + str(args) + ": " +
            result.stderr.decode("utf-8", errors="replace")[:2000])
    require(len(result.stdout) <= 16 * 1024**2 and len(result.stderr) <= ERR_CAP,
            "inspection command output cap")
    return result.stdout


def roster():
    index = 0
    for rep in range(12):
        for step in range(2):
            order = (rep + step) % 2
            for metric_step in range(4):
                metric = (rep + step + metric_step) % 4
                for class_step in range(5):
                    comparison = (2 * rep + step + class_step) % 5
                    for position, side in enumerate(SIDES):
                        pair_index = rep % 4 if position < 2 else (position - 2) // 2
                        q = ((2 * (rep + 12 * pair_index) + 1) * 1000000) // 96
                        yield dict(index=index, rep=rep, order=order, metric=metric,
                                   **{"class": comparison}, position=position,
                                   arm=PAIRS[comparison][side ^ order], q=q)
                        index += 1


def confidence(values):
    require(len(values) == 12 and all(math.isfinite(x) for x in values),
            "t11 requires twelve finite replicate contrasts")
    mean = math.fsum(values) / 12
    radius = T11 * math.sqrt(math.fsum((x - mean)**2 for x in values) / 11 / 12)
    return dict(replicates=12, mean_log=mean, lower95_log=mean - radius,
                upper95_log=mean + radius, ratio=math.exp(mean),
                lower95=math.exp(mean - radius), upper95=math.exp(mean + radius))


def statistics(records):
    """No pooling of order/physical lanes and no deletion of observations."""
    groups = {}
    for start in range(0, len(records), 10):
        panel = records[start:start + 10]
        require(len(panel) == 10, "incomplete panel")
        first = panel[0]
        contrasts = []
        for j in range(4):
            pair = panel[2 + 2*j:4 + 2*j]
            by_side = {}
            for row in pair:
                side = SIDES[row["position"]] ^ row["order"]
                by_side[side] = row["observation"][3] - row["observation"][2]
            require(set(by_side) == {0, 1} and min(by_side.values()) > 0,
                    "complete positive-time adjacent pair")
            contrasts.append(math.log(by_side[1]) - math.log(by_side[0]))
        key = (first["metric"], first["class"], first["order"])
        groups.setdefault(key, []).append(math.fsum(contrasts) / 4)
    require(len(groups) == 40, "all conditional statistical groups")
    results, failed_controls, failed_treatments = [], [], []
    bound = math.log1p(.02)
    for (metric, comparison, order), values in sorted(groups.items()):
        estimate = confidence(values)
        row = dict(metric=METRICS[metric], comparison=("WW", "PP", "AA", "PA", "WA")[comparison],
                   order=order, estimate=estimate, replicate_logs=values)
        if comparison < 3:
            row["control_pass"] = estimate["lower95_log"] > -bound and estimate["upper95_log"] < bound
            if not row["control_pass"]:
                failed_controls.append([metric, comparison, order])
        else:
            limit = bound if comparison == 3 and metric in (0, 3) else 0.0
            row["treatment_pass"] = estimate["upper95_log"] < limit
            row["upper_ratio_limit"] = math.exp(limit)
            if not row["treatment_pass"]:
                failed_treatments.append([metric, comparison, order])
        results.append(row)
    outcome = "CONTROL_FAIL" if failed_controls else "FAIL" if failed_treatments else "PASS"
    return dict(outcome=outcome, statistics=results, failed_controls=failed_controls,
                failed_treatments=failed_treatments, all_K_claimed=False,
                production_promotion_claimed=False)


def hex_bytes(value, size):
    require(type(value) is str and len(value) == 2 * size, "hex byte count")
    result = bytes.fromhex(value)
    require(result.hex() == value, "canonical hex")
    return result


def clocks(observation, previous=None):
    require(type(observation) is list and len(observation) == 8, "clock envelope shape")
    m0, c0, m1, m2, c1, m3 = [integer(x) for x in observation[:6]]
    require(m0 <= m1 < m2 <= m3 and c0 <= c1, "clock ordering")
    for usage in observation[6:]:
        require(type(usage) is list and len(usage) == 4, "usage shape")
        for value in usage:
            integer(value)
    require(all(a <= b for a, b in zip(observation[6], observation[7])), "usage ordering")
    if previous is not None:
        require(previous[5] <= m0 and previous[4] <= c0, "cross-record clocks")
        require(all(a <= b for a, b in zip(previous[7], observation[6])), "cross-record usage")
    return m2 - m1


def verify_identity(identity):
    expected = dict(family=26, model=8, stepping=1, full_apic_id=100,
                    initial_apic_id_8=100, core_id=50, package_id=0, thread_id=0,
                    threads_per_core=2, ccd_id=6, complex_id=6, logical_processors_per_package=128)
    exact(identity["derived"], expected, "frozen physical target")
    for key in ("before_cpu", "after_cpu", "requested_cpu"):
        exact(identity[key], 50, "CPU50 identity")
    for key in ("raw_capture_complete", "semantic_validation_passed", "singleton_affinity_verified"):
        require(identity[key] is True, "target identity gate")
    data = hex_bytes(identity["canonical_hex"], integer(identity["canonical_bytes"], 1, 4096))
    require(sha(data) == identity["canonical_sha256"], "canonical identity hash")
    require(type(identity["capabilities"]) is dict, "CPU capability fields")
    require(set(identity["capabilities"]) == {"leaf1_ecx", "leaf1_edx", "leaf6_eax", "leaf6_ecx",
            "leaf80000001_ecx", "leaf80000001_edx", "leaf80000008_ebx", "leaf80000021_eax",
            "max_basic_leaf", "max_extended_leaf"}, "CPU capability schema")
    for value in identity["capabilities"].values():
        integer(value, 0, (1 << 32) - 1)


def verify_runtime(runtime):
    for key, value in dict(polynomial=0x14d, ssse3=1, avx2=1, gfni=1, avx512=1).items():
        exact(runtime[key], value, "frozen GF256 runtime")
    integer(runtime["address"], 1)


def verify_snapshot(snapshot, arm, rep, lane, source, profile, master):
    exact((snapshot["arm"], snapshot["rep"], snapshot["lane"]), (arm, rep, lane), "snapshot roster")
    integer(snapshot["handle"], 1)
    integer(snapshot["source"], 1, (1 << 63) - 7680)
    require(snapshot["source_sha256"] == sha(source), "snapshot immutable source")
    if arm == 0:
        require(all(snapshot[key] is None for key in ("intermediate", "intermediate_bytes",
                "intermediate_capacity", "source_policy", "profile_hex", "intermediate_sha256")), "WH1 opaque state")
    else:
        integer(snapshot["intermediate"], 1, (1 << 63) - 46080)
        exact(snapshot["intermediate_bytes"], 46080, "retained vector bytes")
        exact(snapshot["intermediate_capacity"], 46080, "retained vector capacity")
        exact(snapshot["source_policy"], 2, "borrowed ownership enum")
        require(snapshot["profile_hex"] == profile.hex(),
                "snapshot profile/borrowed ownership")
        require(snapshot["intermediate_sha256"] == sha(master), "snapshot intermediate bytes")
        require(arm != 2 or snapshot["intermediate"] % 32 == 0, "candidate native32 alignment")


def verify_preflight(header):
    source = hex_bytes(header["source_hex"], 7680)
    require(source == bytes((37*i + i//11) & 255 for i in range(7680)), "frozen source law")
    require(header["source_sha256"] == sha(source), "source hash")
    pre = header["preflight"]
    require(len(pre["packets_hex"]) == len(pre["packet_sha256"]) == 3, "three packet corpora")
    packets = [hex_bytes(value, 15360) for value in pre["packets_hex"]]
    require([sha(value) for value in packets] == pre["packet_sha256"], "packet hashes")
    require(packets[1] == packets[2], "public/candidate exact packets")
    require(all(value[:7680] == source for value in packets), "systematic packet oracle")
    require(len(pre["intermediate_hex"]) == len(pre["intermediate_sha256"]) == 2,
            "two intermediate corpora")
    masters = [hex_bytes(value, 46080) for value in pre["intermediate_hex"]]
    require(masters[0] == masters[1] and [sha(x) for x in masters] == pre["intermediate_sha256"],
            "exact intermediate bytes/hashes")
    require(len(pre["profiles_hex"]) == 3 and pre["profiles_hex"][0] is None, "profile corpora")
    profiles = [hex_bytes(x, 32) for x in pre["profiles_hex"][1:]]
    require(profiles[0] == profiles[1] and profiles[0][:28] == b"WHV2" + struct.pack(
        "<HHQQI", 1, 32, 0x4b295bbb47f4f9c9, 7680, 1280) and profiles[0][29:] == bytes(3),
        "ordinary certified profile binding")
    require(type(pre["columns"]) is list and len(pre["columns"]) == 2, "two original column corpora")
    exact(pre["columns"][0], pre["columns"][1], "same original row columns")
    columns = pre["columns"][0]
    require(len(columns) == 12, "all original systematic/repair rows")
    for block_id in range(12):
        row = columns[block_id]
        require(type(row) is list and 4 <= len(row) <= 9 and len(set(row)) == len(row), "repair row terms")
        require(sum(integer(column, 0, 35) >= 6 for column in row) == 3, "three mix columns")
        expected = bytearray(1280)
        for column in row:
            for offset, value in enumerate(masters[0][column * 1280:(column + 1) * 1280]):
                expected[offset] ^= value
        require(bytes(expected) == packets[1][block_id * 1280:(block_id + 1) * 1280], "independent scalar XOR oracle")
    require(len(pre["decode_statuses"]) == len(pre["first_success"]) == 3, "three decoder references")
    for arm in range(3):
        first = integer(pre["first_success"][arm], 5, 11)
        exact(pre["decode_statuses"][arm], [1] * first + [0], "decode until first success")
    require(pre["decode_statuses"][1] == pre["decode_statuses"][2], "P/A recovery behavior")
    snapshots = pre["snapshots"]
    require(type(snapshots) is list and len(snapshots) == 72, "all natural repair handles")
    by_key, handles, spans = {}, set(), []
    cursor = 0
    for rep in range(12):
        for i in range(3):
            arm = (rep + i) % 3
            for lane in range(2):
                snapshot = snapshots[cursor]
                cursor += 1
                verify_snapshot(snapshot, arm, rep, lane, source, profiles[0], masters[0])
                require(snapshot["handle"] not in handles, "distinct live handles")
                handles.add(snapshot["handle"])
                by_key[(rep, lane, arm)] = snapshot
                if arm:
                    spans.append((snapshot["intermediate"], snapshot["intermediate"] + 46080))
    require(len({s["source"] for s in snapshots}) == 1, "one shared borrowed source")
    source_address = snapshots[0]["source"]
    require(source_address >= 64, "source guard address")
    spans.append((source_address - 64, source_address + 7680 + 64))
    require(len(header["outputs"]) == 2 and header["outputs"][0] != header["outputs"][1], "two common output lanes")
    for address in header["outputs"]:
        integer(address, 1, (1 << 63) - 16896)
        require(address % 4096 == 0, "output allocation alignment")
        spans.append((address, address + 16896))
    for index, left in enumerate(spans):
        for right in spans[index + 1:]:
            require(left[1] <= right[0] or right[1] <= left[0], "disjoint live storage")
    for handle in handles:
        require(all(not (low <= handle < high) for low, high in spans), "handle address overlaps payload allocation")
    return source, profiles[0], masters[0], by_key


def reduce_raw(raw):
    require(type(raw) is bytes and 0 < len(raw) <= RAW_CAP and raw.endswith(b"\n"), "complete bounded JSONL")
    lines = raw.splitlines()
    require(len(lines) == 4802, "header + all4800 callbacks + footer")
    header, *records, footer = [decode(line) for line in lines]
    require(header["type"] == "header" and footer["type"] == "footer" and
            header["protocol"] == footer["protocol"] == PROTOCOL and
            header["schema"] == PROTOCOL + ".raw.v1", "protocol/schema framing")
    require(footer["outcome"] == "COMPLETE" and footer["failure"] is None, "worker complete correctness")
    exact(header["target_cpu"], 50, "frozen CPU")
    hex_bytes(header["claim_sha256"], 32)
    verify_identity(header["identity_before"])
    verify_identity(footer["identity_after"])
    exact(header["identity_before"], footer["identity_after"], "normalized CPU identity stability")
    require(header["identity_before"]["capabilities"] == footer["identity_after"]["capabilities"],
            "CPU capabilities changed")
    verify_runtime(header["runtime_before"])
    verify_runtime(footer["runtime_after"])
    require(header["runtime_before"] == footer["runtime_after"], "shared GF256 runtime changed")
    exact(header["clock_resolutions"], [1, 1], "clock resolutions")
    source, profile, master, snapshots = verify_preflight(header)
    prelude = header["prelude"]
    exact((prelude["iterations"], prelude["seed"], prelude["final_state"]),
          (1 << 20, 0x9e3779b97f4a7c15, 0x43935dad1647741b), "fixed single prelude")
    clocks(prelude["observation"])
    start, start_cpu = integer(header["worker_start_ns"]), integer(header["worker_start_cpu_ns"])
    end, end_cpu = integer(footer["worker_end_ns"]), integer(footer["worker_end_cpu_ns"])
    require(start <= prelude["observation"][0] and start_cpu <= prelude["observation"][1], "prelude lifetime")
    require(0 < end - start < 15 * 10**9 and 0 <= end_cpu - start_cpu <= 14 * 10**9, "worker caps")
    now, t0 = integer(footer["schedule_now_ns"]), integer(footer["t0"])
    require(now >= prelude["observation"][5] and t0 == (now // 2000000 + 2) * 2000000, "fixed schedule epoch")
    previous = prelude["observation"]
    work_sum = wait_wall = wait_cpu = call_sum = 0
    for expected, row in zip(roster(), records):
        require(row["type"] == "record" and all(type(row[key]) is int and row[key] == value
                for key, value in expected.items()), "exact callback roster")
        require(row["checked"] is True and row["work"]["complete"] is True, "all callback correctness checks")
        target = integer(row["target"])
        require(target == t0 + 2000000 * row["index"] + row["q"], "assigned full-phase target")
        ready = integer(row["ready"])
        require(previous[5] <= ready <= target - 100000, "preparation/cleanup deadline")
        wait = row["wait"]
        require(type(wait) is list and len(wait) == 4, "wait accounting shape")
        wm0, wc0, wm1, wc1 = [integer(x) for x in wait]
        require(ready <= wm0 <= target <= wm1 <= row["observation"][0] and
                previous[4] <= wc0 <= wc1 <= row["observation"][1], "wait ordering")
        work_sum += clocks(row["observation"], previous)
        require(target <= row["observation"][2] <= target + 5000, "late-start invalidation")
        wait_wall += wm1 - wm0
        wait_cpu += wc1 - wc0
        previous = row["observation"]
        work = row["work"]
        metric, arm, rep, lane = row["metric"], row["arm"], row["rep"], row["order"]
        integer(work["handle"], 1)
        require(type(work["freed"]) is bool and work["freed"] == (metric in (2, 3)), "timed free boundary")
        exact(work["create_code"], None if metric == 1 else 0, "fresh constructor success")
        integer(work["calls"], 1, 384)
        if metric == 3:
            first = header["preflight"]["first_success"][arm]
            exact(work["codes"], header["preflight"]["decode_statuses"][arm], "decoder status sequence")
            exact(work["lengths"], [1280] * (first + 1), "decoder lengths")
            exact(work["first_success"], first, "first decoder success")
            exact(work["recover_code"], 0, "recovery success")
            exact(work["bytes"], 7680, "recovered message length")
            exact(work["calls"], first + 4, "complete decoder lifecycle")
        else:
            count = (0, 6, 12)[metric]
            exact(work["codes"], [0] * count, "every encode status aggregate")
            exact(work["lengths"], [1280] * count, "every encode length aggregate")
            require(work["first_success"] is None and work["recover_code"] is None and work["bytes"] is None,
                    "nondecoder result fields")
            require(work["calls"] == (1, 384, 14)[metric], "exact full metric call count")
        if metric == 0:
            require(row["snapshot"] is not None, "post-create borrowed snapshot")
            verify_snapshot(row["snapshot"], arm, rep, lane, source, profile, master)
            require(row["snapshot"]["handle"] == work["handle"] and
                    row["snapshot"]["source"] == snapshots[(rep, lane, arm)]["source"], "fresh handle/source binding")
        else:
            require(row["snapshot"] is None, "no fabricated freed-handle snapshot")
        if metric == 1:
            require(work["handle"] == snapshots[(rep, lane, arm)]["handle"], "exact same physical AA handle")
        call_sum += work["calls"]
    require(previous[5] <= end and previous[4] <= end_cpu and work_sum <= 10**9, "complete measured lifetime/work cap")
    for key, expected in dict(records=4800, callbacks=4800, checked_callbacks=4800,
            work_calls=call_sum, sum_work_ns=work_sum, sum_wait_wall_ns=wait_wall,
            sum_wait_cpu_ns=wait_cpu).items():
        require(type(footer[key]) is int and footer[key] == expected, "footer accounting: " + key)
    require(wait_wall <= end - start and wait_cpu <= end_cpu - start_cpu, "wait work is included in process cost")
    result = statistics(records)
    result.update(protocol=PROTOCOL, claim_sha256=header["claim_sha256"], raw_sha256=sha(raw),
                  callbacks=4800, work_calls=call_sum, sum_work_ns=work_sum,
                  conditioning=dict(prelude_wall_ns=prelude["observation"][3] - prelude["observation"][2],
                                    wait_wall_ns=wait_wall, wait_cpu_ns=wait_cpu,
                                    worker_wall_ns=end - start, worker_cpu_ns=end_cpu - start_cpu),
                  memory=dict(intermediate_bytes=46080, baseline_capacity=46080,
                              candidate_requested_allocation_bytes=46112, extra_logical_copy_bytes=46080,
                              publication_two_buffer_requested_bytes=92192,
                              total_codec_or_process_peak_claimed=False))
    return result


def current_receipt(build, deadline):
    """Pin actual fresh build closure, not a historical archive identity."""
    build = Path(build).resolve(strict=True)
    require(build.is_dir() and ROOT not in build.parents and build != ROOT,
            "fresh external build directory")
    require(checked(["ninja", "-C", build, "-n"], deadline).decode().strip().endswith(
        "ninja: no work to do."), "build is dirty")
    cache = read_regular(build / "CMakeCache.txt", 1024**2).decode()
    require("WH2_ALIGNMENT_COST_SANITIZERS:BOOL=OFF" in cache,
            "scientific build must have sanitizers OFF")
    commands = decode(read_regular(build / "compile_commands.json", 2 * 1024**2))
    require(type(commands) is list and len(commands) == 29, "exact29 full variant translation units")
    for item in commands:
        command = item["command"]
        words = shlex.split(command)
        require(Path(words[0]).resolve(strict=True) == Path("/usr/bin/g++").resolve(strict=True), "actual compiler identity")
        require(all(flag in words for flag in
                    ("-std=c++11", "-O3", "-march=native", "-fno-lto")), "compiler mode")
        require(not any(flag in command for flag in
                        ("WIREHAIR_TESTING", "WIREHAIR_V2_ENABLE_TEST_HOOKS", "-fsanitize", "-flto",
                         "WH_ALIGN64", "WH_HUGEPAGE", "-fprofile")), "forbidden compiler instrumentation")
    dependencies = {Path(item["file"]) for item in commands}
    dependency_text = checked(["ninja", "-C", build, "-t", "deps"], deadline).decode()
    for line in dependency_text.splitlines():
        if line.startswith("    "):
            path = Path(line.strip())
            dependencies.add(path if path.is_absolute() else build / path)
    require(any(path.name == "vector.tcc" for path in dependencies), "standard-library dependencies missing")
    require(any(path.name == "WirehairV2Solve.cpp" for path in dependencies), "solver dependency missing")
    require(not any(path.name == "wirehair.hpp" for path in dependencies), "C++ wrapper ODR hazard")
    source_names = set(NEW_SOURCES) | {"AGENTS.md", "bench/Wh2AlignedIntermediateR0/CMakeLists.txt"}
    for path in dependencies:
        if ROOT in path.parents:
            source_names.add(str(path.relative_to(ROOT)))
    for manifest in build.glob("*/alignment-source-manifest.txt"):
        for line in read_regular(manifest, 1024**2).decode().splitlines():
            if line.startswith("input "):
                _, expected_sha, name = line.split(" ", 2)
                require(pin(ROOT / name)["sha256"] == expected_sha, "original mirror input changed")
                source_names.add(name)
    cost_manifest = read_regular(build / "alignment-cost-manifest.txt", 1024**2).decode()
    for line in cost_manifest.splitlines():
        if line.startswith(("input ", "file ", "preappend_manifest ")):
            kind, expected_sha, path_string = line.split(" ", 2)
            path = Path(path_string)
            require(path.is_absolute(), "absolute cost manifest path")
            require(pin(path)["sha256"] == expected_sha, "final mirror/support manifest changed: " + path_string)
            if kind == "input":
                require(ROOT in path.parents, "repository cost input")
                source_names.add(str(path.relative_to(ROOT)))
            else:
                require(build in path.parents, "owned build mirror output")
                dependencies.add(path)
    require("codec/WirehairV2Profile.cpp" in source_names and "gf256.cpp" in source_names,
            "original production inputs missing from manifests")
    head = checked(["git", "rev-parse", "HEAD"], deadline).decode().strip()
    sources = []
    for name in sorted(source_names):
        time_left(deadline)
        value = pin(ROOT / name)
        require(value["sha256"] == sha(checked(["git", "cat-file", "blob", head + ":" + name], deadline)),
                "source is not committed: " + name)
        sources.append(value)
    outputs = {build / name for name in
               (TARGET, "alignment-cost.map", "alignment-cost-manifest.txt", "compile_commands.json",
                "CMakeCache.txt", "build.ninja", "CMakeFiles/rules.ninja")}
    outputs.update(build.rglob("*.o"))
    outputs.update(build.rglob("*.a"))
    outputs.update(build.glob("*/alignment-source-manifest.txt"))
    outputs.update(build.glob("*/alignment-source-owner.txt"))
    installed = set()
    for program in ("/usr/bin/cmake", "/usr/bin/ninja", "/usr/bin/g++", "/usr/bin/as", "/usr/bin/ld",
                    "/usr/bin/git", "/usr/bin/ldd", "/bin/bash"):
        installed.add(Path(program).resolve(strict=True))
    for name in ("cc1plus", "collect2"):
        installed.add(Path(checked(["g++", "-print-prog-name=" + name], deadline).decode().strip()).resolve(strict=True))
    linked = checked(["ldd", build / TARGET], deadline).decode()
    for word in linked.split():
        if word.startswith("/"):
            installed.add(Path(word).resolve(strict=True))
    require(any(path.name.startswith("libc.so") or path.name.startswith("libc-") for path in installed),
            "libc runtime closure missing")
    dependency_pins = []
    for path in sorted(dependencies):
        time_left(deadline)
        dependency_pins.append(pin(path, installed=build not in path.parents and ROOT not in path.parents))
    return dict(protocol=PROTOCOL, schema=PROTOCOL + ".claim.v1", build=str(build), source_head=head,
                sources=sources, dependencies=dependency_pins,
                outputs=[pin(path) for path in sorted(outputs)],
                installed=[pin(path, installed=True) for path in sorted(installed)],
                controller=pin(Path(__file__).resolve()), interpreter=pin(Path(os.sys.executable), installed=True),
                compiler_version=checked(["g++", "--version"], deadline).decode(),
                runtime_ldd=linked, target_cpu=50, raw_cap=RAW_CAP, stderr_cap=ERR_CAP,
                worker_wall_seconds=15, worker_cpu_seconds=14, controller_wall_seconds=20)


def publish(path, data):
    fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    try:
        with os.fdopen(fd, "wb", closefd=False) as output:
            output.write(data)
            output.flush()
            os.fsync(fd)
        os.fchmod(fd, 0o400)
    finally:
        os.close(fd)


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


def controller_run(receipt_path):
    start = time.monotonic()
    deadline = start + 20
    receipt_bytes = read_regular(receipt_path, 1024**2)
    receipt = decode(receipt_bytes)
    require(receipt_bytes == canonical(receipt), "canonical launch receipt")
    require(receipt == current_receipt(receipt["build"], deadline), "prelaunch closure changed")
    # Existing/partial namespace is permanently spent. Never delete or retry it.
    os.mkdir(str(OUTPUT), 0o700)
    publish(OUTPUT / "CLAIM.json", receipt_bytes)
    raw, stderr, returncode, failure = b"", b"", None, None
    analysis = dict(outcome="INVALID", failure="worker not started")
    try:
        raw, stderr, returncode, failure = capture(Path(receipt["build"]) / TARGET, sha(receipt_bytes), deadline,
                                                (OUTPUT / "raw.jsonl", OUTPUT / "stderr.txt"))
        require(receipt == current_receipt(receipt["build"], deadline), "postlaunch closure changed")
        require(failure is None and returncode == 0 and not stderr, failure or "worker exit/stderr")
        analysis = reduce_raw(raw)
        time_left(deadline)
        require(analysis["claim_sha256"] == sha(receipt_bytes), "worker claim binding")
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        failure = str(error)
        analysis = dict(outcome="INVALID", failure=failure)
    # Capture already sealed the durable raw/stderr prefixes.
    for name in ("raw.jsonl", "stderr.txt"):
        if not (OUTPUT / name).exists():
            publish(OUTPUT / name, b"")
    if time.monotonic() >= deadline:
        failure = "controller deadline before sealing"
        analysis = dict(outcome="INVALID", failure=failure)
    publish(OUTPUT / "analysis.json", canonical(analysis))
    summary = dict(protocol=PROTOCOL, outcome=analysis["outcome"], failure=failure,
                   returncode=returncode, elapsed_seconds=time.monotonic() - start,
                   raw_sha256=sha(raw), production_promotion_claimed=False, all_K_claimed=False)
    publish(OUTPUT / "summary.json", canonical(summary))
    members = [pin(OUTPUT / name) for name in ("CLAIM.json", "raw.jsonl", "stderr.txt", "analysis.json", "summary.json")]
    require(sum(item["bytes"] for item in members) < TOTAL_CAP - 65536, "bundle cap")
    publish(OUTPUT / "COMPLETE.json", canonical(dict(protocol=PROTOCOL, files=members, outcome=analysis["outcome"])))
    print(canonical(summary).decode(), end="")


def run(receipt_path):
    """25s emergency observer; the child controller has its own 20s budget."""
    buffers = [bytearray(), bytearray()]
    deadline = time.monotonic() + 25
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
        # Observe exit WITHOUT reaping the leader: its reserved PID protects
        # against group-id reuse until we have cleaned the whole owned group.
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
    require(child.returncode == 0, "controller failed; preserve existing namespace and prefixes")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
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
    if args.command == "receipt":
        publish(args.output, canonical(current_receipt(args.build, time.monotonic() + 60)))
    elif args.command == "run":
        run(args.receipt)
    elif args.command == "_controller":
        controller_run(args.receipt)
    else:
        print(canonical(reduce_raw(read_regular(args.raw, RAW_CAP))).decode(), end="")


if __name__ == "__main__":
    main()
