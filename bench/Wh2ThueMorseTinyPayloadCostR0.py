#!/usr/bin/env python3
"""One-shot .69 encoder-cost reject screen; build/replay never execute a codec."""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import selectors
import shlex
import signal
import subprocess
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Reuse the reviewed bounded reads, JSON/types, clocks, capture and publisher.
# Do not mutate that module's globals or its spent experiment namespace.
A = sibling("_tiny_cost_common", "Wh2AlignedIntermediateCostR0.py")
require, exact, integer = A.require, A.exact, A.integer
canonical, decode, sha, hex_bytes = A.canonical, A.decode, A.sha, A.hex_bytes
read_regular, pin, checked, publish = A.read_regular, A.pin, A.checked, A.publish
time_left = A.time_left
ROOT = A.ROOT
PROTOCOL = "wirehair.wh2.thue-tiny-payload-cost-r0"
OUTPUT = Path("/var/tmp/wh2-thue-tiny-payload-cost-r0")
TARGET = "wh2_thue_tiny_payload_cost_r0"
RAW_CAP, ERR_CAP, TOTAL_CAP = A.RAW_CAP, A.ERR_CAP, A.TOTAL_CAP
LOOKUP_SHA = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d"
WIDTHS = (2, 64, 1280)
PAIRS = ((0, 0), (1, 1), (0, 1))
SIDES = A.SIDES
NEW_SOURCES = tuple("bench/" + name for name in (
    "Wh2ThueMorseTinyPayloadCostR0.py", "test_Wh2ThueMorseTinyPayloadCostR0.py",
    "Wh2ThueMorseTinyPayloadCostR0.cpp", "Wh2ThueMorseTinyPayloadCostR0Bridge.h",
    "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp", "Wh2ThueMorseTinyPayloadCostR0/CMakeLists.txt"))


def lookup_bytes():
    """Reconstruct the already selected companion tables, never score rows/rank."""
    D = sibling("_tiny_cost_tables", "Wh2ThueMorseNativeDataR0.py")
    pair = []
    for first in (124, 125):
        matrix = [[0] * 6 for _ in range(6)]
        for i in range(5):
            matrix[i + 1][i] = 1
        for i, value in enumerate((first, 127, 152, 84, 241, 63)):
            matrix[i][5] = value
        pair.append(matrix)
    packed, _ = D.build_lookup(pair)
    require(len(packed) == 39936 and sha(packed) == LOOKUP_SHA, "fixed selected lookup")
    return packed


def write_lookup(path):
    packed = lookup_bytes()
    path = Path(path)
    if path.exists() or path.is_symlink():
        require(read_regular(path, 39936) == packed, "existing lookup differs; use a fresh build")
    else:
        publish(path, packed)


def roster():
    index = 0
    for rep in range(12):
        for step in range(2):
            order = (rep + step) % 2
            for width_step in range(3):
                width = (rep + step + width_step) % 3
                for metric_step in range(2):
                    metric = (rep + step + width_step + metric_step) % 2
                    for class_step in range(3):
                        comparison = (2 * rep + step + width_step + metric_step + class_step) % 3
                        for position, side in enumerate(SIDES):
                            j = rep % 4 if position < 2 else (position - 2) // 2
                            q = ((2 * (rep + 12 * j) + 1) * 1000000) // 96
                            yield dict(index=index, rep=rep, order=order, width=width,
                                       metric=metric, **{"class": comparison}, position=position,
                                       arm=PAIRS[comparison][side ^ order], q=q)
                            index += 1


def statistics(records):
    groups = {}
    require(len(records) == 4320, "all statistical records")
    for start in range(0, len(records), 10):
        panel = records[start:start + 10]
        first = panel[0]
        contrasts = []
        for j in range(4):
            by_side = {}
            for row in panel[2 + 2*j:4 + 2*j]:
                side = SIDES[row["position"]] ^ row["order"]
                by_side[side] = row["observation"][3] - row["observation"][2]
            require(set(by_side) == {0, 1} and min(by_side.values()) > 0, "complete paired times")
            contrasts.append(math.log(by_side[1]) - math.log(by_side[0]))
        key = (first["width"], first["metric"], first["class"], first["order"])
        groups.setdefault(key, []).append(math.fsum(contrasts) / 4)
    require(len(groups) == 36, "all separate width/metric/class/order groups")
    results, controls, treatments = [], [], []
    bound = math.log1p(.02)
    for (width, metric, comparison, order), values in sorted(groups.items()):
        estimate = A.confidence(values)
        result = dict(block_bytes=WIDTHS[width], metric=("repair", "encoder_lifecycle")[metric],
                      comparison=("BB", "CC", "BC")[comparison], order=order,
                      estimate=estimate, replicate_logs=values)
        if comparison < 2:
            result["control_pass"] = estimate["lower95_log"] > -bound and estimate["upper95_log"] < bound
            if not result["control_pass"]:
                controls.append([width, metric, comparison, order])
        else:
            limit = 0.0 if width == 0 else bound
            result["treatment_pass"] = estimate["upper95_log"] < limit
            result["upper_ratio_limit"] = math.exp(limit)
            if not result["treatment_pass"]:
                treatments.append([width, metric, comparison, order])
        results.append(result)
    return dict(outcome="CONTROL_FAIL" if controls else "FAIL" if treatments else "PASS",
                statistics=results, failed_controls=controls, failed_treatments=treatments,
                WH1_compared=False, public_WH2_compared=False, all_K_claimed=False,
                recovery_rate_claimed=False, production_promotion_claimed=False)


def multiply(a, b):
    product = 0
    for bit in range(8):
        if b & (1 << bit):
            product ^= a << bit
    for bit in range(14, 7, -1):
        if product & (1 << bit):
            product ^= 0x14d << (bit - 8)
    return product


def verify_preflight(header):
    exact(header["lookup_sha256"], LOOKUP_SHA, "lookup identity")
    lookup_address = integer(header["lookup_address"], 1, (1 << 64) - 39936)
    lookup = lookup_bytes()
    fixtures = header["fixtures"]
    require(type(fixtures) is list and len(fixtures) == 3, "three frozen widths")
    handles, spans, seen_handles = {}, [(lookup_address, lookup_address + 39936)], set()
    for width, (B, fixture) in enumerate(zip(WIDTHS, fixtures)):
        exact([fixture["width_index"], fixture["block_bytes"], fixture["message_bytes"], fixture["output_stride"]],
              [width, B, 6 * B, B + 128], "fixture dimensions")
        source = bytes((37 * i + i // 11) & 255 for i in range(6 * B))
        require(hex_bytes(fixture["source_hex"], 6 * B) == source and
                fixture["source_sha256"] == sha(source), "source formula/hash")
        source_address = integer(fixture["source_address"], 64, (1 << 64) - (6 * B + 64))
        require(source_address % 4096 == 64, "guarded source alignment")
        spans.append((source_address - 64, source_address + 6 * B + 64))
        addresses = fixture["output_addresses"]
        require(type(addresses) is list and len(addresses) == 2, "two output lanes")
        for address in addresses:
            integer(address, 1, (1 << 64) - 12 * (B + 128))
            require(address % 4096 == 0, "output base alignment")
            spans.append((address, address + 12 * (B + 128)))
        rows = hex_bytes(fixture["rows_hex"], 72)
        require(rows == lookup[:72], "independent low-ID row oracle")
        packets = bytearray()
        for packet_id in range(12):
            for byte in range(B):
                value = 0
                for column in range(6):
                    value ^= multiply(rows[packet_id * 6 + column], source[column * B + byte])
                packets.append(value)
        require(hex_bytes(fixture["packets_hex"], 12 * B) == bytes(packets) and
                fixture["packets_sha256"] == sha(packets), "independent full packet oracle")
        entries = fixture["handles"]
        require(type(entries) is list and len(entries) == 48, "retained handle count")
        expected_order = [(rep, lane, (rep + width + arm_step) % 2)
                          for rep in range(12) for arm_step in range(2) for lane in range(2)]
        for entry, (rep, lane, arm) in zip(entries, expected_order):
            exact([entry["rep"], entry["lane"], entry["arm"]], [rep, lane, arm], "allocation order")
            address = integer(entry["address"], 1, (1 << 64) - 96)
            require(address not in seen_handles, "retained handles cannot alias")
            seen_handles.add(address)
            spans.append((address, address + 96))
            handles[(width, rep, lane, arm)] = address
            exact(entry["roundtrip"], dict(create_code=0, feed_codes=[1, 1, 1, 1, 1, 0],
                  feed_count=6, recover_code=0, recover_bytes=6 * B), "first repair-only decode/recovery")
    for i, (begin, end) in enumerate(spans):
        require(all(end <= other_begin or other_end <= begin for other_begin, other_end in spans[:i]),
                "shared input/output/handle spans must be disjoint")
    return handles


def reduce_raw(raw):
    require(type(raw) is bytes and 0 < len(raw) <= RAW_CAP and raw.endswith(b"\n"), "bounded complete JSONL")
    lines = raw.splitlines()
    require(len(lines) == 4322, "header plus4320records plusfooter")
    header, *records, footer = [decode(line) for line in lines]
    exact([header["type"], footer["type"], header["protocol"], footer["protocol"], header["schema"]],
          ["header", "footer", PROTOCOL, PROTOCOL, PROTOCOL + ".raw.v1"], "framing")
    exact([footer["outcome"], footer["failure"]], ["COMPLETE", None], "complete correctness")
    exact(header["target_cpu"], 50, "target CPU")
    hex_bytes(header["claim_sha256"], 32)
    A.verify_identity(header["identity_before"])
    A.verify_identity(footer["identity_after"])
    exact(header["identity_before"], footer["identity_after"], "stable target identity")
    A.verify_runtime(header["runtime_before"])
    A.verify_runtime(footer["runtime_after"])
    exact(header["runtime_before"], footer["runtime_after"], "stable GF runtime")
    exact(header["clock_resolutions"], [1, 1], "clock resolutions")
    handles = verify_preflight(header)
    prelude = header["prelude"]
    exact([prelude["iterations"], prelude["seed"], prelude["final_state"]],
          [1 << 20, 0x9e3779b97f4a7c15, 0x43935dad1647741b], "single fixed prelude")
    A.clocks(prelude["observation"])
    start, cpu_start = integer(header["worker_start_ns"]), integer(header["worker_start_cpu_ns"])
    end, cpu_end = integer(footer["worker_end_ns"]), integer(footer["worker_end_cpu_ns"])
    require(start <= prelude["observation"][0] and cpu_start <= prelude["observation"][1], "prelude lifetime")
    require(0 < end - start < 15 * 10**9 and 0 <= cpu_end - cpu_start <= 14 * 10**9, "worker budgets")
    now, t0 = integer(footer["schedule_now_ns"]), integer(footer["t0"])
    require(now >= prelude["observation"][5] and t0 == (now // 2000000 + 2) * 2000000, "schedule epoch")
    previous = prelude["observation"]
    work_sum = wait_wall = wait_cpu = creates = encodes = frees = 0
    for expected, row in zip(roster(), records):
        exact({key: row[key] for key in expected}, expected, "exact callback roster")
        exact([row["type"], row["called"], row["checked"], row["work"]["complete"]],
              ["record", True, True, True], "complete callback")
        target = integer(row["target"])
        require(target == t0 + 2000000 * row["index"] + row["q"], "phase target")
        ready = integer(row["ready"])
        require(previous[5] <= ready <= target - 100000, "preparation margin")
        duration = A.clocks(row["observation"], previous)
        wait = row["wait"]
        require(type(wait) is list and len(wait) == 4, "wait shape")
        wm0, wc0, wm1, wc1 = [integer(x) for x in wait]
        require(ready <= wm0 <= target <= wm1 <= row["observation"][0] and
                previous[4] <= wc0 <= wc1 <= row["observation"][1], "wait ordering")
        work_sum += duration
        require(target <= row["observation"][2] <= target + 5000, "start lateness")
        wait_wall += wm1 - wm0
        wait_cpu += wc1 - wc0
        previous = row["observation"]
        work, metric = row["work"], row["metric"]
        B = WIDTHS[row["width"]]
        exact([work["create_calls"], work["encode_calls"], work["free_calls"]],
              [0, 384, 0] if metric == 0 else [64, 768, 64], "work calls")
        exact(work["create_code"], None if metric == 0 else 0, "creation result")
        count = 6 if metric == 0 else 12
        exact(work["codes"], [0] * count, "all packet status checks")
        exact(work["required"], [B] * count, "required lengths")
        exact(work["written"], [B] * count, "written lengths")
        addresses = work["addresses"]
        if metric == 0:
            exact(work["handle"], handles[(row["width"], row["rep"], row["order"], row["arm"])],
                  "same physical repair handle")
            exact(addresses, dict(count=0, sha256=sha(b""), min=None, max=None), "no fresh repair handles")
        else:
            exact(work["handle"], None, "no fabricated retained lifecycle handle")
            exact(addresses["count"], 64, "all fresh allocations recorded")
            hex_bytes(addresses["sha256"], 32)
            require(integer(addresses["min"], 1) <= integer(addresses["max"], 1), "allocation address range")
        creates += work["create_calls"]
        encodes += work["encode_calls"]
        frees += work["free_calls"]
    require(previous[5] <= end and previous[4] <= cpu_end and work_sum <= 10**9, "final work lifetime/cap")
    for key, value in dict(records=4320, callbacks=4320, checked_callbacks=4320,
            create_calls=138240, encode_calls=2488320, free_calls=138240,
            sum_work_ns=work_sum, sum_wait_wall_ns=wait_wall, sum_wait_cpu_ns=wait_cpu).items():
        exact(footer[key], value, "footer " + key)
    exact([creates, encodes, frees], [138240, 2488320, 138240], "full work ledger")
    exact(footer["preflight"], dict(encoder_creates=144, encoder_frees=144, encodes=1728,
          decoder_creates=144, feeds=864, recovers=144, decoder_frees=144), "untimed API ledger")
    require(wait_wall <= end - start and wait_cpu <= cpu_end - cpu_start, "conditioning included in lifetime")
    result = statistics(records)
    result.update(protocol=PROTOCOL, claim_sha256=header["claim_sha256"], raw_sha256=sha(raw),
                  callbacks=4320, create_calls=creates, encode_calls=encodes, free_calls=frees,
                  sum_work_ns=work_sum, conditioning=dict(
                      prelude_wall_ns=prelude["observation"][3] - prelude["observation"][2],
                      wait_wall_ns=wait_wall, wait_cpu_ns=wait_cpu,
                      worker_wall_ns=end - start, worker_cpu_ns=cpu_end - cpu_start))
    return result


def current_receipt(build, deadline):
    build = Path(build).resolve(strict=True)
    require(build.is_dir() and ROOT not in build.parents and build != ROOT, "external build directory")
    require(checked(["ninja", "-C", build, "-n"], deadline).decode().strip().endswith(
        "ninja: no work to do."), "build must be current")
    cache = read_regular(build / "CMakeCache.txt", 1024**2).decode()
    require("WH2_TINY_COST_SANITIZERS:BOOL=OFF" in cache, "native scientific build")
    require(sha(read_regular(build / "lookup.bin", 39936)) == LOOKUP_SHA, "build lookup identity")
    commands = decode(read_regular(build / "compile_commands.json", 2 * 1024**2))
    require(type(commands) is list and len(commands) == 12, "twelve exact translation units")
    expected_files = ["gf256.cpp", "Wh2ThueMorseNativeCodec.cpp", "Wh2ThueMorseTinyPayloadR0.cpp",
        "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp", "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp",
        "Wh2ThueMorseTinyPayloadCostR0.cpp", "Wh2RdpruTargetIdentityV2.cpp",
        "Wh2PublicBorrowedTargetIdentity.cpp", "Wh2FrozenTrace.cpp",
        "Wh2ThueMorseNativeCodecTest.cpp", "Wh2ThueMorseTinyPayloadR0CodecTest.cpp", "Wh2ThueMorseTinyPayloadR0Test.cpp"]
    expected_paths = [str(ROOT / (name if name == "gf256.cpp" else "bench/" + name)) for name in expected_files]
    exact(sorted(item["file"] for item in commands), sorted(expected_paths), "exact source-path roster")
    bridge_arms = []
    for item in commands:
        words = shlex.split(item["command"])
        require(Path(words[0]).resolve(strict=True) == Path("/usr/bin/g++").resolve(strict=True), "compiler identity")
        require(all(flag in words for flag in ("-std=c++11", "-O3", "-fno-lto")), "compiler flags")
        require("-fPIC" in words or "-fPIE" in words, "position-independent compile")
        arm_definitions = [word for word in words if word.startswith("-DWH2_TINY_COST_CANDIDATE")]
        if Path(item["file"]).name == "Wh2ThueMorseTinyPayloadCostR0Bridge.cpp":
            require(len(arm_definitions) == 1, "one bridge arm definition")
            bridge_arms.extend(arm_definitions)
        else:
            require(not arm_definitions, "bridge macro only belongs to wrappers")
        require(not any(flag in item["command"] for flag in ("WH_COUNT", "WIREHAIR_TESTING",
                    "WIREHAIR_V2_ENABLE_TEST_HOOKS", "-flto", "-fprofile", "-fsanitize", "-march=", "-mavx", "-mgfni")),
                "forbidden compiler instrumentation/ISA")
    exact(sorted(bridge_arms), ["-DWH2_TINY_COST_CANDIDATE=0", "-DWH2_TINY_COST_CANDIDATE=1"], "both bridge arms")
    dependencies = {Path(item["file"]) for item in commands}
    for line in checked(["ninja", "-C", build, "-t", "deps"], deadline).decode().splitlines():
        if line.startswith("    "):
            path = Path(line.strip())
            dependencies.add(path if path.is_absolute() else build / path)
    require(any(path.name == "gf256.h" for path in dependencies) and
            any(path.name == "Wh2ThueMorseNativeCodec.h" for path in dependencies), "actual include dependencies")
    # The configure-once manifest predates every object in this fresh build;
    # reconfiguration refuses changed native inputs instead of blessing stale objects.
    manifest = decode(read_regular(build / "native-inputs.json", 65536))
    require(type(manifest) is dict and manifest.get("schema") == PROTOCOL + ".native-inputs.v1", "native input manifest")
    require(type(manifest.get("sources")) is dict, "native input manifest source mapping")
    native_sources = {str(path.relative_to(ROOT)) for path in dependencies if ROOT in path.parents}
    exact(sorted(manifest["sources"]), sorted(native_sources), "all native inputs captured before build")
    for name, digest in manifest["sources"].items():
        require(pin(ROOT / name)["sha256"] == digest, "native source changed since configure: " + name)
    source_names = set(NEW_SOURCES) | {"AGENTS.md", "bench/Wh2AlignedIntermediateCostR0.py"}
    source_names.update(str(path.relative_to(ROOT)) for path in dependencies if ROOT in path.parents)
    source_names.update("bench/" + name for name in (
        "Wh2ThueMorseNativeDataR0.py", "Wh2ThueMorseRecoveryHistoryR0.py",
        "Wh2NoncommutingRadixRunR0.py", "Wh2ThueMorseR0.py", "Wh2NoncommutingRadixR0.py"))
    head = checked(["git", "rev-parse", "HEAD"], deadline).decode().strip()
    sources = []
    for name in sorted(source_names):
        time_left(deadline)
        value = pin(ROOT / name)
        require(value["sha256"] == sha(checked(["git", "cat-file", "blob", head + ":" + name], deadline)),
                "uncommitted source: " + name)
        sources.append(value)
    outputs = {build / name for name in (TARGET, "tiny_cost_baseline_test", "tiny_cost_candidate_test",
        "tiny_cost_kernel_test", "lookup.bin", "tiny-cost.map", "compile_commands.json", "CMakeCache.txt",
        "build.ninja", "CMakeFiles/rules.ninja", "native-inputs.json")}
    outputs.update(build.rglob("*.o"))
    outputs.update(build.rglob("*.a"))
    installed = {Path(program).resolve(strict=True) for program in (
        "/usr/bin/cmake", "/usr/bin/ninja", "/usr/bin/g++", "/usr/bin/as", "/usr/bin/ld", "/usr/bin/git",
        "/usr/bin/ldd", "/bin/bash", "/usr/bin/python3")}
    for name in ("cc1plus", "collect2"):
        installed.add(Path(checked(["g++", "-print-prog-name=" + name], deadline).decode().strip()).resolve(strict=True))
    linked = A.runtime_libraries(checked(["ldd", build / TARGET], deadline).decode())
    installed.update(Path(entry["path"]) for entry in linked if entry["path"] is not None)
    dependency_pins = []
    for path in sorted(dependencies):
        time_left(deadline)
        dependency_pins.append(pin(path, installed=build not in path.parents and ROOT not in path.parents))
    return dict(protocol=PROTOCOL, schema=PROTOCOL + ".claim.v1", build=str(build), source_head=head,
                sources=sources, dependencies=dependency_pins, outputs=[pin(path) for path in sorted(outputs)],
                installed=[pin(path, installed=True) for path in sorted(installed)],
                controller=pin(Path(__file__).resolve()), interpreter=pin(Path(sys.executable), installed=True),
                compiler_version=checked(["g++", "--version"], deadline).decode(), runtime_ldd=linked,
                target_cpu=50, raw_cap=RAW_CAP, stderr_cap=ERR_CAP, worker_wall_seconds=15,
                worker_cpu_seconds=14, controller_wall_seconds=20)


def controller_run(receipt_path):
    start = time.monotonic()
    deadline = start + 20
    receipt_bytes = read_regular(receipt_path, 1024**2)
    receipt = decode(receipt_bytes)
    require(receipt_bytes == canonical(receipt), "canonical launch receipt")
    exact(receipt, current_receipt(receipt["build"], deadline), "prelaunch closure")
    os.mkdir(str(OUTPUT), 0o700)
    publish(OUTPUT / "CLAIM.json", receipt_bytes)
    raw, stderr, returncode, failure = b"", b"", None, None
    analysis = dict(outcome="INVALID", failure="worker not started")
    try:
        raw, stderr, returncode, failure = A.capture(Path(receipt["build"]) / TARGET, sha(receipt_bytes), deadline,
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
                   returncode=returncode, elapsed_seconds=time.monotonic() - start,
                   raw_sha256=sha(raw), production_promotion_claimed=False, all_K_claimed=False)
    publish(OUTPUT / "summary.json", canonical(summary))
    members = [pin(OUTPUT / name) for name in ("CLAIM.json", "raw.jsonl", "stderr.txt", "analysis.json", "summary.json")]
    require(sum(item["bytes"] for item in members) < TOTAL_CAP - 65536, "bundle cap")
    publish(OUTPUT / "COMPLETE.json", canonical(dict(protocol=PROTOCOL, files=members, outcome=analysis["outcome"])))
    print(canonical(summary).decode(), end="")


def run(receipt_path):
    """Same reviewed 25s process-group observer as .68, launching THIS module."""
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
