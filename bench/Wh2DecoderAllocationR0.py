#!/usr/bin/env python3
"""One bounded allocation diagnostic; never repeats or scores codec timing."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time

SPEC = importlib.util.spec_from_file_location(
    "allocation_common", Path(__file__).with_name("Wh2AlignedIntermediateCostR0.py"))
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("/var/tmp/wh2-decoder-allocation-r0")
OLD = Path("/var/tmp/wh2-thue-decoder-cost-r0")
PROTOCOL = "wirehair.wh2.decoder-allocation-r0"
SOURCE_NAMES = ("bench/Wh2DecoderAllocationR0.cpp", "bench/Wh2DecoderAllocationR0.py")
REUSED = Path("/tmp/wh2-decoder-cost-neutral.utKXCA/native")


def pin_input(path):
    # Installed tools/headers may legitimately have multiple hard links. Keep
    # the stricter single-link rule for repository and temporary artifacts.
    path = Path(path)
    owned = any(base == path or base in path.parents for base in (ROOT, Path("/tmp"), Path("/var/tmp")))
    return A.pin(path, installed=not owned)


def build_commands(build):
    return [
        ["/usr/bin/g++", "-std=c++11", "-O3", "-g1", "-fPIC", "-fno-lto",
         "-Wall", "-Wextra", "-Wpedantic", "-Werror", "-I", str(ROOT / "include"),
         "-MD", "-MF", str(build / "worker.d"), "-c", str(ROOT / SOURCE_NAMES[0]),
         "-o", str(build / "worker.o")],
        ["/usr/bin/g++", "-fno-lto", "-no-pie", str(build / "worker.o"),
         str(REUSED / "codec-neutral/CMakeFiles/tiny_gfni_runtime.dir/home/catid/wirehair/bench/Wh2ThueMorseTinyGfniR0.cpp.o"),
         str(REUSED / "production/libwirehair.a"), "-pthread",
         "-Wl," + ",".join("--wrap=" + name for name in
                           ("_Znwm", "_Znam", "_ZdlPv", "_ZdaPv", "_ZdlPvm", "_ZdaPvm",
                            "malloc", "calloc", "realloc", "free")),
         "-Wl,-Map," + str(build / "worker.map"), "-o", str(build / "worker")]]


def build_worker(build):
    build = build.resolve(strict=True)
    A.require(build.parent == Path("/tmp") and build.name.startswith("wh2-decoder-allocation-build."), "owned build directory")
    commands = build_commands(build)
    for command in commands:
        subprocess.run(command, check=True, timeout=30, cwd=ROOT)
    record = A.canonical(commands)
    if (build / "build.json").exists():
        A.exact(A.read_regular(build / "build.json", 65536), record, "unchanged build commands")
    else:
        A.publish(build / "build.json", record)


def verify(raw, claim, original=None):
    A.require(0 < len(raw) <= 2*1024**2 and raw.endswith(b"\n"), "bounded complete worker output")
    rows = [A.decode(line) for line in raw.splitlines()]
    header, footer = rows[0], rows[-1]
    A.require(header["protocol"] == PROTOCOL and header["type"] == "header" and
              header["claim_sha256"] == claim and header["speed_claimed"] is False, "header")
    A.integer(header["pid"], 1)
    A.require(footer["type"] == "footer" and footer["outcome"] == "COMPLETE", "footer")
    A.require(all(r["type"] in ("fixture", "phase", "allocation") for r in rows[1:-1]), "row type")
    fixtures = [r for r in rows if r["type"] == "fixture"]
    phases = [r for r in rows if r["type"] == "phase"]
    events = [r for r in rows if r["type"] == "allocation"]
    A.require(len(fixtures) == 3 and 162 <= len(phases) <= 270 and len(events) <= 8192, "fixed cohort/caps")
    A.exact((len(phases), len(events)), (footer["phases"], footer["events"]), "footer counts")
    # This compares deterministic fixture bytes, never replays the spent codec.
    if original is None:
        with (OLD / "raw.jsonl").open("rb") as stream:
            original = A.decode(stream.readline(1024*1024))
    A.exact(len(original["fixtures"]), 3, "original fixture count")
    for f, old, width in zip(fixtures, original["fixtures"], (2, 64, 1280)):
        A.exact(f["width"], width, "fixture width")
        profile = next(s["profile_hex"] for s in old["handles"] if s["arm"] == 3)
        A.exact((f["profile_hex"], f["packets_hex"]), (profile, old["packets_hex"][3]), "exact .76 public corpus")
    p = 0
    def consume(operation, width, family, cycle, packet):
        nonlocal p
        r = phases[p]; p += 1
        A.exact((r["index"], r["operation"], r["width"], r["family"], r["cycle"], r["packet"]),
                (p, operation, width, family, cycle, packet), "exact phase chronology")
        return r
    for width in (2, 64, 1280):
        for family in range(2):
            for cycle in range(3):
                A.exact(consume("create", width, family, cycle, 0)["result"], 0, "create success")
                for step in range(12):
                    packet = (6+step if family == 0 else 0xffffffff-2*step) if step < 6 else step-6
                    r = consume("feed", width, family, cycle, packet)
                    A.integer(r["result"], 0, 1)
                    if r["result"] == 0:
                        A.require(step >= 5, "no premature success")
                        break
                    A.exact(r["result"], 1, "need more")
                else:
                    raise ValueError("no decode endpoint")
                A.exact(consume("recover", width, family, cycle, 0)["result"], 0, "recover success")
                A.exact(consume("free", width, family, cycle, 0)["result"], 0, "free success")
    A.exact(p, len(phases), "no extra phases")
    cursor = 0
    counters = [0]*4
    live = {}
    live_bytes = peak_bytes = 0
    allocations = 0
    last_counters = [0]*4
    for r in phases:
        A.exact(r["first_event"], cursor, "event cursor")
        end = A.integer(r["end_event"], cursor, len(events))
        for i in range(cursor, end):
            e = events[i]
            A.exact((e["index"], e["phase"]), (i, r["index"]), "event attribution")
            A.integer(e["kind"], 0, 9)
            for key in ("pointer", "prior", "caller", "bytes"):
                A.integer(e[key], 0, (1 << 64)-1)
            A.require(e["caller"] != 0, "allocation callsite")
            kind, pointer, prior, size = (e[k] for k in ("kind", "pointer", "prior", "bytes"))
            A.require(kind == 8 or prior == 0, "only realloc has a prior pointer")
            if kind == 8 and prior:
                A.require(prior in live and live[prior][1] == 6, "realloc of observed C allocation")
            if kind in (2, 3, 4, 5, 9) or (kind == 8 and prior and (pointer or size == 0)):
                freed = prior if kind == 8 else pointer
                if freed:
                    A.require(freed in live, "free of unobserved allocation")
                    old_size, old_kind = live.pop(freed)
                    expected_kind = 0 if kind in (2, 4) else 1 if kind in (3, 5) else 6
                    A.exact(old_kind, expected_kind, "matching allocation/free family")
                    live_bytes -= old_size
            if kind in (0, 1, 6, 7, 8) and pointer:
                A.require(pointer not in live, "overlapping live allocation")
                live[pointer] = (size, kind if kind < 2 else 6)
                live_bytes += size
                peak_bytes = max(peak_bytes, live_bytes)
                allocations += 1
        cursor = end
        A.require(len(r["counters_before"]) == len(r["counters_after"]) == 4, "counter shape")
        for i, (before, after) in enumerate(zip(r["counters_before"], r["counters_after"])):
            A.integer(before, last_counters[i])
            counters[i] += A.integer(after, A.integer(before)) - before
            last_counters[i] = after
        if r["operation"] == "free":
            A.require(not live and live_bytes == 0, "complete observed lifecycle ownership")
    A.exact(cursor, len(events), "all events attributed")
    return dict(outcome="COMPLETE", phases=len(phases), events=len(events), counters=counters,
                allocations=allocations, peak_live_requested_bytes=peak_bytes,
                speed_claimed=False, reproduces_old_heap_layout_claimed=False)


def pins_current(receipt):
    A.exact(receipt["protocol"], PROTOCOL, "receipt protocol")
    A.exact(receipt["build_commands"], build_commands(Path(receipt["build"])), "exact diagnostic link/compile")
    A.exact(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            receipt["source_head"], "source HEAD")
    for pin in receipt["pins"]:
        A.exact(pin_input(Path(pin["path"])), pin, "unchanged input/build/helper")


def make_receipt(build):
    build = build.resolve(strict=True)
    A.require(ROOT not in build.parents and build != ROOT, "external build")
    old = A.decode(A.read_regular(OLD / "CLAIM.json", 1024*1024))
    # Every reused original source/object/dependency retains its original bytes.
    old_pins = old["sources"] + old["dependencies"] + old["outputs"] + old["installed"]
    paths = set()
    for pin in old_pins:
        path = Path(pin["path"])
        A.exact(hashlib.sha256(path.read_bytes()).hexdigest(), pin["sha256"], "unchanged .76 input")
        paths.add(path.resolve(strict=True))
    paths.update((build / name) for name in ("worker", "worker.o", "worker.d", "worker.map", "build.json"))
    commands = A.decode(A.read_regular(build / "build.json", 65536))
    A.exact(commands, build_commands(build), "recorded diagnostic build")
    dependency_text = A.read_regular(build / "worker.d", 1024*1024).decode().replace("\\\n", " ")
    paths.update(Path(name).resolve(strict=True) for name in shlex.split(dependency_text.split(":", 1)[1]))
    paths.update((ROOT / name) for name in SOURCE_NAMES)
    paths.update(Path(name).resolve(strict=True) for name in ("/usr/bin/strace", "/usr/bin/prlimit", sys.executable))
    paths.update((OLD / "CLAIM.json", OLD / "raw.jsonl"))
    for executable in (build / "worker", Path("/usr/bin/strace"), Path("/usr/bin/prlimit"), Path(sys.executable)):
        linked = subprocess.check_output(["/usr/bin/ldd", str(executable)], text=True, timeout=10)
        paths.update(Path(word).resolve(strict=True) for word in linked.split() if word.startswith("/"))
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    for name in SOURCE_NAMES:
        A.exact((ROOT / name).read_bytes(), subprocess.check_output(["git", "cat-file", "blob", head+":"+name], cwd=ROOT),
                "committed diagnostic source")
    return dict(protocol=PROTOCOL, source_head=head, build=str(build), build_commands=commands,
                pins=[pin_input(path) for path in sorted(paths)])


def run(receipt_path):
    receipt_bytes = A.read_regular(receipt_path, 1024*1024); receipt = A.decode(receipt_bytes)
    A.exact(receipt_bytes, A.canonical(receipt), "canonical receipt")
    pins_current(receipt)
    os.mkdir(str(OUTPUT), 0o700)
    A.publish(OUTPUT / "CLAIM.json", receipt_bytes)
    # Ten CPU seconds each for the tracer and its sole codec child: the
    # prospective twenty-second group budget is not a per-process twenty.
    command = ["/usr/bin/prlimit", "--fsize=4194304", "--cpu=10", "--as=268435456", "--core=0",
        "/usr/bin/strace", "-f", "-qq", "-k", "--stack-trace-frame-limit=12", "-s", "80",
        "-e", "trace=brk,mmap,munmap,mremap,madvise,write", "-o", str(OUTPUT / "syscalls.txt"),
        str(Path(receipt["build"]) / "worker"), "--worker", A.sha(receipt_bytes)]
    child = None; raw = error = b""; failure = None; start = time.monotonic()
    try:
        child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True, close_fds=True)
        try:
            raw, error = child.communicate(timeout=35)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            raw, error = child.communicate()
            raise ValueError("35-second observer deadline")
        A.require(child.returncode == 0 and len(error) <= 65536, "worker/tracer result or stderr cap")
        pins_current(receipt)
        analysis = verify(raw, A.sha(receipt_bytes))
        expected_markers = b"".join(("WH2_ALLOC_PHASE %d %s\n" % (i, side)).encode()
            for i in range(1, analysis["phases"]+1) for side in ("BEGIN", "END"))
        A.exact(error, expected_markers, "complete stderr phase markers")
        A.require(0 < (OUTPUT / "syscalls.txt").stat().st_size <= 4*1024**2, "syscall artifact cap")
        A.require(time.monotonic()-start < 35, "35-second observer deadline")
    except Exception as exc:
        failure = str(exc); analysis = dict(outcome="INVALID", failure=failure, speed_claimed=False)
    finally:
        if child is not None and child.returncode is None:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
    A.publish(OUTPUT / "raw.jsonl", raw)
    A.publish(OUTPUT / "stderr.txt", error)
    if not (OUTPUT / "syscalls.txt").exists():
        A.publish(OUTPUT / "syscalls.txt", b"")
    os.chmod(str(OUTPUT / "syscalls.txt"), 0o400)
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-start,
                    returncode=None if child is None else child.returncode)
    A.publish(OUTPUT / "analysis.json", A.canonical(analysis))
    members = [A.pin(OUTPUT / name) for name in ("CLAIM.json", "raw.jsonl", "stderr.txt", "syscalls.txt", "analysis.json")]
    A.require(sum(p["bytes"] for p in members) < 8*1024**2-65536, "bundle cap")
    A.publish(OUTPUT / "COMPLETE.json", A.canonical(dict(protocol=PROTOCOL, outcome=analysis["outcome"], files=members)))
    print(json.dumps(analysis, sort_keys=True))


def selftest():
    """Synthetic records and mocked children only: no selected codec calls."""
    import copy
    import contextlib
    import io
    import tempfile
    from unittest import mock
    with mock.patch.object(A, "pin", return_value={}) as pinned:
        for path, installed in ((Path("/usr/bin/bash"), True), (ROOT / SOURCE_NAMES[0], False),
                                (Path("/tmp/worker"), False), (Path("/var/tmp/raw"), False)):
            pin_input(path)
            pinned.assert_called_with(path, installed=installed)
    claim = "a" * 64
    rows = [dict(type="header", protocol=PROTOCOL, claim_sha256=claim, pid=1, speed_claimed=False)]
    original = dict(fixtures=[])
    events = []
    phase = 0
    for width in (2, 64, 1280):
        profile, packets = "00"*32, "00"*(18*width)
        rows.append(dict(type="fixture", width=width, profile_hex=profile, packets_hex=packets))
        original["fixtures"].append(dict(handles=[dict(arm=3, profile_hex=profile)],
                                         packets_hex=[None, None, None, packets]))
        for family in range(2):
            for cycle in range(3):
                pointer = 1000 + phase
                operations = [("create", 0, 0)] + [
                    ("feed", 6+i if family == 0 else 0xffffffff-2*i, int(i != 5)) for i in range(6)
                ] + [("recover", 0, 0), ("free", 0, 0)]
                for operation, packet, result in operations:
                    phase += 1
                    first = len(events)
                    if operation in ("create", "free"):
                        events.append(dict(type="allocation", index=first, phase=phase,
                            kind=0 if operation == "create" else 2, bytes=17 if operation == "create" else 0,
                            pointer=pointer, prior=0, caller=0x400000))
                    rows.append(dict(type="phase", index=phase, operation=operation, width=width,
                        family=family, cycle=cycle, packet=packet, result=result, first_event=first,
                        end_event=len(events), counters_before=[0]*4, counters_after=[0]*4))
    rows += events + [dict(type="footer", outcome="COMPLETE", phases=phase, events=len(events))]
    serialize = lambda value: b"".join(A.canonical(r) for r in value)
    raw = serialize(rows)
    result = verify(raw, claim, original)
    A.exact((result["phases"], result["events"], result["allocations"], result["peak_live_requested_bytes"]),
            (162, 36, 18, 17), "synthetic independent ledger")
    rejected = 0
    def reject(payload, frozen=original):
        nonlocal rejected
        try:
            verify(payload, claim, frozen)
        except (ValueError, KeyError, IndexError, TypeError, StopIteration):
            rejected += 1
            return
        raise AssertionError("corrupt synthetic record accepted")
    mutations = [
        (0, "claim_sha256", "b"*64), (0, "pid", False), (-1, "phases", 163),
        (1, "packets_hex", "ff"), (2, "index", True), (3, "packet", 42),
        (8, "result", False), (2, "end_event", 0),
        (2, "counters_before", [1, 0, 0, 0]), (2, "counters_after", [-1, 0, 0, 0]),
        (-2, "pointer", 999999), (-2, "kind", 9), (-2, "caller", 0),
        (-2, "prior", 1), (-2, "phase", 1), (-2, "bytes", True),
    ]
    for index, key, value in mutations:
        damaged = copy.deepcopy(rows); damaged[index][key] = value
        reject(serialize(damaged))
    reject(raw[:-1]); reject(raw[:100]); reject(b"x"*(2*1024**2+1))
    reject(serialize(rows[:-2] + rows[-1:]))
    reject(raw, dict(fixtures=original["fixtures"][:2]))
    # Both terminal paths publish retained prefixes and never launch a codec.
    for timed_out in (False, True):
        with tempfile.TemporaryDirectory(prefix="wh2-allocation-neutral.") as tmp:
            base = Path(tmp); output = base / "artifact"
            receipt_path = base / "receipt.json"
            A.publish(receipt_path, A.canonical(dict(build=str(base))))
            fake = mock.Mock(pid=123456789, returncode=None if timed_out else 1)
            fake.communicate.side_effect = [subprocess.TimeoutExpired("neutral", 35), (b"prefix", b"error")] if timed_out else [(b"prefix", b"error")]
            with mock.patch.dict(globals(), OUTPUT=output), mock.patch.object(sys.modules[__name__], "pins_current"), \
                    mock.patch.object(subprocess, "Popen", return_value=fake), \
                    mock.patch.object(os, "killpg") as killed, contextlib.redirect_stdout(io.StringIO()):
                run(receipt_path)
                A.exact(A.read_regular(output / "raw.jsonl", 1024), b"prefix", "retained failure prefix")
                A.exact(A.decode(A.read_regular(output / "analysis.json", 65536))["outcome"], "INVALID", "terminal failure")
                if timed_out:
                    A.require(killed.call_count >= 1, "timeout cleanup")
                else:
                    killed.assert_not_called()
                try:
                    run(receipt_path)
                except FileExistsError:
                    pass
                else:
                    raise AssertionError("spent namespace accepted")
    print("PASS synthetic ledger, %d corruptions, failure/timeout retention and spent namespace; no codec work" % rejected)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "selftest", "receipt", "run"))
    parser.add_argument("path", nargs="?", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    args = parser.parse_args()
    if args.command == "selftest":
        A.require(args.path is None and args.output is None, "selftest has no paths")
        selftest()
    elif args.command == "build":
        A.require(args.path is not None and args.output is None, "build takes one path")
        build_worker(args.path)
    elif args.command == "receipt":
        A.require(args.path is not None, "build path required")
        A.require(args.output is not None, "receipt output required")
        A.publish(args.output, A.canonical(make_receipt(args.path)))
    else:
        A.require(args.path is not None, "receipt path required")
        A.require(args.output is None, "run takes only a receipt")
        run(args.path)
