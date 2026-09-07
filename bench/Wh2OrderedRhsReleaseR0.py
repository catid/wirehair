#!/usr/bin/env python3
"""One fixed B,C,C,B mechanism screen. No timing or footprint claims."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

SPEC = importlib.util.spec_from_file_location("release_common", Path(__file__).with_name("Wh2AlignedIntermediateCostR0.py"))
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("/var/tmp/wh2-ordered-rhs-release-r0")
OLD = Path("/var/tmp/wh2-thue-decoder-cost-r0")
PROTOCOL = "wirehair.wh2.ordered-rhs-release-r0"
ARMS = ("baseline", "candidate", "candidate", "baseline")
SOURCES = ("bench/Wh2OrderedRhsReleaseR0.cpp", "bench/Wh2OrderedRhsReleaseR0.py",
           "bench/Wh2OrderedRhsReleaseR0/CMakeLists.txt", "bench/test_Wh2OrderedRhsReleaseR0.py",
           "bench/Wh2AlignedIntermediateCostR0.py", "CMakeLists.txt")
ENVIRONMENT_KEYS = ("MALLOC_TRIM_THRESHOLD_", "MALLOC_MMAP_THRESHOLD_", "MALLOC_TOP_PAD_",
                    "MALLOC_PERTURB_", "GLIBC_TUNABLES", "LD_PRELOAD", "LD_LIBRARY_PATH")


def pin_input(path):
    path = Path(path)
    owned = any(base == path or base in path.parents for base in (ROOT, Path("/tmp"), Path("/var/tmp")))
    return A.pin(path, installed=not owned)


def current(receipt):
    A.exact(receipt["protocol"], PROTOCOL, "receipt protocol")
    A.exact(receipt["arms"], list(ARMS), "fixed process order")
    A.exact(receipt["allocator_environment"], {k: os.environ.get(k) for k in ENVIRONMENT_KEYS}, "unchanged allocator environment")
    A.exact(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), receipt["head"], "source HEAD")
    for pin in receipt["pins"]:
        A.exact(pin_input(Path(pin["path"])), pin, "unchanged input")


def receipt(build):
    build = build.resolve(strict=True)
    A.require(ROOT not in build.parents and build != ROOT, "external build")
    cache = (build / "CMakeCache.txt").read_text()
    for name in ("WH2_RELEASE_SANITIZERS", "WH2_RELEASE_SCALAR", "MARCH_NATIVE"):
        A.require(name+":BOOL=OFF\n" in cache, "native scientific build")
    paths = {ROOT / name for name in SOURCES}
    paths.update(build / name for name in ("CMakeCache.txt", "build.ninja", "compile_commands.json", "WirehairV2SolveOrdered.cpp",
        "baseline.map", "candidate.map", "release_baseline_worker", "release_candidate_worker"))
    paths.update(build.rglob("*.o")); paths.update(build.rglob("*.a"))
    deps = subprocess.check_output(["/usr/bin/ninja", "-C", str(build), "-t", "deps"], text=True, timeout=10)
    A.require(len(deps) < 4*1024**2, "dependency output cap")
    paths.update(Path(line.strip()).resolve(strict=True) for line in deps.splitlines() if line.startswith("    "))
    paths.update((OLD / "CLAIM.json", OLD / "raw.jsonl"))
    programs = [Path("/usr/bin") / name for name in ("cmake", "ninja", "strace", "prlimit", "c++", "as", "ld")]
    programs.append(Path(sys.executable))
    paths.update(p.resolve(strict=True) for p in programs)
    for program in (Path(sys.executable), Path("/usr/bin/strace"), Path("/usr/bin/prlimit"), build / "release_baseline_worker", build / "release_candidate_worker"):
        linked = subprocess.check_output(["/usr/bin/ldd", str(program)], text=True, timeout=10)
        paths.update(Path(word).resolve(strict=True) for word in linked.split() if word.startswith("/"))
    for arm in ("baseline", "candidate"):
        mapping = (build / (arm+".map")).read_text()
        A.require("libwirehair.a(WirehairV2Solve.cpp.o)" not in mapping, "no archive solve extraction")
        A.require("CMakeFiles/release_"+arm+".dir/" in mapping, "explicit solve object")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    paths = {p.resolve(strict=True) for p in paths}
    for path in sorted(paths):
        if ROOT in path.parents:
            name = str(path.relative_to(ROOT))
            A.exact(path.read_bytes(), subprocess.check_output(["git", "cat-file", "blob", head+":"+name], cwd=ROOT), "committed source")
    return dict(protocol=PROTOCOL, arms=list(ARMS), head=head, build=str(build),
                allocator_environment={k: os.environ.get(k) for k in ENVIRONMENT_KEYS},
                pins=[pin_input(p) for p in sorted(paths)])


def verify(raw, stderr, trace, claim, original=None):
    A.require(0 < len(raw) <= 256*1024 and raw.endswith(b"\n"), "complete bounded output")
    A.require(len(stderr) <= 65536 and 0 < len(trace) <= 4*1024**2, "trace/stderr caps")
    rows = [A.decode(line) for line in raw.splitlines()]
    header, footer = rows[0], rows[-1]
    A.require(header["type"] == "header" and header["protocol"] == PROTOCOL and
              header["claim_sha256"] == claim and header["speed_claimed"] is False, "header")
    pid = A.integer(header["pid"], 1)
    A.require(footer["type"] == "footer" and footer["outcome"] == "COMPLETE", "complete footer")
    A.require(all(r["type"] in ("fixture", "phase") for r in rows[1:-1]), "record types")
    fixtures = [r for r in rows if r["type"] == "fixture"]
    phases = [r for r in rows if r["type"] == "phase"]
    A.require(len(fixtures) == 3 and 162 <= len(phases) <= 270, "fixed cohort/cap")
    A.exact(footer["phases"], len(phases), "phase count")
    if original is None:
        with (OLD / "raw.jsonl").open("rb") as stream:
            original = A.decode(stream.readline(1024*1024))
    A.exact(len(original["fixtures"]), 3, "old fixtures")
    for f, old, width in zip(fixtures, original["fixtures"], (2, 64, 1280)):
        profile = next(h["profile_hex"] for h in old["handles"] if h["arm"] == 3)
        A.exact((f["width"], f["profile_hex"], f["packets_hex"]), (width, profile, old["packets_hex"][3]), "exact public fixture")
    index = 0
    chronology = []
    def consume(op, width, family, cycle, packet):
        nonlocal index
        p = phases[index]; index += 1
        want = (index, op, width, family, cycle, packet)
        A.exact(tuple(p[k] for k in ("index", "operation", "width", "family", "cycle", "packet")), want, "phase chronology")
        result = A.integer(p["result"], 0, 1)
        chronology.append(list(want)+[result])
        return result
    for width in (2, 64, 1280):
        for family in range(2):
            for cycle in range(3):
                A.exact(consume("create", width, family, cycle, 0), 0, "create success")
                for step in range(12):
                    packet = (6+step if family == 0 else 0xffffffff-2*step) if step < 6 else step-6
                    if consume("feed", width, family, cycle, packet) == 0:
                        A.require(step >= 5, "no premature success")
                        break
                else:
                    raise ValueError("no endpoint")
                A.exact(consume("recover", width, family, cycle, 0), 0, "recover success")
                A.exact(consume("free", width, family, cycle, 0), 0, "free success")
    A.exact(index, len(phases), "no extra phases")
    counters = [0]*4; last = [0]*4
    for p in phases:
        A.require(len(p["before"]) == len(p["after"]) == 4, "counter shape")
        for i, (before, after) in enumerate(zip(p["before"], p["after"])):
            A.integer(before, last[i]); A.integer(after, before)
            counters[i] += after-before; last[i] = after
    expected_markers = [(i, side) for i in range(1, len(phases)+1) for side in ("BEGIN", "END")]
    A.exact(stderr, b"".join(("WH2_RELEASE_PHASE %d %s\n" % x).encode() for x in expected_markers), "stderr markers")
    active = None; found = []; program_break = None
    calls = [[] for _ in phases]
    for line in trace.decode().splitlines():
        marker = re.fullmatch(r'(\d+) +write\(2, "WH2_RELEASE_PHASE (\d+) (BEGIN|END)\\n", (\d+)\) += (\d+)', line)
        if marker:
            A.exact(int(marker[1]), pid, "marker PID")
            i, side = int(marker[2]), marker[3]
            A.require(int(marker[4]) == int(marker[5]) == len("WH2_RELEASE_PHASE %d %s\n" % (i, side)), "complete marker write")
            found.append((i, side))
            if side == "BEGIN":
                A.require(active is None and 1 <= i <= len(phases), "marker begin")
                active = i
            else:
                A.exact(active, i, "marker end"); active = None
            continue
        call = re.match(r'(\d+) +(brk|mmap|munmap|mremap|madvise)\(', line)
        if not call:
            continue
        A.exact(int(call[1]), pid, "syscall PID")
        shrink = False
        if call[2] == "brk":
            brk = re.fullmatch(r'\d+ +brk\((NULL|0x[0-9a-f]+)\) += (0x[0-9a-f]+)', line)
            A.require(brk is not None, "successful brk syntax")
            next_break = int(brk[2], 16)
            if brk[1] != "NULL":
                A.exact(int(brk[1], 16), next_break, "successful brk")
                A.require(program_break is not None, "known previous break")
                shrink = next_break < program_break
            program_break = next_break
        if active is not None:
            calls[active-1].append(dict(syscall=call[2], shrink=shrink))
    A.require(active is None, "closed final phase")
    A.exact(found, expected_markers, "complete ordered syscall markers")
    cells = []
    for family in range(2):
        for cycle in range(3):
            chosen = [(p, calls[i]) for i, p in enumerate(phases) if (p["width"], p["family"], p["cycle"]) == (1280, family, cycle)]
            cells.append(dict(family=family, cycle=cycle,
                solve_shrinks=sum(int(c["shrink"]) for p, cs in chosen if p["operation"] == "feed" for c in cs),
                lifecycle_vm_calls=sum(len(cs) for _, cs in chosen)))
    return dict(cells=cells, chronology=chronology, counters=counters, phases=len(phases), speed_claimed=False)


def decide(results):
    A.exact(len(results), 4, "four processes")
    for r in results[1:]:
        A.exact(r["chronology"], results[0]["chronology"], "unchanged decoder endpoints")
    b, c = results[0]["cells"], results[1]["cells"]
    if b != results[3]["cells"] or c != results[2]["cells"] or any(x["solve_shrinks"] < 2 for x in b):
        return "INCONCLUSIVE"
    for family in range(2):
        bb = [x for x in b if x["family"] == family]; cc = [x for x in c if x["family"] == family]
        if sum(x["solve_shrinks"] for x in cc) >= sum(x["solve_shrinks"] for x in bb) or \
                sum(x["lifecycle_vm_calls"] for x in cc) > sum(x["lifecycle_vm_calls"] for x in bb):
            return "FAIL"
    return "PASS"


def capture(command):
    child = None; raw = error = b""; failure = None
    try:
        child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 start_new_session=True, close_fds=True)
        try:
            raw, error = child.communicate(timeout=35)
        except subprocess.TimeoutExpired:
            try: os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            raw, error = child.communicate()
            raise ValueError("35-second process observer")
        A.require(child.returncode == 0, "worker/tracer exit")
    except Exception as exc:
        failure = str(exc)
    finally:
        if child is not None and child.returncode is None:
            try: os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            child.wait()
    return raw, error, failure


def run(path):
    frozen = A.read_regular(path, 1024*1024); claimed = A.decode(frozen)
    A.exact(frozen, A.canonical(claimed), "canonical receipt"); current(claimed)
    os.mkdir(str(OUTPUT), 0o700); A.publish(OUTPUT / "CLAIM.json", frozen)
    results = []; failure = None; outcome = "INVALID"; start = time.monotonic()
    try:
        for i, arm in enumerate(ARMS):
            A.require(time.monotonic()-start < 150, "whole observer")
            trace_path = OUTPUT / (str(i)+"-trace.txt")
            command = ["/usr/bin/prlimit", "--cpu=10", "--as=268435456", "--fsize=4194304", "--core=0",
                "/usr/bin/strace", "-f", "-qq", "-k", "--stack-trace-frame-limit=12", "-s", "80",
                "-e", "trace=brk,mmap,munmap,mremap,madvise,write", "-o", str(trace_path),
                str(Path(claimed["build"]) / ("release_"+arm+"_worker")), "--worker", A.sha(frozen)]
            process_start = time.monotonic()
            raw, error, worker_failure = capture(command)
            A.publish(OUTPUT / (str(i)+"-raw.jsonl"), raw); A.publish(OUTPUT / (str(i)+"-stderr.txt"), error)
            if not trace_path.exists(): A.publish(trace_path, b"")
            os.chmod(str(trace_path), 0o400)
            A.require(worker_failure is None, str(worker_failure))
            result = verify(raw, error, A.read_regular(trace_path, 4*1024**2), A.sha(frozen))
            A.require(time.monotonic()-process_start < 35, "process observer")
            result["arm"] = arm; results.append(result)
        current(claimed)
        A.require(time.monotonic()-start < 150, "whole observer")
        outcome = decide(results)
    except Exception as exc:
        failure = str(exc)
    analysis = dict(protocol=PROTOCOL, outcome=outcome, failure=failure, results=results,
                    elapsed_seconds=time.monotonic()-start, speed_claimed=False)
    A.publish(OUTPUT / "analysis.json", A.canonical(analysis))
    members = [A.pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p["bytes"] for p in members) < 24*1024**2-65536, "bundle cap")
    A.publish(OUTPUT / "COMPLETE.json", A.canonical(dict(protocol=PROTOCOL, outcome=outcome, files=members)))
    print(json.dumps(dict(outcome=outcome, failure=failure, processes=len(results), speed_claimed=False), sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("receipt", "run"))
    parser.add_argument("path", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    args = parser.parse_args()
    if args.command == "receipt":
        A.require(args.output is not None, "receipt output required")
        A.publish(args.output, A.canonical(receipt(args.path)))
    else:
        A.require(args.output is None, "run takes one receipt")
        run(args.path)
