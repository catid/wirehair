#!/usr/bin/env python3
"""One-shot, rank-only .61 launcher. No codec or performance claims.

The scientific contract is frozen in wirehair-sxvz.16.1.20.61. Receipt creation
does not import or execute the mathematical worker. Historical experiments are
never opened for writing. Only this controller's freshly claimed directory is
published, and an existing directory is never reused.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import selectors
import stat
import subprocess
import sys
import time


PROTOCOL = "wirehair.wh2.noncommuting-radix-r0"
OUTPUT = Path("/var/tmp/wh2-noncommuting-radix-r0")
BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
SOURCES = (
    "bench/Wh2NoncommutingRadixR0.py",
    "bench/Wh2NoncommutingRadixRunR0.py",
    "bench/test_Wh2NoncommutingRadixR0.py",
    "bench/test_Wh2NoncommutingRadixRunR0.py",
)
MAX_RECEIPTS = 8 * 1024 * 1024
MAX_STDOUT = 4 * 1024 * 1024
MAX_STDERR = 1024 * 1024
ENV = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C", "TZ": "UTC"}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("ascii")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def identity(info):
    return (info.st_dev, info.st_ino, info.st_mode, info.st_nlink,
            info.st_uid, info.st_gid, info.st_size, info.st_mtime_ns,
            info.st_ctime_ns)


def time_left(deadline, maximum=5):
    if deadline is None:
        return maximum
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("shared controller deadline")
    return min(maximum, remaining)


def read_regular(path, cap, single_link=True, deadline=None):
    """Read stable bounded regular bytes; installed interpreters may hardlink."""
    path = Path(path)
    time_left(deadline)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode) or before.st_size > cap:
        raise ValueError("not a bounded regular file: " + str(path))
    if before.st_nlink < 1 or (single_link and before.st_nlink != 1):
        raise ValueError("unexpected link count: " + str(path))
    fd = os.open(str(path), os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        if identity(os.fstat(fd)) != identity(before):
            raise ValueError("file identity changed before read")
        chunks = []
        remaining = cap + 1
        while remaining:
            time_left(deadline)
            chunk = os.read(fd, min(65536, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        data = b"".join(chunks)
        if len(data) != before.st_size or len(data) > cap:
            raise ValueError("file size changed or exceeded bound")
        if identity(os.fstat(fd)) != identity(before):
            raise ValueError("file identity changed during read")
        if identity(path.lstat()) != identity(before):
            raise ValueError("named file identity changed")
        return data
    finally:
        os.close(fd)


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def invalid_constant(value):
        raise ValueError("non-finite JSON number: " + value)

    return json.loads(data.decode("ascii"), object_pairs_hook=pairs,
                      parse_constant=invalid_constant)


def git_bytes(*args, deadline=None, cap=2 * 1024 * 1024):
    result = subprocess.run(["/usr/bin/git"] + list(args), cwd=str(ROOT),
                            env=ENV, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            timeout=time_left(deadline), check=True)
    if result.stderr or len(result.stdout) > cap:
        raise ValueError("unexpected git diagnostic/output")
    return result.stdout


def git(*args, deadline=None):
    return git_bytes(*args, deadline=deadline, cap=65536).decode("ascii").strip()


def current_receipt(deadline=None):
    if sys.version_info[:2] not in ((3, 8), (3, 12)):
        raise ValueError("unvalidated interpreter version")
    if git("status", "--porcelain", "--", *SOURCES, deadline=deadline):
        raise ValueError("screen sources must be committed and clean")
    tracked = git("ls-files", "--", *SOURCES, deadline=deadline).splitlines()
    if sorted(tracked) != sorted(SOURCES):
        raise ValueError("screen sources must all be tracked")
    commit = git("rev-parse", "HEAD", deadline=deadline)
    if commit != git("rev-parse", "@{upstream}", deadline=deadline):
        raise ValueError("screen commit must already be pushed")
    source_hashes = {}
    for name in SOURCES:
        source = read_regular(ROOT / name, 2 * 1024 * 1024, deadline=deadline)
        committed = git_bytes("show", commit + ":" + name, deadline=deadline)
        if source != committed:
            raise ValueError("source differs from its exact committed blob: " + name)
        source_hashes[name] = sha(source)
    executable = str(Path(sys.executable).resolve(strict=True))
    return {
        "protocol": PROTOCOL,
        "source_commit": commit,
        "sources": source_hashes,
        "interpreter": executable,
        "interpreter_sha256": sha(read_regular(executable, 32 * 1024 * 1024,
                                               single_link=False, deadline=deadline)),
        "python_version": list(sys.version_info[:3]),
        "worker_argv": [executable, "-I", "-B", "-S",
                        str(ROOT / SOURCES[0]), "--worker"],
        "environment": ENV,
        "worker_seconds": 60,
        "outer_seconds": 70,
        "max_receipt_bytes": MAX_RECEIPTS,
        "max_stdout_bytes": MAX_STDOUT,
        "max_stderr_bytes": MAX_STDERR,
    }


def write_new(path, data):
    fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                 os.O_NOFOLLOW, 0o600)
    try:
        with os.fdopen(fd, "wb", closefd=False) as stream:
            stream.write(data)
            stream.flush()
            os.fsync(fd)
        os.fchmod(fd, 0o400)
    finally:
        os.close(fd)


def capture_worker(argv, deadline, stdout_cap=MAX_STDOUT,
                   stderr_cap=MAX_STDERR):
    """Only one owned child; live output bounds and an absolute deadline."""
    started = time.monotonic()
    streams = {"stdout": bytearray(), "stderr": bytearray()}
    failure = None
    selector = selectors.DefaultSelector()
    child = None
    try:
        child = subprocess.Popen(argv, cwd=str(ROOT), env=ENV,
                                 stdin=subprocess.DEVNULL,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for stream, name, cap in ((child.stdout, "stdout", stdout_cap),
                                  (child.stderr, "stderr", stderr_cap)):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, (name, cap))
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = "worker deadline"
                break
            for key, _ in selector.select(min(remaining, 0.1)):
                name, cap = key.data
                try:
                    chunk = os.read(key.fileobj.fileno(), 65536)
                except BlockingIOError:
                    continue
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                room = cap - len(streams[name])
                streams[name].extend(chunk[:room])
                if len(chunk) > room:
                    failure = name + " cap"
                    break
            if failure:
                break
        if failure:
            child.kill()
        try:
            child.wait(timeout=max(0.001, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            failure = failure or "worker deadline"
            child.kill()
            child.wait(timeout=2)
    finally:
        try:
            selector.close()
        finally:
            if child is not None:
                try:
                    if child.poll() is None:
                        child.kill()
                        child.wait(timeout=2)
                finally:
                    try:
                        child.stdout.close()
                    finally:
                        child.stderr.close()
    return {"stdout": bytes(streams["stdout"]),
            "stderr": bytes(streams["stderr"]),
            "returncode": child.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "failure": failure}


def run(receipt_path, started):
    # Source/interpreter checks and namespace claim precede the only worker.
    receipt_bytes = read_regular(receipt_path, 1024 * 1024, deadline=started + 8)
    receipt = strict_json(receipt_bytes)
    if receipt_bytes != canonical(receipt) or receipt != current_receipt(deadline=started + 8):
        raise ValueError("receipt does not match exact pushed inputs")
    OUTPUT.mkdir(mode=0o700, parents=False, exist_ok=False)
    members = {}
    total = 0

    def publish(name, data):
        nonlocal total
        if total + len(data) > MAX_RECEIPTS:
            raise ValueError("aggregate receipt cap")
        write_new(OUTPUT / name, data)
        total += len(data)
        members[name] = {"sha256": sha(data), "bytes": len(data)}

    claim = {"protocol": PROTOCOL, "receipt_sha256": sha(receipt_bytes),
             "receipt": receipt}
    publish("CLAIM.json", canonical(claim))
    stdout, stderr = b"", b""
    result = None
    observation = {}
    failure = None
    outcome = "INVALID"
    try:
        worker_start = time.monotonic()
        if worker_start - started >= 10:
            raise ValueError("insufficient time for the frozen worker budget")
        observation = capture_worker(receipt["worker_argv"],
                                     min(worker_start + 60, started + 68))
        stdout = observation.pop("stdout")
        stderr = observation.pop("stderr")
        if observation["failure"] or observation["returncode"] != 0 or stderr:
            raise ValueError("worker failed, timed out, or emitted diagnostics")
        if observation["elapsed_seconds"] > 60:
            raise ValueError("worker exceeded scientific wall-time bound")
        result = strict_json(stdout)
        if stdout != canonical(result) + b"\n":
            raise ValueError("worker output is not canonical JSON plus newline")
        if not isinstance(result, dict) or result.get("protocol") != PROTOCOL:
            raise ValueError("wrong worker result protocol/type")
        if result.get("outcome") not in ("PASS", "FAIL", "EXHAUSTED"):
            raise ValueError("worker has no complete scientific verdict")
        outcome = result["outcome"]
    except Exception as error:
        failure = type(error).__name__ + ": " + str(error)[:1000]
    # Failures also receive an explicit bounded post-pin observation.
    try:
        if current_receipt(deadline=started + 68) != receipt:
            raise ValueError("source/interpreter identity changed after worker")
        observation["post_pins"] = "match"
    except Exception as error:
        observation["post_pins"] = "failed_or_unobserved"
        failure = (failure + "; " if failure else "") + type(error).__name__ + ": " + str(error)[:1000]
        outcome = "INVALID"
    publish("raw.json", stdout)
    publish("stderr.txt", stderr)
    if time.monotonic() - started >= 69:
        failure = "outer publication deadline"
        outcome = "INVALID"
    summary = {"protocol": PROTOCOL, "outcome": outcome,
               "failure": failure, "worker": observation,
               "elapsed_before_publication": time.monotonic() - started,
               "raw_sha256": sha(stdout), "WH1_compared": False,
               "promotion_claimed": False, "all_K_claimed": False,
               "global_injectivity_claimed": False}
    publish("summary.json", canonical(summary))
    manifest = canonical({"protocol": PROTOCOL, "outcome": outcome,
                          "files": members})
    if total + len(manifest) > MAX_RECEIPTS:
        raise ValueError("terminal manifest exceeds aggregate receipt cap")
    write_new(OUTPUT / "COMPLETE.json", manifest)
    directory_fd = os.open(str(OUTPUT), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    if time.monotonic() - started >= 70:
        # Immutable evidence remains; the independent CLI observer must reject.
        raise ValueError("whole CLI exceeded scientific deadline")
    print(canonical(summary).decode("ascii"))
    return 0 if outcome != "INVALID" else 1


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--make-receipt", metavar="PATH")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--receipt", metavar="PATH")
    args = parser.parse_args()
    if args.make_receipt:
        if args.receipt:
            parser.error("--receipt is only valid with --run")
        write_new(Path(args.make_receipt), canonical(current_receipt()))
        return 0
    if not args.receipt:
        parser.error("--run requires --receipt")
    return run(Path(args.receipt), started)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
