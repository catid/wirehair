#!/usr/bin/env python3
"""Fixed-map, fresh-root public MIX3 portability falsifier; never certification."""
import argparse
import contextlib
import ctypes as ct
import importlib.util
import math
from pathlib import Path
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
WORKER_SOURCE = Path(__file__).resolve()
PROTOCOL = "wirehair.wh2.production-mix3-k3k6-portability-r0"
WIDTHS, SELECTED = (2, 64, 1280), (5, 3)
MAP_PATH = Path("/var/tmp/wh2-production-mix3-k3k6-broad-r0/selection.json")
MAP_SHA256 = "41d0cfbdb0a22ccbc88d615ea028154798449d4ff1e19ca2eb1811fca9f43543"
ROSTER_SHA256 = "a739277bae151083192f2cf9c9c56fe5c321467a639775b4821ba4f03293c53d"
SOURCE_PATHS = tuple(sorted((
    "bench/Wh2ProductionMix3RecoveryPortabilityR0.py",
    "bench/test_Wh2ProductionMix3RecoveryPortabilityR0.py",
    "bench/Wh2ProductionMix3RecoveryBroadR0.py", "bench/Wh2ProductionMix3RecoveryR0.py",
    "include/wirehair/wirehair.h")))
spec = importlib.util.spec_from_file_location("_wh2_portability_broad_private",
                                             ROOT / "bench/Wh2ProductionMix3RecoveryBroadR0.py")
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)
c = b.c
HOLDOUT_ROOTS = tuple("0x" + c.digest((PROTOCOL + ":holdout/" + str(i)).encode("ascii"))[:16]
                     for i in range(32))
c.require(c.digest(c.canonical(HOLDOUT_ROOTS)) == ROSTER_SHA256, "root roster drift")
c.PROTOCOL, c.HOLDOUT_ROOTS = PROTOCOL, HOLDOUT_ROOTS
_trace_b2 = c.trace


def source_bytes(K, B):
    return bytes((73 * i + 19 * K + 11) & 255 for i in range(K * B))


def profile_bytes(K, B, attempt):
    c.require(K in c.KS and B in WIDTHS and type(attempt) is int and 0 <= attempt <= 255,
              "profile dimensions/attempt")
    return b"WHV2" + struct.pack("<HHQQIB3s", 1, 32, 0x4b295bbb47f4f9c9,
                                K * B, B, attempt, bytes(3))


def trace(K, B, root, schedule):
    # FrozenTraceSeed = root XOR K*increment XOR B*multiplier. Cancel the
    # reviewed helper's B2 term and substitute actual B, modulo uint64.
    adjusted = (int(root, 16) ^ (2 * 0xbf58476d1ce4e5b9) ^
                (B * 0xbf58476d1ce4e5b9)) & c.MASK64
    return _trace_b2(K, "0x{:016x}".format(adjusted), schedule)


def configure_width(B):
    c.require(type(B) is int and B in WIDTHS, "width")
    c.source_bytes = lambda K: source_bytes(K, B)
    c.profile_bytes = lambda K, attempt: profile_bytes(K, B, attempt)
    c.trace = lambda K, root, schedule: trace(K, B, root, schedule)
    c.validate_cell = lambda *args: validate_cell(B, *args)


def validate_cell(B, row, phase, arm, K, attempt, root_index, root, schedule):
    c.require(type(row) is dict and set(row) == c.CELL_KEYS, "cell schema")
    for key, expected in (("type", "cell"), ("phase", phase), ("arm", arm), ("K", K),
                          ("attempt", attempt), ("root_index", root_index),
                          ("root", root), ("schedule", schedule)):
        c.require(type(row[key]) is type(expected) and row[key] == expected, "cell " + key)
    profile, source = profile_bytes(K, B, attempt), source_bytes(K, B)
    source_hash = c.digest(source)
    c.require(row["profile_hex"] == profile.hex() and row["profile_sha256"] == c.digest(profile),
              "profile bytes/hash")
    c.require(row["source_sha256"] == source_hash, "source hash")
    ids, attempted = trace(K, B, root, schedule)
    c.require(type(row["ids"]) is list and all(type(x) is int for x in row["ids"]) and
              row["ids"] == ids, "raw packet IDs")
    c.require(type(row["attempted_candidates"]) is int and row["attempted_candidates"] == attempted,
              "attempted candidates")
    c.require(row["trace_sha256"] == c.digest(struct.pack("<" + "I" * len(ids), *ids)), "trace hash")
    outcomes = row["outcomes"]
    c.require(type(outcomes) is list and len(outcomes) == 4 and
              all(type(x) is int and x in (0, 1, 10) for x in outcomes), "outcome schema")
    if 10 in outcomes:
        c.require(arm == "candidate" and outcomes == [10] * 4 and row["packets_hex"] is None,
                  "construction failure")
    else:
        packets = row["packets_hex"]
        c.require(type(packets) is str and len(packets) == (K + 4) * B * 2 and
                  all(x in "0123456789abcdef" for x in packets), "packet payload receipt")
        c.require(all(outcomes[i] >= outcomes[i + 1] for i in range(3)), "non-nested outcomes")
        payload = bytes.fromhex(packets)
        for i, packet_id in enumerate(ids):
            if packet_id < K:
                c.require(payload[B * i:B * (i + 1)] == source[B * packet_id:B * (packet_id + 1)],
                          "systematic packet receipt")
    c.require(row["recovered_sha256"] == [source_hash if x == 0 else None for x in outcomes],
              "recovered byte oracle")


class Profile(ct.Structure):
    _fields_ = [("struct_bytes", ct.c_uint32), ("profile_version", ct.c_uint32),
                ("profile_id", ct.c_uint64), ("message_bytes", ct.c_uint64),
                ("block_bytes", ct.c_uint32), ("seed_attempt", ct.c_uint8),
                ("reserved", ct.c_uint8 * 3)]


class Native:
    def __init__(self):
        c.require(c.file_digest(b.PINNED_LIBRARY) == b.LIBRARY_SHA256, "pinned library drift")
        self.lib = ct.CDLL(str(b.PINNED_LIBRARY))
        u32, u64, ptr = ct.c_uint32, ct.c_uint64, ct.c_void_p
        signatures = {
            "wirehair_init_": ([ct.c_int], ct.c_int),
            "wirehair_v2_encoder_create": ([ptr, u64, u32, ptr, u32, ct.POINTER(u32), ct.POINTER(ptr)], ct.c_int),
            "wirehair_v2_encoder_create_profile": ([ptr, ptr, u32, ct.POINTER(ptr)], ct.c_int),
            "wirehair_v2_profile_deserialize": ([ptr, u32, ct.POINTER(Profile)], ct.c_int),
            "wirehair_v2_profile_serialize": ([ct.POINTER(Profile), ptr, u32, ct.POINTER(u32)], ct.c_int),
            "wirehair_v2_decoder_create": ([ptr, u32, ct.POINTER(ptr)], ct.c_int),
            "wirehair_v2_encode": ([ptr, u32, ptr, u32, ct.POINTER(u32)], ct.c_int),
            "wirehair_v2_decode": ([ptr, u32, ptr, u32], ct.c_int),
            "wirehair_v2_recover": ([ptr, ptr, u64, ct.POINTER(u64)], ct.c_int),
            "wirehair_v2_free": ([ptr], None)}
        for name, (arguments, result) in signatures.items():
            function = getattr(self.lib, name)
            function.argtypes, function.restype = arguments, result
        c.require(ct.sizeof(Profile) == 32 and Profile.seed_attempt.offset == 28, "profile ABI")
        c.require(self.lib.wirehair_init_(2) == 0, "initialization")

    @contextlib.contextmanager
    def encoder(self, K, B, attempt):
        source = source_bytes(K, B)
        message, profile, handle = ct.create_string_buffer(source, len(source)), ct.create_string_buffer(32), ct.c_void_p()
        try:
            written = ct.c_uint32()
            if attempt is None:
                result = self.lib.wirehair_v2_encoder_create(message, len(source), B, profile, 32,
                                                            ct.byref(written), ct.byref(handle))
                c.require(result == 0 and written.value == 32, "baseline create")
            else:
                profile.raw = profile_bytes(K, B, attempt)
                result = self.lib.wirehair_v2_encoder_create_profile(message, profile, 32, ct.byref(handle))
                c.require(result in (0, 10), "exact constructor result")
            c.require(bool(handle.value) == (result == 0) and message.raw == source, "constructor handle/source")
            decoded = Profile()
            c.require(self.lib.wirehair_v2_profile_deserialize(profile, 32, ct.byref(decoded)) == 0,
                      "deserialize profile")
            actual = decoded.seed_attempt
            c.require(profile.raw == profile_bytes(K, B, actual) and (attempt is None or actual == attempt),
                      "realized profile")
            serialized = ct.create_string_buffer(32)
            c.require(self.lib.wirehair_v2_profile_serialize(ct.byref(decoded), serialized, 32,
                                                           ct.byref(written)) == 0 and written.value == 32 and
                      serialized.raw == profile.raw, "canonical profile roundtrip")
            if result == 0:
                for packet_id in range(K):
                    c.require(self.encode(handle, packet_id, B) == source[B * packet_id:B * (packet_id + 1)],
                              "systematic bytes")
            yield handle, profile.raw, message, source, result == 10, actual
        finally:
            self.lib.wirehair_v2_free(handle)

    def encode(self, handle, packet_id, B):
        packet, written = ct.create_string_buffer(B), ct.c_uint32()
        c.require(self.lib.wirehair_v2_encode(handle, packet_id, packet, B, ct.byref(written)) == 0 and
                  written.value == B, "packet encode")
        return packet.raw

    def decode(self, profile_bytes_value, source, B, ids, packets, count):
        profile = ct.create_string_buffer(profile_bytes_value, 32)
        data, handle = ct.create_string_buffer(packets, len(packets)), ct.c_void_p()
        try:
            c.require(self.lib.wirehair_v2_decoder_create(profile, 32, ct.byref(handle)) == 0 and handle.value,
                      "decoder create")
            result = 1
            for i in range(count):
                result = self.lib.wirehair_v2_decode(handle, ids[i], ct.byref(data, B * i), B)
                c.require(result in (0, 1), "unexpected decode result")
                if result == 0:
                    break
            if result == 0:
                recovered, written = ct.create_string_buffer(len(source)), ct.c_uint64()
                c.require(self.lib.wirehair_v2_recover(handle, recovered, len(source), ct.byref(written)) == 0 and
                          written.value == len(source) and recovered.raw == source, "recovered bytes")
            return result
        finally:
            self.lib.wirehair_v2_free(handle)


def emit(value):
    print(c.canonical(value).decode("ascii"), flush=True)


def emit_arm(native, K, B, arm, attempt, roots, deadline):
    with native.encoder(K, B, attempt) as (handle, profile, message, source, bad, actual):
        for root_index, root in enumerate(roots):
            for schedule in c.SCHEDULES:
                c.require(time.monotonic() < deadline, "worker deadline")
                ids, attempted = trace(K, B, root, schedule)
                packets = None if bad else b"".join(native.encode(handle, packet_id, B) for packet_id in ids)
                outcomes = []
                for overhead in c.OVERHEADS:
                    c.require(time.monotonic() < deadline, "worker deadline")
                    outcomes.append(10 if bad else native.decode(profile, source, B, ids, packets, K + overhead))
                c.require(message.raw == source, "source mutated during packet/decode work")
                source_hash = c.digest(source)
                emit({"type": "cell", "phase": "holdout", "arm": arm, "K": K, "attempt": actual,
                      "root_index": root_index, "root": root, "schedule": schedule,
                      "profile_hex": profile.hex(), "profile_sha256": c.digest(profile), "source_sha256": source_hash,
                      "ids": ids, "attempted_candidates": attempted,
                      "trace_sha256": c.digest(struct.pack("<" + "I" * len(ids), *ids)),
                      "packets_hex": None if bad else packets.hex(), "outcomes": outcomes,
                      "recovered_sha256": [source_hash if result == 0 else None for result in outcomes]})
        return actual


def worker(B, seconds, commit):
    configure_width(B)
    c.require(0 < seconds <= 60 and len(commit) == 40 and all(x in "0123456789abcdef" for x in commit),
              "worker arguments")
    deadline = time.monotonic() + seconds
    native = Native()
    emit({"type": "begin", "protocol": PROTOCOL, "phase": "holdout", "source_commit": commit,
          "library_path": str(b.PINNED_LIBRARY)})
    for K, selected in zip(c.KS, SELECTED):
        emit_arm(native, K, B, "baseline", None, HOLDOUT_ROOTS, deadline)
        emit_arm(native, K, B, "candidate", selected, HOLDOUT_ROOTS, deadline)
    c.require(time.monotonic() < deadline, "worker deadline")
    emit({"type": "terminal", "phase": "holdout", "rows": 384, "selected_attempts": list(SELECTED)})


def selftest():
    native = Native()
    for B in WIDTHS:
        for K in c.KS:
            # Only the neutral loss root and actual baseline attempt are decoded;
            # never probe fixed candidate5/3 at new widths before the pilot.
            roots = ("0x1234567890abcdef",)
            actual = emit_arm(native, K, B, "baseline", None, roots, math.inf)
            emit_arm(native, K, B, "candidate", actual, roots, math.inf)
    emit({"type": "selftest", "pass": True})


def source_receipt():
    return {path: c.file_digest(ROOT / path) for path in SOURCE_PATHS}


def validate_receipt(receipt):
    b.validate_basis(receipt["production_build"])
    c.require(receipt["protocol"] == PROTOCOL and receipt["production_source_commit"] == b.PRODUCTION_COMMIT and
              receipt["library"] == str(b.PINNED_LIBRARY) and receipt["library_sha256"] == b.LIBRARY_SHA256,
              "pinned production identity")
    c.require(c.digest(c.canonical(receipt["map_basis"]) + b"\n") == MAP_SHA256 and
              receipt["map_basis"]["selected_attempts"] == list(SELECTED), "fixed map basis")


def check_build(receipt):
    validate_receipt(receipt)
    c.require(c.command(["git", "rev-parse", "HEAD"]).strip() == receipt["source_commit"], "build commit drift")
    c.require(source_receipt() == receipt["source_files"], "build source drift")
    for key in ("worker", "library", "interpreter"):
        c.require(c.file_digest(receipt[key]) == receipt[key + "_sha256"], key + " changed")
    for path in set(SOURCE_PATHS) & set(receipt["production_build"]["source_files"]):
        c.require(receipt["source_files"][path] == receipt["production_build"]["source_files"][path],
                  "inherited source drift")
    c.require(c.file_digest(MAP_PATH) == MAP_SHA256, "fixed map file changed")


def build(output):
    output = Path(output).resolve()
    c.require(not output.exists(), "build output already exists")
    basis = b.pinned_library_receipt()
    c.require(c.file_digest(MAP_PATH) == MAP_SHA256, "fixed map file changed")
    interpreter = str(Path(sys.executable).resolve())
    receipt = {"protocol": PROTOCOL, "source_commit": c.command(["git", "rev-parse", "HEAD"]).strip(),
               "source_files": source_receipt(), "worker": str(WORKER_SOURCE), "worker_sha256": c.file_digest(WORKER_SOURCE),
               "interpreter": interpreter, "interpreter_sha256": c.file_digest(interpreter),
               "interpreter_version": sys.version, "interpreter_flags": ["-I", "-B", "-S"],
               "library": str(b.PINNED_LIBRARY), "library_sha256": b.LIBRARY_SHA256,
               "production_source_commit": b.PRODUCTION_COMMIT, "production_build": basis,
               "map_basis": c.parse_json(MAP_PATH.read_bytes())}
    check_build(receipt)
    output.mkdir()
    c.write_json(output / "build.json", receipt)
    return receipt


def frozen_protocol(receipt):
    validate_receipt(receipt)
    return {"protocol": PROTOCOL, "build": receipt, "K": list(c.KS), "widths": list(WIDTHS),
            "selected_attempts": list(SELECTED), "holdout_roots": list(HOLDOUT_ROOTS),
            "roster_sha256": ROSTER_SHA256, "schedules": list(c.SCHEDULES), "loss_ppm": 500000,
            "overheads": list(c.OVERHEADS), "training": False, "worker_budget_seconds": 60,
            "outer_budget_seconds": 70, "max_prefix_decodes": 4608}


def summarize(parsed):
    c.require(len(parsed) == len(WIDTHS), "width roster")
    widths = [{"block_bytes": B, "baseline_attempts": phase["baseline_attempts"],
               "paired": c.paired_counts(phase, SELECTED)} for B, phase in zip(WIDTHS, parsed)]
    return {"outcome": "FAIL" if any(w["paired"]["candidate_failures"][0] for w in widths) else "PASS",
            "selected_attempts": list(SELECTED), "widths": widths,
            "known_prior_counterexample_resolved": False}


def parse_width(bundle, B, receipt):
    configure_width(B)
    return c.parse_phase((bundle / ("b%d.jsonl" % B)).read_text(), "holdout",
                         receipt["source_commit"], receipt["library"], SELECTED)


def run_once(build_path, bundle_path):
    receipt = c.parse_json(Path(build_path).read_bytes())
    check_build(receipt)
    c.require(not c.command(["git", "status", "--porcelain", "--untracked-files=no", "--"] +
                            list(SOURCE_PATHS)).strip(), "harness tree dirty")
    for path in SOURCE_PATHS:
        c.command(["git", "ls-files", "--error-unmatch", "--", path])
    c.require(c.command(["git", "rev-parse", "@{upstream}"]).strip() == receipt["source_commit"], "source not pushed")
    bundle = Path(bundle_path).resolve()
    bundle.mkdir(exist_ok=False)
    c.write_json(bundle / "freeze.json", frozen_protocol(receipt))
    c.write_json(bundle / "selection.json", receipt["map_basis"])
    started = time.monotonic()
    try:
        parsed = []
        for B in WIDTHS:
            remaining = 60 - (time.monotonic() - started)
            c.require(remaining >= 1, "aggregate worker budget exhausted")
            argv = [receipt["interpreter"], "-I", "-B", "-S", receipt["worker"], "--worker", str(B),
                    "--seconds", str(int(remaining)), "--source-commit", receipt["source_commit"]]
            c.execute_worker(argv, bundle / ("b%d.jsonl" % B), bundle / ("b%d.stderr" % B), remaining)
            parsed.append(parse_width(bundle, B, receipt))
        c.require(c.file_digest(bundle / "selection.json") == MAP_SHA256, "sealed map changed")
        result = summarize(parsed)
        check_build(receipt)
        c.require(time.monotonic() - started < 70, "outer deadline")
    except Exception as error:
        result = {"outcome": "INVALID", "error": str(error)}
    result.update({"protocol": PROTOCOL, "elapsed_seconds": time.monotonic() - started,
                   "promotion_claimed": False, "all_K_claimed": False, "timing_claimed": False})
    return c.publish_result(bundle, result, started)


def replay(bundle_path):
    bundle = Path(bundle_path)
    complete = c.parse_json((bundle / "COMPLETE").read_bytes())
    c.require(type(complete) is dict and set(complete) == {p.name for p in bundle.iterdir() if p.name != "COMPLETE"},
              "bundle files")
    for name, expected in complete.items():
        c.require(Path(name).name == name and c.file_digest(bundle / name) == expected, "bundle digest")
    freeze = c.parse_json((bundle / "freeze.json").read_bytes())
    c.require(c.canonical(freeze) == c.canonical(frozen_protocol(freeze["build"])), "frozen protocol")
    result = c.parse_json((bundle / "summary.json").read_bytes())
    c.require(result["protocol"] == PROTOCOL and all(result[key] is False for key in
              ("promotion_claimed", "all_K_claimed", "timing_claimed")), "summary scope")
    if result["outcome"] == "INVALID":
        return result
    c.require(result["outcome"] in ("PASS", "FAIL"), "outcome")
    c.require(type(result["elapsed_seconds"]) in (int, float) and 0 <= result["elapsed_seconds"] < 70,
              "elapsed bound")
    expected_files = {"freeze.json", "selection.json", "summary.json"}
    for B in WIDTHS:
        expected_files.update(("b%d.jsonl" % B, "b%d.stderr" % B))
        c.require((bundle / ("b%d.stderr" % B)).read_bytes() == b"", "worker stderr")
    c.require(set(complete) == expected_files and c.file_digest(bundle / "selection.json") == MAP_SHA256,
              "scientific roster/map")
    expected = summarize([parse_width(bundle, B, freeze["build"]) for B in WIDTHS])
    c.require(set(result) == set(expected) | {"protocol", "elapsed_seconds", "promotion_claimed", "all_K_claimed",
                                             "timing_claimed"} and
              all(c.canonical(result[key]) == c.canonical(value) for key, value in expected.items()), "summary mismatch")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build")
    modes.add_argument("--run")
    modes.add_argument("--replay")
    modes.add_argument("--worker", type=int, choices=WIDTHS)
    modes.add_argument("--selftest", action="store_true")
    parser.add_argument("--build-receipt")
    parser.add_argument("--seconds", type=int)
    parser.add_argument("--source-commit")
    args = parser.parse_args()
    if args.worker:
        c.require(args.seconds is not None and args.source_commit is not None, "worker needs deadline/source")
        worker(args.worker, args.seconds, args.source_commit)
    elif args.selftest:
        selftest()
    elif args.build:
        emit(build(args.build))
    elif args.run:
        c.require(args.build_receipt is not None, "run needs build-receipt")
        emit(run_once(args.build_receipt, args.run))
    else:
        emit(replay(args.replay))


if __name__ == "__main__":
    try:
        main()
    except (c.ValidationError, OSError, subprocess.SubprocessError) as error:
        print("production MIX3 portability r0: " + str(error), file=sys.stderr)
        sys.exit(2)
