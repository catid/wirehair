#!/usr/bin/env python3
"""Fresh-root finite seed screen using the unchanged, receipt-pinned MIX3 DSO."""
import argparse
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = "wirehair.wh2.production-mix3-k3k6-broad-r0"
WORKER_SOURCE = ROOT / "bench/Wh2ProductionMix3RecoveryBroadR0.cpp"
PINNED_LIBRARY = Path("/tmp/wh2-page-aa.uWaCGu/libwirehair.so.2.0.0")
LIBRARY_SHA256 = "c0bf0666e3cc51523e18847b8a07384f2eac312518d0b4f9ac48cd14e63fa038"
BASIS_FREEZE = Path("/var/tmp/wh2-production-mix3-k3k6-r0/freeze.json")
BASIS_FREEZE_SHA256 = "559050987e6b683c85f8b81dc69d33c1d4e2b1e54f88b364a703ae8c17384de0"
BASIS_BUILD_SHA256 = "beae9fbb37d15e5412810840c13eb1a29ac12d792627a1170eca1b94f16267f9"
PRODUCTION_COMMIT = "c848277ab5d2202b71edf576838b61697a7951eb"
ROSTER_SHA256 = "2a8ae1333952f4cbd28e4af15cdd689e2314506b18526c377f0a5bb6af79e202"
SOURCE_PATHS = tuple(sorted((
    "bench/Wh2ProductionMix3RecoveryBroadR0.cpp", "bench/Wh2ProductionMix3RecoveryBroadR0.py",
    "bench/test_Wh2ProductionMix3RecoveryBroadR0.py", "bench/Wh2ProductionMix3RecoveryR0.cpp",
    "bench/Wh2ProductionMix3RecoveryR0.py", "bench/test_Wh2ProductionMix3RecoveryR0.py",
    "bench/Wh2FrozenTrace.cpp", "bench/Wh2FrozenTrace.h", "include/wirehair/wirehair.h")))

# Private execution namespace: the frozen r0 module and its replay are untouched.
spec = importlib.util.spec_from_file_location("_wh2_mix3_broad_r0_private",
                                             ROOT / "bench/Wh2ProductionMix3RecoveryR0.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
TRAINING_ROOTS = tuple("0x" + hashlib.sha256((PROTOCOL + ":train/" + str(i)).encode("ascii"))
                      .hexdigest()[:16] for i in range(16))
HOLDOUT_ROOTS = tuple("0x" + hashlib.sha256((PROTOCOL + ":holdout/" + str(i)).encode("ascii"))
                     .hexdigest()[:16] for i in range(64))
c.require(c.digest(c.canonical({"holdout": HOLDOUT_ROOTS, "train": TRAINING_ROOTS})) == ROSTER_SHA256,
          "root roster drift")
c.PROTOCOL, c.TRAINING_ROOTS, c.HOLDOUT_ROOTS = PROTOCOL, TRAINING_ROOTS, HOLDOUT_ROOTS
c.WORKER_SOURCE = WORKER_SOURCE
_original_command, _original_check_build = c.command, c.check_build


def source_receipt():
    # These are the actual new harness inputs, not today's unrelated production tree.
    return {path: c.file_digest(ROOT / path) for path in SOURCE_PATHS}


def validate_basis(basis):
    c.require(c.digest(c.canonical(basis)) == BASIS_BUILD_SHA256, "production basis drift")
    c.require(basis["source_commit"] == PRODUCTION_COMMIT and
              basis["library"] == str(PINNED_LIBRARY) and basis["library_sha256"] == LIBRARY_SHA256,
              "production basis identity")


def pinned_library_receipt():
    c.require(c.file_digest(BASIS_FREEZE) == BASIS_FREEZE_SHA256, "baseline freeze drift")
    basis = c.parse_json(BASIS_FREEZE.read_bytes())["build"]
    validate_basis(basis)
    c.require(c.file_digest(PINNED_LIBRARY) == LIBRARY_SHA256, "pinned library drift")
    return basis


def check_build(receipt):
    validate_basis(receipt["production_build"])
    c.require(receipt["production_source_commit"] == PRODUCTION_COMMIT and
              receipt["library"] == str(PINNED_LIBRARY) and receipt["library_sha256"] == LIBRARY_SHA256,
              "pinned production identity")
    _original_check_build(receipt)
    # Included helpers and the public header must match the old DSO's source basis.
    for path in set(SOURCE_PATHS) & set(receipt["production_build"]["source_files"]):
        c.require(receipt["source_files"][path] == receipt["production_build"]["source_files"][path],
                  "inherited source drift: " + path)


def build(output, compiler="/usr/bin/g++-13"):
    output, compiler = Path(output).resolve(), Path(compiler).resolve(strict=True)
    c.require(not output.exists(), "build output already exists")
    basis, sources = pinned_library_receipt(), source_receipt()
    commit = c.command(["git", "rev-parse", "HEAD"]).strip()
    output.mkdir()
    worker = output / "worker"
    argv = [str(compiler), "-std=c++11", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(commit), "-I" + str(ROOT / "include"),
            str(WORKER_SOURCE), str(ROOT / "bench/Wh2FrozenTrace.cpp"), str(PINNED_LIBRARY),
            "-Wl,-rpath," + str(PINNED_LIBRARY.parent), "-ldl", "-o", str(worker)]
    output_text = c.command(argv, timeout=60)
    c.require(source_receipt() == sources, "sources changed during build")
    description = c.parse_json(c.command([str(worker), "--describe"]))
    c.require(description == {"type": "begin", "protocol": PROTOCOL, "phase": "describe",
                              "source_commit": commit, "library_path": str(PINNED_LIBRARY)}, "worker binding")
    receipt = {"protocol": PROTOCOL, "source_commit": commit, "source_files": sources,
               "compiler": str(compiler), "compiler_sha256": c.file_digest(compiler),
               "compiler_version": c.command([str(compiler), "--version"]), "command": argv,
               "build_output": output_text, "library": str(PINNED_LIBRARY), "library_sha256": LIBRARY_SHA256,
               "worker": str(worker), "worker_sha256": c.file_digest(worker),
               "production_source_commit": PRODUCTION_COMMIT, "production_build": basis}
    check_build(receipt)
    c.write_json(output / "build.json", receipt)
    return receipt


def frozen_protocol(receipt):
    # Offline replay verifies the embedded original production receipt without
    # requiring its old build directory or inspecting today's production source.
    validate_basis(receipt["production_build"])
    c.require(receipt["production_source_commit"] == PRODUCTION_COMMIT and
              receipt["library"] == str(PINNED_LIBRARY) and receipt["library_sha256"] == LIBRARY_SHA256,
              "frozen production identity")
    return {"protocol": PROTOCOL, "build": receipt, "training_roots": list(TRAINING_ROOTS),
            "holdout_roots": list(HOLDOUT_ROOTS), "K": list(c.KS), "block_bytes": 2,
            "schedules": list(c.SCHEDULES), "loss_ppm": 500000, "overheads": list(c.OVERHEADS),
            "worker_budget_seconds": 60, "outer_budget_seconds": 70, "max_prefix_decodes": 101760,
            "roster_sha256": ROSTER_SHA256, "production_source_commit": PRODUCTION_COMMIT}


def harness_command(argv, timeout=30):
    # The inherited run gate is intentionally scoped to compiled/imported harness
    # inputs. Unrelated production edits cannot alter the separately pinned DSO.
    if argv == ["git", "status", "--porcelain", "--untracked-files=no"]:
        argv = argv + ["--"] + list(SOURCE_PATHS)
    return _original_command(argv, timeout)


c.source_receipt, c.check_build, c.frozen_protocol = source_receipt, check_build, frozen_protocol
c.command = harness_command
run_once, replay = c.run_once, c.replay


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build", metavar="NEW_DIRECTORY")
    modes.add_argument("--run", metavar="NEW_BUNDLE")
    modes.add_argument("--replay", metavar="BUNDLE")
    parser.add_argument("--compiler", default="/usr/bin/g++-13")
    parser.add_argument("--build-receipt")
    args = parser.parse_args()
    if args.build:
        result = build(args.build, args.compiler)
    elif args.run:
        c.require(args.build_receipt, "run needs build-receipt")
        result = run_once(args.build_receipt, args.run)
    else:
        result = replay(args.replay)
    print(c.canonical(result).decode("ascii"))


if __name__ == "__main__":
    try:
        main()
    except (c.ValidationError, OSError, subprocess.SubprocessError) as error:
        print("production MIX3 recovery broad r0: " + str(error), file=sys.stderr)
        sys.exit(2)
