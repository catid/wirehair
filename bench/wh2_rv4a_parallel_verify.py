#!/usr/bin/env python3
"""Source-pinned bounded verifier for a completed WH2 rv4a campaign.

The pin record is external, canonical JSON and contains every benchmark,
native, parser, runner, test and build/configuration digest bound by the
campaign manifest.  Its expected SHA-256 is a required command-line input.
This avoids a runner/verifier/pin-record self-hash cycle while still requiring
an out-of-band trusted digest before any campaign-provided Python is executed.

Only exact, retained-descriptor source snapshots are compiled.  The frozen
campaign runner performs the bounded multi-process replay and exact completion
reproduction; this wrapper independently re-hashes the complete runtime
binding before and after that replay.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import types


PIN_SCHEMA = "wirehair.wh2.rv4a.verifier-pins.v1"
CAMPAIGN_SCHEMA = "wirehair.wh2.rv4a.repair-campaign.v1"
REPAIRTIMING_PROTOCOL = "wirehair-v2-bench:repairtiming:repair-v1:v3"
REPAIRTIMING_SCHEMA = "wirehair.wh2.repairtiming.v3"
MAX_REPLAY_WORKERS = 32
PIN_BYTE_LIMIT = 2 * 1024 * 1024
MANIFEST_BYTE_LIMIT = 64 * 1024 * 1024

SOURCE_CLASSES = {
    "native_benchmark": "native",
    "native_repairtiming": "native",
    "repair_contract": "native",
    "repair_contract_header": "native",
    "native_cli_test": "native",
    "native_repairtiming_cli_test": "native",
    "build_policy_e2e": "build",
    "codec_build": "build",
    "root_build": "build",
    "parser": "python",
    "context_tool": "python",
    "runner": "python",
    "parser_test": "python",
    "runner_test": "python",
    "parallel_verifier": "python",
    "parallel_verifier_test": "python",
    "build_configuration": "configuration",
}
FILE_BYTE_LIMITS = {
    "binary": 512 * 1024 * 1024,
    "native": 8 * 1024 * 1024,
    "python": 2 * 1024 * 1024,
    "build": 4 * 1024 * 1024,
    "configuration": 16 * 1024 * 1024,
}
RUNTIME_FILE_NAMES = ("benchmark", *SOURCE_CLASSES)
_LOADED_CAMPAIGNS = {}
_EXECUTING_VERIFIER_SOURCE_SHA256 = globals().get(
    "_EXECUTING_VERIFIER_SOURCE_SHA256")


class PinError(RuntimeError):
    """A source pin, runtime binding, or frozen replay was invalid."""


def canonical_json_bytes(value):
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ) + "\n"
    ).encode("ascii")


def sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def _is_sha256(value):
    return (
        type(value) is str and len(value) == 64 and
        all(character in "0123456789abcdef" for character in value)
    )


def _reject_duplicate_pairs(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"duplicate JSON key {key!r}")
        value[key] = item
    return value


def _strict_json(payload, label):
    try:
        value = json.loads(
            payload.decode("ascii"),
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=lambda token: (
                (_ for _ in ()).throw(
                    ValueError(f"invalid JSON constant {token}"))
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise PinError(f"{label} is not strict JSON: {error}")
    if canonical_json_bytes(value) != payload:
        raise PinError(f"{label} is not canonical JSON")
    return value


def _stable_file(path, *, byte_limit, retain_payload):
    """Hash one regular file through one retained, no-follow descriptor."""
    if type(byte_limit) is not int or byte_limit < 1:
        raise PinError("stable-file byte limit is invalid")
    path = Path(path)
    if not path.is_absolute():
        raise PinError(f"runtime path is not absolute: {path}")
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = None
    try:
        descriptor = os.open(str(path), flags)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode) or
            before.st_size > byte_limit
        ):
            raise PinError(f"runtime file is invalid or oversized: {path}")
        digest = hashlib.sha256()
        chunks = [] if retain_payload else None
        size = 0
        while True:
            chunk = os.read(
                descriptor, min(1024 * 1024, byte_limit - size + 1))
            if not chunk:
                break
            size += len(chunk)
            if size > byte_limit:
                raise PinError(f"runtime file exceeds byte limit: {path}")
            digest.update(chunk)
            if chunks is not None:
                chunks.append(chunk)
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise PinError(f"could not read runtime file {path}: {error}")
    finally:
        if descriptor is not None:
            os.close(descriptor)
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        size != after.st_size
    ):
        raise PinError(f"runtime file changed while read: {path}")
    return {
        "path": str(path),
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": after.st_mode,
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
        "sha256": digest.hexdigest(),
        "payload": None if chunks is None else b"".join(chunks),
    }


def _read_canonical(path, *, byte_limit):
    snapshot = _stable_file(
        Path(path).absolute(),
        byte_limit=byte_limit,
        retain_payload=True,
    )
    payload = snapshot.pop("payload")
    return _strict_json(payload, str(path)), sha256_bytes(payload)


def _runtime_files(bindings):
    return {
        "benchmark": bindings["benchmark"],
        **bindings["sources"],
    }


def _file_class(name):
    return "binary" if name == "benchmark" else SOURCE_CLASSES[name]


def _validate_runtime_binding_layout(bindings):
    file_fields = {
        "path", "device", "inode", "mode", "size",
        "mtime_ns", "ctime_ns", "sha256", "byte_limit",
    }
    if (
        not isinstance(bindings, dict) or
        set(bindings) != {"benchmark", "sources", "thermal"} or
        not isinstance(bindings.get("sources"), dict) or
        set(bindings["sources"]) != set(SOURCE_CLASSES)
    ):
        raise PinError("runtime binding layout is invalid")
    for name, item in _runtime_files(bindings).items():
        limit = FILE_BYTE_LIMITS[_file_class(name)]
        if (
            not isinstance(item, dict) or set(item) != file_fields or
            type(item.get("path")) is not str or
            not Path(item["path"]).is_absolute() or
            any(type(item.get(field)) is not int or item[field] < 0
                for field in (
                    "device", "inode", "mode", "size",
                    "mtime_ns", "ctime_ns", "byte_limit")) or
            item["size"] > item["byte_limit"] or
            item["byte_limit"] > limit or
            not _is_sha256(item.get("sha256"))
        ):
            raise PinError(f"runtime file binding is invalid: {name}")
    thermal = bindings["thermal"]
    if (
        not isinstance(thermal, dict) or
        set(thermal) != {"path", "device", "inode", "mode"} or
        type(thermal.get("path")) is not str or
        not Path(thermal["path"]).is_absolute() or
        any(type(thermal.get(field)) is not int or thermal[field] < 0
            for field in ("device", "inode", "mode"))
    ):
        raise PinError("thermal runtime binding is invalid")
    return bindings


def _pin_projection(bindings):
    _validate_runtime_binding_layout(bindings)
    return {
        name: {
            "class": _file_class(name),
            "sha256": item["sha256"],
            "size": item["size"],
            "byte_limit": item["byte_limit"],
        }
        for name, item in sorted(_runtime_files(bindings).items())
    }


def make_pin_record(manifest):
    """Derive the result-free external pin record from a frozen manifest."""
    if (
        not isinstance(manifest, dict) or
        manifest.get("schema") != CAMPAIGN_SCHEMA or
        not _is_sha256(manifest.get("frozen_roster_sha256")) or
        not _is_sha256(manifest.get("policy_sha256")) or
        not _is_sha256(manifest.get("result_free_plan_sha256")) or
        not isinstance(manifest.get("frozen_roster"), dict) or
        manifest["frozen_roster"].get("native_protocol") !=
            REPAIRTIMING_PROTOCOL or
        manifest["frozen_roster"].get("native_schema") !=
            REPAIRTIMING_SCHEMA
    ):
        raise PinError("manifest cannot produce a v3 pin record")
    return {
        "schema": PIN_SCHEMA,
        "campaign_schema": CAMPAIGN_SCHEMA,
        "repairtiming_protocol": REPAIRTIMING_PROTOCOL,
        "repairtiming_schema": REPAIRTIMING_SCHEMA,
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "policy_sha256": manifest["policy_sha256"],
        "result_free_plan_sha256": manifest["result_free_plan_sha256"],
        "runtime_files": _pin_projection(manifest["runtime_bindings"]),
        "thermal": dict(manifest["runtime_bindings"]["thermal"]),
    }


def _validate_pin_record(pin, manifest):
    expected = make_pin_record(manifest)
    if (
        type(pin) is not dict or
        canonical_json_bytes(pin) != canonical_json_bytes(expected)
    ):
        raise PinError("pin record does not match the frozen manifest")
    return pin


def _verify_runtime_files(bindings, *, retain_python=False):
    """Independently reproduce every manifest file binding."""
    _validate_runtime_binding_layout(bindings)
    snapshots = {}
    for name, expected in _runtime_files(bindings).items():
        snapshot = _stable_file(
            expected["path"],
            byte_limit=FILE_BYTE_LIMITS[_file_class(name)],
            retain_payload=retain_python and name in (
                "runner", "parser", "context_tool"),
        )
        actual = {
            "path": snapshot["path"],
            "device": snapshot["device"],
            "inode": snapshot["inode"],
            "mode": snapshot["mode"],
            "size": snapshot["size"],
            "mtime_ns": snapshot["mtime_ns"],
            "ctime_ns": snapshot["ctime_ns"],
            "sha256": snapshot["sha256"],
            "byte_limit": expected["byte_limit"],
        }
        if actual != expected:
            raise PinError(f"runtime file changed: {name}")
        if name == "benchmark" and snapshot["mode"] & 0o111 == 0:
            raise PinError("benchmark binding is not executable")
        snapshots[name] = snapshot
    thermal = bindings["thermal"]
    try:
        current = os.stat(thermal["path"], follow_symlinks=False)
    except OSError as error:
        raise PinError(f"thermal binding disappeared: {error}")
    if {
        "path": thermal["path"],
        "device": current.st_dev,
        "inode": current.st_ino,
        "mode": current.st_mode,
    } != thermal or not stat.S_ISREG(current.st_mode):
        raise PinError("thermal binding changed")
    return snapshots


def _module_from_runner_snapshot(snapshot, dependency_snapshots=None):
    path = Path(snapshot["path"])
    module_name = (
        "_wirehair_pinned_rv4a_"
        + snapshot["sha256"][:16]
    )
    cached = _LOADED_CAMPAIGNS.get(module_name)
    if cached is not None:
        if cached["snapshot"]["sha256"] != snapshot["sha256"]:
            raise PinError("cached runner source attestation changed")
        return cached["module"]
    occupied = [
        name for name in (
            module_name, "peel_codec", "repair_timing_codec")
        if name in sys.modules
    ]
    if occupied:
        raise PinError(
            "refusing preloaded campaign dependency modules: "
            + ",".join(occupied))
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    module._BOOTSTRAP_RUNNER_SOURCE_SHA256 = snapshot["sha256"]
    module.__rv4a_source_sha256__ = snapshot["sha256"]
    if dependency_snapshots is not None:
        if set(dependency_snapshots) != {"parser", "context_tool"}:
            raise PinError("runner dependency snapshots are incomplete")
        module._BOOTSTRAP_SOURCE_PINS = {
            "repair_timing_codec":
                dependency_snapshots["parser"]["sha256"],
            "peel_codec":
                dependency_snapshots["context_tool"]["sha256"],
        }
    sys.modules[module_name] = module
    try:
        code = compile(
            snapshot["payload"], str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
    except BaseException:
        for name in (
                module_name, "peel_codec", "repair_timing_codec"):
            sys.modules.pop(name, None)
        raise
    _LOADED_CAMPAIGNS[module_name] = {
        "module": module,
        "snapshot": {
            key: value for key, value in snapshot.items()
            if key != "payload"
        },
    }
    return module


def _assert_loaded_sources(campaign, snapshots, bindings):
    expected = {
        "runner": campaign,
        "parser": getattr(campaign, "repair", None),
        "context_tool": getattr(campaign, "peel_codec", None),
    }
    module_names = {
        "parser": "repair_timing_codec",
        "context_tool": "peel_codec",
    }
    for name, module in expected.items():
        if (
            not isinstance(module, types.ModuleType) or
            Path(module.__file__) != Path(bindings["sources"][name]["path"]) or
            getattr(module, "__rv4a_source_sha256__", None) !=
                snapshots[name]["sha256"]
        ):
            raise PinError(f"campaign did not source-load pinned {name}")
        if name != "runner" and sys.modules.get(module_names[name]) is not module:
            raise PinError(f"campaign substituted pinned {name}")
        if snapshots[name]["sha256"] != \
                bindings["sources"][name]["sha256"]:
            raise PinError(f"executed source hash changed: {name}")
    if (
        getattr(campaign, "CAMPAIGN_SCHEMA", None) != CAMPAIGN_SCHEMA or
        getattr(campaign, "REQUIRED_REPAIRTIMING_PROTOCOL", None) !=
            REPAIRTIMING_PROTOCOL or
        getattr(campaign, "REQUIRED_REPAIRTIMING_SCHEMA", None) !=
            REPAIRTIMING_SCHEMA
    ):
        raise PinError("loaded campaign protocol is not frozen v3")


def make_preflight_pin(runner_path, bench_path, thermal_source):
    """Source-load the prospective runner and pin its complete runtime."""
    runner_path = Path(runner_path).absolute()
    repository = runner_path.parent.parent
    snapshots = {
        "runner": _stable_file(
            runner_path,
            byte_limit=FILE_BYTE_LIMITS["python"],
            retain_payload=True,
        ),
        "parser": _stable_file(
            repository / "tools/repair_timing_codec.py",
            byte_limit=FILE_BYTE_LIMITS["python"],
            retain_payload=True,
        ),
        "context_tool": _stable_file(
            repository / "tools/peel_codec.py",
            byte_limit=FILE_BYTE_LIMITS["python"],
            retain_payload=True,
        ),
    }
    campaign = _module_from_runner_snapshot(
        snapshots["runner"], {
            "parser": snapshots["parser"],
            "context_tool": snapshots["context_tool"],
        })
    bindings = campaign.runtime_bindings(bench_path, thermal_source)
    independently_replayed = _verify_runtime_files(
        bindings, retain_python=True)
    _assert_loaded_sources(campaign, independently_replayed, bindings)
    record = campaign.make_runtime_pin_record(bindings)
    validation = campaign.validate_frozen_contract()
    expected = {
        "schema": PIN_SCHEMA,
        "campaign_schema": CAMPAIGN_SCHEMA,
        "repairtiming_protocol": REPAIRTIMING_PROTOCOL,
        "repairtiming_schema": REPAIRTIMING_SCHEMA,
        "frozen_roster_sha256":
            validation["frozen_roster_sha256"],
        "policy_sha256":
            validation["policy_sha256"],
        "result_free_plan_sha256":
            validation["result_free_plan_sha256"],
        "runtime_files": _pin_projection(bindings),
        "thermal": dict(bindings["thermal"]),
    }
    if canonical_json_bytes(record) != canonical_json_bytes(expected):
        raise PinError("prospective runner generated a mismatched pin")
    return record


def verify_with_pins(
        directory, pin_path, expected_pin_sha256, *,
        replay_workers=MAX_REPLAY_WORKERS):
    """Strictly reproduce one completion with exact external source pins."""
    if (
        not _is_sha256(expected_pin_sha256) or
        type(replay_workers) is not int or
        not 1 <= replay_workers <= MAX_REPLAY_WORKERS
    ):
        raise PinError("trusted pin digest or replay worker count is invalid")
    directory = Path(directory).absolute()
    manifest, unused_manifest_sha256 = _read_canonical(
        directory / "manifest.json", byte_limit=MANIFEST_BYTE_LIMIT)
    pin, pin_sha256 = _read_canonical(
        pin_path, byte_limit=PIN_BYTE_LIMIT)
    if (
        pin_sha256 != expected_pin_sha256 or
        manifest.get("preflight_pin_sha256") != pin_sha256
    ):
        raise PinError("pin record does not match its trusted SHA-256")
    _validate_pin_record(pin, manifest)
    bindings = manifest["runtime_bindings"]
    snapshots = _verify_runtime_files(bindings, retain_python=True)
    verifier_binding = bindings["sources"]["parallel_verifier"]
    if (
        Path(__file__).absolute() != Path(verifier_binding["path"]) or
        not _is_sha256(_EXECUTING_VERIFIER_SOURCE_SHA256) or
        _EXECUTING_VERIFIER_SOURCE_SHA256 != verifier_binding["sha256"]
    ):
        raise PinError("executing verifier is not manifest-bound")
    campaign = _module_from_runner_snapshot(
        snapshots["runner"], {
            "parser": snapshots["parser"],
            "context_tool": snapshots["context_tool"],
        })
    _assert_loaded_sources(campaign, snapshots, bindings)
    try:
        campaign._validate_manifest(manifest)
        campaign.check_runtime_bindings(bindings, full_hash=True)
        verified = campaign.verify_campaign(
            directory, replay_workers=replay_workers)
        campaign.check_runtime_bindings(bindings, full_hash=True)
    except Exception as error:
        if isinstance(error, PinError):
            raise
        raise PinError(f"frozen campaign replay failed: {error}")
    after = _verify_runtime_files(bindings, retain_python=True)
    _assert_loaded_sources(campaign, after, bindings)
    return {
        "schema": "wirehair.wh2.rv4a.pinned-verification.v1",
        "pin_sha256": pin_sha256,
        "runtime_pin_sha256":
            sha256_bytes(canonical_json_bytes({
                "runtime_files": pin["runtime_files"],
                "thermal": pin["thermal"],
            })),
        "verified": verified,
    }


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    make = subparsers.add_parser(
        "make-pin",
        help="print the canonical result-free pin record for a manifest",
    )
    make.add_argument("--manifest", required=True, type=Path)
    preflight = subparsers.add_parser(
        "preflight-pin",
        help="pin a complete prospective runtime before any workers",
    )
    preflight.add_argument("--runner", required=True, type=Path)
    preflight.add_argument("--bench", required=True, type=Path)
    preflight.add_argument(
        "--thermal-source", required=True, type=Path)
    verify = subparsers.add_parser(
        "verify",
        help="source-pin and exactly replay a completed campaign",
    )
    verify.add_argument("--directory", required=True, type=Path)
    verify.add_argument("--pin-record", required=True, type=Path)
    verify.add_argument("--pin-sha256", required=True)
    verify.add_argument(
        "--replay-workers", type=int, default=MAX_REPLAY_WORKERS)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.command == "make-pin":
        manifest, unused_digest = _read_canonical(
            args.manifest, byte_limit=MANIFEST_BYTE_LIMIT)
        print(
            canonical_json_bytes(make_pin_record(manifest)).decode("ascii"),
            end="",
        )
    elif args.command == "preflight-pin":
        print(
            canonical_json_bytes(make_preflight_pin(
                args.runner,
                args.bench,
                args.thermal_source,
            )).decode("ascii"),
            end="",
        )
    else:
        print(
            canonical_json_bytes(verify_with_pins(
                args.directory,
                args.pin_record,
                args.pin_sha256,
                replay_workers=args.replay_workers,
            )).decode("ascii"),
            end="",
        )


def _source_forced_cli(argv):
    """Run the verifier from one retained exact self-source snapshot."""
    path = Path(__file__).absolute()
    snapshot = _stable_file(
        path,
        byte_limit=FILE_BYTE_LIMITS["python"],
        retain_payload=True,
    )
    digest = snapshot["sha256"]
    module_name = f"_wirehair_rv4a_verifier_cli_{digest[:16]}"
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    module._EXECUTING_VERIFIER_SOURCE_SHA256 = digest
    sys.modules[module_name] = module
    try:
        code = compile(
            snapshot["payload"], str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
        current = module._stable_file(
            path,
            byte_limit=module.FILE_BYTE_LIMITS["python"],
            retain_payload=False,
        )
        if current["sha256"] != digest:
            raise RuntimeError(
                "parallel verifier changed during source-forced launch")
        module.main(argv)
    except Exception as error:
        if isinstance(error, (
                getattr(module, "PinError", RuntimeError),
                ValueError,
        )):
            print(f"error: {error}", file=sys.stderr)
            return 2
        raise
    finally:
        sys.modules.pop(module_name, None)
    return 0


if __name__ == "__main__":
    raise SystemExit(_source_forced_cli(sys.argv[1:]))
