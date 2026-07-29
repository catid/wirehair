#!/usr/bin/env python3
"""Build and authenticate an exact WH2 hook-on/off A/B matrix.

This tool only builds binaries and records provenance.  It deliberately does
not execute the timing target or interpret benchmark output.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import subprocess
import sys


SCHEMA = "wirehair.wh2.hook-build-matrix.v1"
HOOK_MACRO = "WIREHAIR_V2_ENABLE_TEST_HOOKS"
MATRIX = (
    ("hook-on-a", "on", "a"),
    ("hook-off-a", "off", "a"),
    ("hook-on-b", "on", "b"),
    ("hook-off-b", "off", "b"),
)
MAX_EVIDENCE_BYTES = 128 * 1024 * 1024
RESERVED_DEFINITIONS = frozenset({
    "BUILD_CODEC_V2",
    "BUILD_TESTS",
    "BUILD_SHARED_LIBS",
    "CMAKE_BUILD_TYPE",
    "CMAKE_C_FLAGS",
    "CMAKE_CXX_FLAGS",
    "CMAKE_EXE_LINKER_FLAGS",
    "CMAKE_EXPORT_COMPILE_COMMANDS",
    "CMAKE_FIND_USE_PACKAGE_REGISTRY",
    "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY",
    "CMAKE_MODULE_LINKER_FLAGS",
    "CMAKE_SHARED_LINKER_FLAGS",
    "MARCH_NATIVE",
    "WH_LTO",
    "WH_PGO_DIR",
    "WH_PGO_MODE",
    "WIREHAIR_BUILD_BENCHMARKS",
    "WIREHAIR_STATIC_PIC",
})
RETAINED_BUILD_ENVIRONMENT_KEYS = ("HOME", "PATH", "TEMP", "TMP", "TMPDIR")
CMAKE_TOOL_KEYS = (
    "CMAKE_AR", "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER",
    "CMAKE_C_COMPILER_AR", "CMAKE_C_COMPILER_RANLIB",
    "CMAKE_CXX_COMPILER_AR", "CMAKE_CXX_COMPILER_RANLIB",
    "CMAKE_LINKER", "CMAKE_MAKE_PROGRAM", "CMAKE_NM",
    "CMAKE_OBJDUMP", "CMAKE_RANLIB", "CMAKE_STRIP",
)
class MatrixError(RuntimeError):
    """Raised when build provenance or comparability cannot be proven."""


def canonical_json_bytes(value):
    return (
        json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        ) + "\n"
    ).encode("ascii")


def canonical_sha256(value):
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _is_within(path, parent):
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _stable_bytes(path, *, byte_limit=MAX_EVIDENCE_BYTES):
    # Do not resolve the final component before O_NOFOLLOW; doing so would
    # silently turn a symlink substitution into an ordinary-file read.
    path = Path(os.path.abspath(os.fspath(path)))
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
        try:
            before = os.fstat(descriptor)
            if not stat.S_ISREG(before.st_mode):
                raise MatrixError(f"not a regular file: {path}")
            if before.st_size > byte_limit:
                raise MatrixError(f"file exceeds evidence limit: {path}")
            chunks = []
            size = 0
            while True:
                chunk = os.read(descriptor, min(1024 * 1024, byte_limit + 1))
                if not chunk:
                    break
                size += len(chunk)
                if size > byte_limit:
                    raise MatrixError(f"file exceeds evidence limit: {path}")
                chunks.append(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise MatrixError(f"could not read {path}: {error}")
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    if (
        any(getattr(before, key) != getattr(after, key) for key in stable) or
        any(getattr(after, key) != getattr(current, key) for key in stable)
    ):
        raise MatrixError(f"file changed while read: {path}")
    payload = b"".join(chunks)
    if len(payload) != after.st_size:
        raise MatrixError(f"short read: {path}")
    return payload


def file_binding(path, *, byte_limit=MAX_EVIDENCE_BYTES):
    payload = _stable_bytes(path, byte_limit=byte_limit)
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _run(command, *, cwd, environment, timeout, output_limit):
    command = [str(part) for part in command]
    try:
        result = subprocess.run(
            command,
            cwd=str(cwd),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise MatrixError(f"command could not run: {command!r}: {error}")
    if (
        len(result.stdout) > output_limit or
        len(result.stderr) > output_limit
    ):
        raise MatrixError(f"command output exceeded limit: {command!r}")
    if result.returncode != 0:
        stderr = result.stderr[-4096:].decode("utf-8", "replace")
        raise MatrixError(
            f"command failed with {result.returncode}: {command!r}: "
            f"{stderr}")
    return {
        "command": command,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def _output_binding(payload):
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _command_receipt(result):
    return {
        "command": result["command"],
        "stdout": _output_binding(result["stdout"]),
        "stderr": _output_binding(result["stderr"]),
    }


def _git(repository, arguments, *, timeout=60):
    environment, unused_evidence = _minimal_environment()
    environment.update({
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_COUNT": "2",
        "GIT_CONFIG_KEY_0": "core.fsmonitor",
        "GIT_CONFIG_VALUE_0": "false",
        "GIT_CONFIG_KEY_1": "core.untrackedCache",
        "GIT_CONFIG_VALUE_1": "false",
        "GIT_NO_REPLACE_OBJECTS": "1",
        "LC_ALL": "C",
        "LANG": "C",
    })
    return _run(
        ["git", "-C", str(repository), *arguments],
        cwd=repository,
        environment=environment,
        timeout=timeout,
        output_limit=MAX_EVIDENCE_BYTES,
    )["stdout"]


def _decode_hex_identifier(payload, length, description):
    value = payload.decode("ascii").strip()
    if (
        len(value) != length or
        any(character not in "0123456789abcdef" for character in value)
    ):
        raise MatrixError(f"invalid {description}: {value!r}")
    return value


def _parse_ls_tree(payload):
    entries = []
    seen = set()
    for raw in payload.split(b"\0"):
        if not raw:
            continue
        try:
            header, raw_path = raw.split(b"\t", 1)
            mode, kind, object_id = header.decode("ascii").split(" ")
            path = raw_path.decode("utf-8")
        except (UnicodeDecodeError, ValueError) as error:
            raise MatrixError(f"malformed git tree entry: {error}")
        relative = Path(path)
        if (
            kind != "blob" or
            mode not in ("100644", "100755") or
            not re.fullmatch(r"[0-9a-f]{40}", object_id) or
            relative.is_absolute() or
            ".." in relative.parts or
            path in seen
        ):
            raise MatrixError(f"unsupported git tree entry: {path!r}")
        seen.add(path)
        entries.append({
            "path": path,
            "mode": mode,
            "git_object": object_id,
        })
    if not entries:
        raise MatrixError("source commit has no tracked files")
    return entries


def _worktree_payload(repository, entry):
    path = repository / entry["path"]
    try:
        current = os.lstat(path)
    except OSError as error:
        raise MatrixError(f"could not stat tracked source {path}: {error}")
    if not stat.S_ISREG(current.st_mode):
        raise MatrixError(f"tracked regular file changed type: {path}")
    executable = bool(current.st_mode & (
        stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH))
    if executable != (entry["mode"] == "100755"):
        raise MatrixError(f"tracked executable mode changed: {path}")
    return _stable_bytes(path)


def source_receipt(repository):
    repository = Path(repository).resolve()
    top = Path(
        _git(repository, ["rev-parse", "--show-toplevel"])
        .decode("utf-8").strip()
    ).resolve()
    if top != repository:
        raise MatrixError(
            f"source must be the repository root: {repository} != {top}")
    status = _git(
        repository,
        ["status", "--porcelain=v2", "-z", "--untracked-files=all"],
    )
    if status:
        raise MatrixError(
            "source repository is dirty or contains untracked files")
    ignored = _git(
        repository,
        ["ls-files", "-z", "--others", "--ignored", "--exclude-standard"],
    )
    if ignored:
        raise MatrixError(
            "source repository contains ignored untracked files; use a "
            "fresh detached worktree")
    commit = _decode_hex_identifier(
        _git(repository, ["rev-parse", "--verify", "HEAD^{commit}"]),
        40, "commit",
    )
    tree = _decode_hex_identifier(
        _git(repository, ["rev-parse", "--verify", f"{commit}^{{tree}}"]),
        40, "tree",
    )
    entries = _parse_ls_tree(_git(
        repository, ["ls-tree", "-r", "-z", "--full-tree", commit]))

    # Hash worktree files through Git in bounded groups and compare the object
    # IDs with the committed tree.  This catches assume-unchanged/index flags
    # that a porcelain cleanliness check alone can miss.
    for offset in range(0, len(entries), 100):
        group = entries[offset:offset + 100]
        objects = _git(
            repository,
            ["hash-object", "--no-filters", "--",
             *[entry["path"] for entry in group]],
        ).decode("ascii").splitlines()
        if len(objects) != len(group):
            raise MatrixError("git hash-object returned the wrong row count")
        for entry, object_id in zip(group, objects):
            if object_id != entry["git_object"]:
                raise MatrixError(
                    f"worktree differs from commit: {entry['path']}")

    files = []
    total_bytes = 0
    for entry in entries:
        payload = _worktree_payload(repository, entry)
        object_payload = (
            b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload)
        if hashlib.sha1(object_payload).hexdigest() != entry["git_object"]:
            raise MatrixError(
                f"worktree payload differs from commit: {entry['path']}")
        total_bytes += len(payload)
        files.append({
            **entry,
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        })
    if (
        _git(
            repository,
            ["status", "--porcelain=v2", "-z", "--untracked-files=all"],
        ) or
        _decode_hex_identifier(
            _git(repository, ["rev-parse", "--verify", "HEAD^{commit}"]),
            40, "commit",
        ) != commit
    ):
        raise MatrixError("source repository changed while pinned")
    receipt = {
        "repository": str(repository),
        "commit": commit,
        "tree": tree,
        "tracked_files": files,
        "tracked_file_count": len(files),
        "tracked_bytes": total_bytes,
    }
    receipt["tracked_manifest_sha256"] = canonical_sha256(files)
    return receipt


def _parse_definition(value):
    if "\0" in value or "\n" in value or "=" not in value:
        raise MatrixError(
            f"CMake definition must be NAME=VALUE: {value!r}")
    name, setting = value.split("=", 1)
    if (
        not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name) or
        name in RESERVED_DEFINITIONS or name.startswith("CMAKE_") or
        HOOK_MACRO in setting
    ):
        raise MatrixError(f"reserved or invalid CMake definition: {value!r}")
    raise MatrixError(
        "v1 does not admit custom CMake definitions because they can add "
        f"unbound build inputs: {name}")


def _validate_flag_string(value, description):
    if "\0" in value or "\n" in value or HOOK_MACRO in value:
        raise MatrixError(
            f"{description} contains invalid data or the hook macro")
    try:
        shlex.split(value)
    except ValueError as error:
        raise MatrixError(f"{description} is not shell-tokenizable: {error}")
    if value.strip():
        raise MatrixError(
            f"v1 requires empty custom {description}; use the pinned "
            "Release/MARCH_NATIVE/WH_LTO controls")
    return ""


def _state_flags(common, state):
    macro = f"-D{HOOK_MACRO}=1"
    return " ".join(part for part in (
        common,
        macro if state == "on" else "",
    ) if part)


def _normalize_text(value, *, source, build):
    value = value.replace(str(build), "$BUILD")
    value = value.replace(str(source), "$SOURCE")
    return value


def _macro_token(token):
    prefixes = ("-D", "/D", "-U", "/U")
    for prefix in prefixes:
        if token.startswith(prefix):
            body = token[len(prefix):]
            if body == HOOK_MACRO or body.startswith(HOOK_MACRO + "="):
                return prefix, body
    return None


def _response_file_token(token):
    if token.startswith("@"):
        return True
    for prefix in ("-Wl,", "-Wp,", "-Wa,"):
        if token.startswith(prefix):
            return any(
                segment.startswith("@")
                for segment in token[len(prefix):].split(",")
            )
    return False


def _normalize_tokens(tokens, *, state, source, build):
    normalized = []
    definitions = 0
    for token in tokens:
        if _response_file_token(token):
            raise MatrixError(
                "response-file compile commands cannot prove macro state")
        macro = _macro_token(token)
        if macro is not None:
            prefix, body = macro
            if prefix in ("-U", "/U"):
                raise MatrixError("hook macro is undefined in compile command")
            if body != f"{HOOK_MACRO}=1":
                raise MatrixError("hook macro has an unexpected value")
            definitions += 1
            continue
        if HOOK_MACRO in token:
            raise MatrixError("hook macro is embedded in an unknown token")
        normalized.append(
            _normalize_text(token, source=source, build=build))
    if state == "on" and definitions != 1:
        raise MatrixError(
            "each hook-on compile command must define the hook exactly once")
    if state == "off" and definitions != 0:
        raise MatrixError("hook-off compile command defines the hook")
    return normalized


def compile_commands_evidence(
        path, *, state, source, build, allowed_compilers=None,
        environment=None):
    payload = _stable_bytes(path)
    try:
        document = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise MatrixError(f"invalid compile_commands.json: {error}")
    if not isinstance(document, list) or not document:
        raise MatrixError("compile_commands.json is empty or malformed")
    projection = []
    source_paths = set()
    for row in document:
        if not isinstance(row, dict):
            raise MatrixError("compile command row is not an object")
        if ("arguments" in row) == ("command" in row):
            raise MatrixError(
                "compile command must contain exactly one command form")
        if "arguments" in row:
            tokens = row["arguments"]
            if (
                not isinstance(tokens, list) or not tokens or
                any(not isinstance(token, str) for token in tokens)
            ):
                raise MatrixError("compile arguments are malformed")
        else:
            if not isinstance(row["command"], str):
                raise MatrixError("compile command is malformed")
            try:
                tokens = shlex.split(row["command"])
            except ValueError as error:
                raise MatrixError(f"compile command cannot be split: {error}")
        directory = Path(row.get("directory", "")).resolve()
        if not _is_within(directory, build):
            raise MatrixError(
                f"compile directory escapes fresh build: {directory}")
        if allowed_compilers is not None:
            compiler = shutil.which(
                tokens[0],
                path=(environment or {}).get("PATH"),
            )
            if compiler is None:
                compiler = tokens[0]
            compiler_paths = {
                str(Path(value).resolve()) for value in (
                    compiler, tokens[0])
                if Path(value).is_absolute()
            }
            if not compiler_paths.intersection(allowed_compilers):
                raise MatrixError(
                    f"compile command uses an unbound launcher: {tokens[0]}")
        file_path = Path(row.get("file", ""))
        if not file_path.is_absolute():
            file_path = (directory / file_path).resolve()
        else:
            file_path = file_path.resolve()
        if not _is_within(file_path, source):
            raise MatrixError(
                f"compile source is outside pinned repository: {file_path}")
        relative_source = file_path.relative_to(source).as_posix()
        source_paths.add(relative_source)
        normalized = {
            "directory": _normalize_text(
                str(directory), source=source, build=build),
            "file": "$SOURCE/" + relative_source,
            "tokens": _normalize_tokens(
                tokens, state=state, source=source, build=build),
        }
        if "output" in row:
            if not isinstance(row["output"], str):
                raise MatrixError("compile output is malformed")
            output = Path(row["output"])
            if not output.is_absolute():
                output = (directory / output).resolve()
            else:
                output = output.resolve()
            if not _is_within(output, build):
                raise MatrixError(
                    f"compile output escapes fresh build: {output}")
            normalized["output"] = _normalize_text(
                str(output), source=source, build=build)
        projection.append(normalized)
    projection.sort(key=canonical_json_bytes)
    return {
        "file": {
            "path": str(Path(path).resolve()),
            **_output_binding(payload),
        },
        "entries": len(projection),
        "source_paths": sorted(source_paths),
        "normalized_sha256": canonical_sha256(projection),
        "normalized": projection,
    }


def _parse_cache(payload):
    values = {}
    for raw_line in payload.decode("utf-8").splitlines():
        if not raw_line or raw_line.startswith(("#", "//")):
            continue
        match = re.fullmatch(r"([^:=]+):([^=]+)=(.*)", raw_line)
        if match is None:
            raise MatrixError(f"malformed CMake cache line: {raw_line!r}")
        name, kind, value = match.groups()
        if name in values:
            raise MatrixError(f"duplicate CMake cache key: {name}")
        values[name] = {"type": kind, "value": value}
    return values


def cache_evidence(path, *, state, source, build, expected=None):
    payload = _stable_bytes(path)
    values = _parse_cache(payload)
    required = {
        "BUILD_CODEC_V2": "ON",
        "BUILD_TESTS": "OFF",
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_FIND_USE_PACKAGE_REGISTRY": "OFF",
        "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY": "OFF",
        "WH_PGO_MODE": "OFF",
        "WIREHAIR_BUILD_BENCHMARKS": "ON",
    }
    if expected:
        required.update(expected)
    required.update({
        "CMAKE_EXE_LINKER_FLAGS": "",
        "CMAKE_MODULE_LINKER_FLAGS": "",
        "CMAKE_SHARED_LINKER_FLAGS": "",
    })
    for name, expected in required.items():
        if values.get(name, {}).get("value") != expected:
            raise MatrixError(f"CMake cache does not pin {name}={expected}")
    forbidden_external_cache = (
        "CMAKE_C_COMPILER_LAUNCHER",
        "CMAKE_C_LINKER_LAUNCHER",
        "CMAKE_CXX_COMPILER_LAUNCHER",
        "CMAKE_CXX_LINKER_LAUNCHER",
        "CMAKE_PROJECT_INCLUDE",
        "CMAKE_PROJECT_INCLUDE_BEFORE",
        "CMAKE_PROJECT_TOP_LEVEL_INCLUDES",
        "CMAKE_TOOLCHAIN_FILE",
        "CMAKE_USER_MAKE_RULES_OVERRIDE",
        "CMAKE_USER_MAKE_RULES_OVERRIDE_C",
        "CMAKE_USER_MAKE_RULES_OVERRIDE_CXX",
    )
    for name, row in values.items():
        dynamic_project_include = re.fullmatch(
            r"CMAKE_PROJECT_.+_INCLUDE(?:_BEFORE)?", name)
        if (
            (name in forbidden_external_cache or dynamic_project_include) and
            row["value"]
        ):
            raise MatrixError(
                f"unbound external CMake input is forbidden: {name}")
    projection = {}
    flag_keys = {"CMAKE_C_FLAGS", "CMAKE_CXX_FLAGS"}
    macro_cache_definitions = 0
    for name, row in values.items():
        value = row["value"]
        if name in flag_keys:
            try:
                tokens = shlex.split(value)
            except ValueError as error:
                raise MatrixError(f"cache flag parsing failed: {error}")
            normalized = _normalize_tokens(
                tokens, state=state, source=source, build=build)
            if normalized:
                raise MatrixError(
                    f"v1 requires no custom global flags in {name}")
            if state == "on":
                macro_cache_definitions += 1
            value = shlex.join(normalized)
        elif HOOK_MACRO in value:
            raise MatrixError(
                f"hook macro escaped into unexpected cache key {name}")
        projection[name] = {
            "type": row["type"],
            "value": _normalize_text(
                value, source=source, build=build),
        }
    if state == "on" and macro_cache_definitions != len(flag_keys):
        raise MatrixError("hook macro is not pinned in both language flags")
    return {
        "file": {
            "path": str(Path(path).resolve()),
            **_output_binding(payload),
        },
        "normalized_sha256": canonical_sha256(projection),
        "normalized": projection,
        "raw": values,
    }


def _write_log(path, payload):
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    try:
        descriptor = os.open(str(path), flags, 0o644)
        try:
            view = memoryview(payload)
            while view:
                written = os.write(descriptor, view)
                if written <= 0:
                    raise OSError("short write")
                view = view[written:]
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except (OSError, MatrixError) as error:
        raise MatrixError(f"could not publish log {path}: {error}")
    return {"path": str(path.resolve()), **_output_binding(payload)}


def _atomic_write_receipt(path, payload):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    if temporary.exists():
        raise MatrixError(f"stale receipt temporary exists: {temporary}")
    try:
        _write_log(temporary, payload)
        os.link(temporary, path, follow_symlinks=False)
        os.unlink(temporary)
        directory_flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_DIRECTORY"):
            directory_flags |= os.O_DIRECTORY
        descriptor = os.open(str(path.parent), directory_flags)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except (OSError, MatrixError) as error:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise MatrixError(
            f"could not atomically publish receipt {path}: {error}")
    return {"path": str(path.resolve()), **_output_binding(payload)}


def _extract_link_evidence(build_result, build, binary, source, state):
    text = (build_result["stdout"] + b"\n" + build_result["stderr"]).decode(
        "utf-8", "replace")
    candidates = []
    binary_strings = (str(binary), binary.name)
    for line in text.splitlines():
        if not any(value in line for value in binary_strings):
            continue
        if " -c " in f" {line} " or " -o " not in f" {line} ":
            continue
        try:
            tokens = shlex.split(line)
        except ValueError:
            continue
        normalized = _normalize_tokens(
            tokens, state=state, source=source, build=build)
        candidates.append(normalized)
    if not candidates:
        for path in build.rglob("link.txt"):
            payload = _stable_bytes(path)
            text = payload.decode("utf-8", "replace")
            if not any(value in text for value in binary_strings):
                continue
            for line in text.splitlines():
                try:
                    tokens = shlex.split(line)
                except ValueError:
                    continue
                candidates.append(_normalize_tokens(
                    tokens, state=state, source=source, build=build))
    unique = {
        canonical_json_bytes(candidate): candidate
        for candidate in candidates
    }
    if len(unique) != 1:
        raise MatrixError(
            f"could not identify one exact link command for {binary}")
    command = next(iter(unique.values()))
    return {
        "normalized_sha256": canonical_sha256(command),
        "normalized_tokens": command,
    }


def _tool_receipt(path, *, cwd, environment, timeout):
    configured = str(path)
    resolved = shutil.which(configured, path=environment.get("PATH"))
    if resolved is None:
        candidate = Path(configured)
        if candidate.exists():
            resolved = str(candidate)
    if resolved is None:
        raise MatrixError(f"tool not found: {configured}")
    resolved_path = Path(resolved).resolve()
    before = file_binding(resolved_path)
    result = _run(
        [str(resolved_path), "--version"],
        cwd=cwd,
        environment=environment,
        timeout=min(timeout, 60),
        output_limit=1024 * 1024,
    )
    after = file_binding(resolved_path)
    if after != before:
        raise MatrixError(f"tool changed while versioned: {resolved_path}")
    return {
        "configured_path": configured,
        "resolved_path": str(resolved_path),
        **before,
        "version": {
            "stdout": {
                **_output_binding(result["stdout"]),
                "text": result["stdout"].decode("utf-8", "replace"),
            },
            "stderr": {
                **_output_binding(result["stderr"]),
                "text": result["stderr"].decode("utf-8", "replace"),
            },
        },
    }


def _require_tool_unchanged(
        description, path, expected, *, cwd, environment, timeout):
    current = _tool_receipt(
        path, cwd=cwd, environment=environment, timeout=timeout)
    if current != expected:
        raise MatrixError(f"{description} changed during build matrix")


def _build_id(
        binary, readelf, *, expected_binary, cwd, environment, timeout):
    readelf = Path(readelf).resolve()
    readelf_before = file_binding(readelf)
    binary_before = file_binding(binary)
    if binary_before != expected_binary:
        raise MatrixError("binary changed before build-ID extraction")
    result = _run(
        [str(readelf), "--wide", "--notes", str(binary)],
        cwd=cwd,
        environment=environment,
        timeout=min(timeout, 60),
        output_limit=16 * 1024 * 1024,
    )
    if (
        file_binding(readelf) != readelf_before or
        file_binding(binary) != binary_before
    ):
        raise MatrixError(
            "binary or readelf changed during build-ID extraction")
    text = result["stdout"].decode("utf-8", "replace")
    matches = re.findall(r"Build ID:\s*([0-9A-Fa-f]+)", text)
    if len(matches) != 1:
        raise MatrixError(f"binary lacks one unambiguous build ID: {binary}")
    return {
        "value": matches[0].lower(),
        "readelf_tool": {
            "path": str(readelf),
            **readelf_before,
        },
        "readelf_stdout": _output_binding(result["stdout"]),
        "readelf_stderr": _output_binding(result["stderr"]),
    }


def _minimal_environment():
    source = os.environ
    environment = {}
    for key in RETAINED_BUILD_ENVIRONMENT_KEYS:
        if key in source:
            if "\0" in source[key]:
                raise MatrixError(
                    f"environment variable contains NUL: {key}")
            environment[key] = source[key]
    if "PATH" not in environment:
        environment["PATH"] = os.defpath
    environment.update({"LC_ALL": "C", "LANG": "C", "TZ": "UTC"})
    retained = {
        key: environment[key]
        for key in RETAINED_BUILD_ENVIRONMENT_KEYS
        if key in environment
    }
    retained.update({"LANG": "C", "LC_ALL": "C", "TZ": "UTC"})
    return environment, {
        "policy": "minimal-allowlist-v1",
        "allowlist": list(RETAINED_BUILD_ENVIRONMENT_KEYS),
        "removed_count": sum(
            key not in environment for key in source),
        "retained": retained,
    }


def _sanitized_environment():
    return _minimal_environment()


def _build_one(
        name, state, replicate, *, source, build_root, target,
        binary_relative, cmake, readelf, generator, definitions,
        tracked_sources,
        c_flags, cxx_flags, exe_linker_flags, shared_linker_flags,
        module_linker_flags, build_shared, static_pic, march_native,
        lto, pgo, jobs, timeout, environment,
):
    build = build_root / name
    try:
        build.mkdir(mode=0o755)
    except OSError as error:
        raise MatrixError(f"could not create fresh build tree: {error}")
    cmake_before = _tool_receipt(
        cmake, cwd=source, environment=environment, timeout=timeout)
    readelf_before = _tool_receipt(
        readelf, cwd=source, environment=environment, timeout=timeout)
    configure = [
        str(cmake), "-S", str(source), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF",
        "-DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF",
        "-DBUILD_TESTS=OFF",
        "-DBUILD_CODEC_V2=ON",
        "-DWIREHAIR_BUILD_BENCHMARKS=ON",
        f"-DBUILD_SHARED_LIBS={build_shared}",
        f"-DWIREHAIR_STATIC_PIC={static_pic}",
        f"-DMARCH_NATIVE={march_native}",
        f"-DWH_LTO={lto}",
        f"-DWH_PGO_MODE={pgo}",
        f"-DCMAKE_C_FLAGS={_state_flags(c_flags, state)}",
        f"-DCMAKE_CXX_FLAGS={_state_flags(cxx_flags, state)}",
        f"-DCMAKE_EXE_LINKER_FLAGS={exe_linker_flags}",
        f"-DCMAKE_SHARED_LINKER_FLAGS={shared_linker_flags}",
        f"-DCMAKE_MODULE_LINKER_FLAGS={module_linker_flags}",
    ]
    if generator:
        configure.extend(["-G", generator])
    configure.extend(
        f"-D{name_}={value}" for name_, value in definitions)
    configure_result = _run(
        configure, cwd=source, environment=environment, timeout=timeout,
        output_limit=MAX_EVIDENCE_BYTES)
    configure_logs = {
        "stdout": _write_log(
            build / "matrix-configure.stdout", configure_result["stdout"]),
        "stderr": _write_log(
            build / "matrix-configure.stderr", configure_result["stderr"]),
    }
    expected_cache = {
        "BUILD_SHARED_LIBS": build_shared,
        "MARCH_NATIVE": march_native,
        "WH_LTO": lto,
        "WIREHAIR_STATIC_PIC": static_pic,
    }
    cache = cache_evidence(
        build / "CMakeCache.txt",
        state=state, source=source, build=build,
        expected=expected_cache)
    cache_tools_before = {
        key: _tool_receipt(
            cache["raw"][key]["value"],
            cwd=source, environment=environment, timeout=timeout)
        for key in CMAKE_TOOL_KEYS
        if key in cache["raw"] and cache["raw"][key]["value"]
    }
    if "CMAKE_CXX_COMPILER" not in cache_tools_before:
        raise MatrixError("CMake cache lacks the C++ compiler")
    if _tool_receipt(
            cmake, cwd=source, environment=environment,
            timeout=timeout) != cmake_before:
        raise MatrixError("CMake changed during configuration")

    build_command = [
        str(cmake), "--build", str(build), "--config", "Release",
        "--target", target, "--parallel", str(jobs), "--verbose",
    ]
    build_result = _run(
        build_command, cwd=source, environment=environment, timeout=timeout,
        output_limit=MAX_EVIDENCE_BYTES)
    build_logs = {
        "stdout": _write_log(
            build / "matrix-build.stdout", build_result["stdout"]),
        "stderr": _write_log(
            build / "matrix-build.stderr", build_result["stderr"]),
    }
    cache_after = cache_evidence(
        build / "CMakeCache.txt",
        state=state, source=source, build=build,
        expected=expected_cache)
    if (
        cache_after["file"] != cache["file"] or
        cache_after["normalized_sha256"] !=
            cache["normalized_sha256"]
    ):
        raise MatrixError("CMake cache changed during target build")
    cache_tools_after = {
        key: _tool_receipt(
            cache_after["raw"][key]["value"],
            cwd=source, environment=environment, timeout=timeout)
        for key in CMAKE_TOOL_KEYS
        if key in cache_after["raw"] and cache_after["raw"][key]["value"]
    }
    if cache_tools_after != cache_tools_before:
        raise MatrixError("compiler or build tool changed during target build")
    if _tool_receipt(
            cmake, cwd=source, environment=environment,
            timeout=timeout) != cmake_before:
        raise MatrixError("CMake changed during target build")
    if _tool_receipt(
            readelf, cwd=source, environment=environment,
            timeout=timeout) != readelf_before:
        raise MatrixError("readelf changed during target build")
    binary_path = build / binary_relative
    try:
        binary_stat = os.lstat(binary_path)
    except OSError as error:
        raise MatrixError(f"target binary is missing: {error}")
    if (
        not stat.S_ISREG(binary_stat.st_mode) or
        not binary_stat.st_mode & (
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    ):
        raise MatrixError("target artifact is not a regular executable")
    binary = binary_path.resolve()
    if not _is_within(binary, build):
        raise MatrixError("binary path escapes its fresh build tree")
    binary_binding = file_binding(binary)
    binary_evidence = {
        "path": str(binary),
        "mode": stat.S_IMODE(binary_stat.st_mode),
        **binary_binding,
        "build_id": _build_id(
            binary, readelf, expected_binary=binary_binding,
            cwd=source, environment=environment,
            timeout=timeout),
    }
    compile_evidence = compile_commands_evidence(
        build / "compile_commands.json",
        state=state, source=source, build=build,
        allowed_compilers={
            path
            for key in ("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER")
            if key in cache_tools_before
            for path in (
                cache_tools_before[key]["configured_path"],
                cache_tools_before[key]["resolved_path"],
            )
            if Path(path).is_absolute()
        },
        environment=environment)
    untracked_compile_sources = sorted(
        set(compile_evidence["source_paths"]) - set(tracked_sources))
    if untracked_compile_sources:
        raise MatrixError(
            "compile commands consume sources outside the pinned Git tree: "
            + ", ".join(untracked_compile_sources[:8]))
    link_evidence = _extract_link_evidence(
        build_result, build, binary, source, state)

    tools = {
        "CMAKE_COMMAND": cmake_before,
        "READELF_COMMAND": readelf_before,
        **cache_tools_before,
    }
    return {
        "name": name,
        "hook_state": state,
        "replicate": replicate,
        "build_directory": str(build),
        "configure": {
            **_command_receipt(configure_result),
            "logs": configure_logs,
        },
        "build": {
            **_command_receipt(build_result),
            "logs": build_logs,
        },
        "cache": {
            key: value for key, value in cache.items() if key != "raw"
        },
        "compile_commands": compile_evidence,
        "link": link_evidence,
        "tools": tools,
        "binary": binary_evidence,
    }


def validate_comparability(builds):
    names = [name for name, unused_state, unused_rep in MATRIX]
    if sorted(builds) != sorted(names):
        raise MatrixError("build matrix is incomplete")
    for name, state, replicate in MATRIX:
        if (
            builds[name].get("name") != name or
            builds[name].get("hook_state") != state or
            builds[name].get("replicate") != replicate
        ):
            raise MatrixError(f"build matrix identity changed: {name}")
    comparable_fields = (
        ("cache", "normalized_sha256"),
        ("compile_commands", "normalized_sha256"),
        ("link", "normalized_sha256"),
    )
    comparisons = {}
    for outer, inner in comparable_fields:
        values = {builds[name][outer][inner] for name in names}
        key = f"{outer}_{inner}_identical"
        comparisons[key] = len(values) == 1
        if not comparisons[key]:
            raise MatrixError(
                f"builds differ beyond the hook macro: {outer}.{inner}")
    tool_hashes = {
        canonical_sha256(builds[name]["tools"]) for name in names
    }
    comparisons["toolchains_identical"] = len(tool_hashes) == 1
    if not comparisons["toolchains_identical"]:
        raise MatrixError("builds used different compiler/tool provenance")
    for state in ("on", "off"):
        left = builds[f"hook-{state}-a"]["binary"]
        right = builds[f"hook-{state}-b"]["binary"]
        identical = (
            left["sha256"] == right["sha256"] and
            left["bytes"] == right["bytes"] and
            left.get("mode") == right.get("mode") and
            left["build_id"]["value"] == right["build_id"]["value"]
        )
        comparisons[f"hook_{state}_ab_binary_identical"] = identical
        if not identical:
            raise MatrixError(
                f"hook-{state} A/B binaries are not reproducible")
    comparisons["hook_on_off_binary_identical"] = (
        builds["hook-on-a"]["binary"]["sha256"] ==
        builds["hook-off-a"]["binary"]["sha256"]
    )
    if comparisons["hook_on_off_binary_identical"]:
        raise MatrixError(
            "hook-on/off binaries are identical; macro distinction was "
            "not linked into the requested target")
    comparisons["only_allowed_compile_difference"] = (
        f"-D{HOOK_MACRO}=1")
    return comparisons


def derive_matrix(args):
    source = Path(args.source).resolve()
    build_root = Path(args.build_root).resolve()
    output = Path(args.output).resolve()
    if not source.is_dir():
        raise MatrixError(f"source repository does not exist: {source}")
    if (
        args.target.startswith("-") or
        not re.fullmatch(r"[A-Za-z0-9_.+-]+", args.target)
    ):
        raise MatrixError("target must be a non-option CMake target name")
    if not 1 <= args.jobs <= 256:
        raise MatrixError("--jobs must be in [1, 256]")
    if not 1 <= args.timeout <= 86400:
        raise MatrixError("--timeout must be in [1, 86400]")
    if build_root.exists():
        raise MatrixError(
            f"fresh build root already exists: {build_root}")
    if _is_within(build_root, source) or _is_within(output, source):
        raise MatrixError("build root and receipt must be outside source")
    if output.exists():
        raise MatrixError(f"receipt already exists: {output}")
    if "\0" in args.binary_relative or "\n" in args.binary_relative:
        raise MatrixError("binary-relative contains invalid characters")
    binary_relative = Path(args.binary_relative)
    if (
        binary_relative.is_absolute() or
        ".." in binary_relative.parts or
        not binary_relative.parts
    ):
        raise MatrixError("binary-relative must be a safe relative path")
    definitions = []
    seen_definitions = set()
    for value in args.cmake_definition:
        parsed = _parse_definition(value)
        if parsed[0] in seen_definitions:
            raise MatrixError(f"duplicate CMake definition: {parsed[0]}")
        seen_definitions.add(parsed[0])
        definitions.append(parsed)
    flag_values = {
        "c_flags": _validate_flag_string(args.c_flags, "C flags"),
        "cxx_flags": _validate_flag_string(args.cxx_flags, "C++ flags"),
        "exe_linker_flags": _validate_flag_string(
            args.exe_linker_flags, "executable linker flags"),
        "shared_linker_flags": _validate_flag_string(
            args.shared_linker_flags, "shared linker flags"),
        "module_linker_flags": _validate_flag_string(
            args.module_linker_flags, "module linker flags"),
    }
    cmake = Path(shutil.which(args.cmake) or args.cmake).resolve()
    readelf = Path(shutil.which(args.readelf) or args.readelf).resolve()
    environment, removed_environment = _sanitized_environment()
    if args.pgo != "OFF":
        raise MatrixError(
            "v1 requires WH_PGO_MODE=OFF; profile state is not bound")
    if not sys.executable:
        raise MatrixError("Python executable path is unavailable")
    python = Path(os.path.abspath(sys.executable))
    python_before = _tool_receipt(
        python, cwd=source, environment=environment, timeout=args.timeout)
    git_before = _tool_receipt(
        "git", cwd=source, environment=environment, timeout=args.timeout)
    initial_source = source_receipt(source)
    _require_tool_unchanged(
        "Python", python, python_before,
        cwd=source, environment=environment, timeout=args.timeout)
    if _tool_receipt(
            "git", cwd=source, environment=environment,
            timeout=args.timeout) != git_before:
        raise MatrixError("Git changed while source was initially pinned")
    tracked_sources = {
        row["path"] for row in initial_source["tracked_files"]
    }
    try:
        build_root.mkdir(mode=0o755)
    except OSError as error:
        raise MatrixError(f"could not create fresh build root: {error}")
    builds = {}
    for name, state, replicate in MATRIX:
        builds[name] = _build_one(
            name, state, replicate,
            source=source,
            build_root=build_root,
            target=args.target,
            binary_relative=binary_relative,
            cmake=cmake,
            readelf=readelf,
            generator=args.generator,
            definitions=definitions,
            tracked_sources=tracked_sources,
            build_shared=args.build_shared,
            static_pic=args.static_pic,
            march_native=args.march_native,
            lto=args.lto,
            pgo=args.pgo,
            jobs=args.jobs,
            timeout=args.timeout,
            environment=environment,
            **flag_values,
        )
        if source_receipt(source) != initial_source:
            raise MatrixError(f"source changed during build {name}")
        _require_tool_unchanged(
            "Python", python, python_before,
            cwd=source, environment=environment, timeout=args.timeout)
        if _tool_receipt(
                "git", cwd=source, environment=environment,
                timeout=args.timeout) != git_before:
            raise MatrixError(f"Git changed during build {name}")
    comparisons = validate_comparability(builds)
    final_source = source_receipt(source)
    if final_source != initial_source:
        raise MatrixError("source changed across build matrix")
    if _tool_receipt(
            "git", cwd=source, environment=environment,
            timeout=args.timeout) != git_before:
        raise MatrixError("Git changed while source was finally pinned")
    try:
        helper_relative = \
            Path(__file__).resolve().relative_to(source).as_posix()
    except ValueError:
        raise MatrixError(
            "the executing helper is outside the pinned source repository")
    helper_rows = [
        row for row in initial_source["tracked_files"]
        if row["path"] == helper_relative
    ]
    if len(helper_rows) != 1:
        raise MatrixError(
            "build-matrix helper is not committed in the pinned source")
    top_tools = {
        "python": python_before,
        "git": git_before,
        "cmake": _tool_receipt(
            cmake, cwd=source, environment=environment,
            timeout=args.timeout),
        "readelf": _tool_receipt(
            readelf, cwd=source, environment=environment,
            timeout=args.timeout),
    }
    if _tool_receipt(
            "git", cwd=source, environment=environment,
            timeout=args.timeout) != git_before:
        raise MatrixError("Git changed during build matrix")
    for name, unused_state, unused_replicate in MATRIX:
        if (
            builds[name]["tools"]["CMAKE_COMMAND"] !=
                top_tools["cmake"] or
            builds[name]["tools"]["READELF_COMMAND"] !=
                top_tools["readelf"]
        ):
            raise MatrixError(
                f"build driver provenance changed for {name}")
    if source_receipt(source) != initial_source:
        raise MatrixError("source changed before receipt publication")
    _require_tool_unchanged(
        "Python", python, python_before,
        cwd=source, environment=environment, timeout=args.timeout)
    if _tool_receipt(
            "git", cwd=source, environment=environment,
            timeout=args.timeout) != git_before:
        raise MatrixError("Git changed before receipt publication")
    receipt = {
        "schema": SCHEMA,
        "source": initial_source,
        "provenance_helper": helper_rows[0],
        "hook_macro": {
            "name": HOOK_MACRO,
            "on_token": f"-D{HOOK_MACRO}=1",
            "off_contract": "macro-absent",
        },
        "configuration": {
            "target": args.target,
            "binary_relative": binary_relative.as_posix(),
            "build_type": "Release",
            "build_order": [
                name for name, unused_state, unused_replicate in MATRIX
            ],
            "generator": args.generator or "cmake-default",
            "definitions": [
                {"name": name, "value": value}
                for name, value in definitions
            ],
            "flags": flag_values,
            "build_shared": args.build_shared,
            "static_pic": args.static_pic,
            "march_native": args.march_native,
            "lto": args.lto,
            "pgo": args.pgo,
            "jobs": args.jobs,
            "build_environment":
                removed_environment,
        },
        "tools": top_tools,
        "builds": builds,
        "comparability": comparisons,
        "timing_executed": False,
    }
    payload = canonical_json_bytes(receipt)
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise MatrixError(f"could not create receipt directory: {error}")
    _atomic_write_receipt(output, payload)
    return receipt, hashlib.sha256(payload).hexdigest()


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).parents[1])
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--binary-relative", required=True)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--readelf", default="readelf")
    parser.add_argument("--generator")
    parser.add_argument(
        "--cmake-definition", action="append", default=[],
        metavar="NAME=VALUE")
    parser.add_argument("--c-flags", default="")
    parser.add_argument("--cxx-flags", default="")
    parser.add_argument("--exe-linker-flags", default="")
    parser.add_argument("--shared-linker-flags", default="")
    parser.add_argument("--module-linker-flags", default="")
    parser.add_argument(
        "--build-shared", choices=("ON", "OFF"), default="OFF")
    parser.add_argument(
        "--static-pic", choices=("ON", "OFF"), default="ON")
    parser.add_argument(
        "--march-native", choices=("ON", "OFF"), default="OFF")
    parser.add_argument(
        "--lto", choices=("OFF", "ON", "AUTO", "THIN"), default="OFF")
    parser.add_argument(
        "--pgo", choices=("OFF",), default="OFF")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    unused_receipt, digest = derive_matrix(args)
    print(canonical_json_bytes({
        "output": str(Path(args.output).resolve()),
        "schema": SCHEMA,
        "sha256": digest,
    }).decode("ascii"), end="")


if __name__ == "__main__":
    try:
        main()
    except MatrixError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
