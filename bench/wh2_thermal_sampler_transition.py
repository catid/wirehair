#!/usr/bin/env python3
"""Source/test-only V6 successor for the bounded WH2 thermal handoff.

This retired operational model remains specific to the reviewed g8iv.6.22.1
live dry so its recovery semantics and evidence contracts can be verified.
Its frozen old identity is historical and does not describe the separately
rescued sampler that now owns thermal monitoring.

V6 exposes no operational staging, preparation, live-helper, or execution
entry point.  Every former CLI mode and every public operational function is
hard-retired.  Explicitly synthetic lower models remain importable only for
source assurance and reject the production V6 namespace before mutation.
"""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import select
import signal
import stat
import subprocess
import sys
import time
from types import ModuleType
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

TRANSITION_ID = "wirehair-g8iv.6.22.1-live-dry-3c291790-v6"
EXECUTE_CONFIRMATION = "STOP-EXACT-OLD-SAMPLER-AND-RESTORE-IT"
SOURCE_ONLY_RETIREMENT = (
    "V6 is source/test-only: its frozen original identity does not match "
    "the live rescue-r1 sampler; staging, preparation, live execution, and "
    "privileged live helper modes are disabled")
SOURCE_TEST_TRANSITION_PREFIX = "wirehair-wh2-v6-source-test-"
TRANSITION_SCHEMA = "wirehair.wh2.thermal_transition.v8"
PHASE_SCHEMA = "wirehair.wh2.thermal_transition.phase.v8"
AUDIT_BINDING_SCHEMA = "wirehair.wh2.thermal_transition.audit_binding.v8"
PLAN_SCHEMA = "wirehair.wh2.thermal_transition.plan.v9"
CODE_SEAL_STAGE_SCHEMA = "wirehair.wh2.thermal_transition.code_stage.v6"
ROOT_CODE_STAGE_RECEIPT_SCHEMA = (
    "wirehair.wh2.thermal_transition.root_code_stage_receipt.v6")
CODE_SEAL_RECEIPT_SCHEMA = "wirehair.wh2.thermal_transition.code_seal.v6"
FD_INSPECTION_SCHEMA = "wirehair.wh2.thermal_transition.fd_inspection.v6"
IDENTITY_INSPECTION_SCHEMA = (
    "wirehair.wh2.thermal_transition.identity_inspection.v6")
PROCESS_SIGNAL_SCHEMA = (
    "wirehair.wh2.thermal_transition.process_signal.v6")
RECOVERY_ACTION_EVIDENCE_SCHEMA = (
    "wirehair.wh2.thermal_transition.recovery_action_evidence.v6")
# Hashes of the last reviewed V6 authority-shaped literals.  This tuple is
# evidence only: it contains no argv/program bytes, every public authority
# emitter is retired, and the private stage model rejects the production V6
# namespace.  Never interpret these digests as executable authorization.
RETIRED_NONEXECUTABLE_AUTHORITY_HASHES = (
    ("root_stage_perl_opener_sha256",
     "cbccd80040d8d3ea50f2404a6185b6e04836c8a652c1b03da825157a0c7a01ef"),
    ("root_stage_bash_program_sha256",
     "3f17a3fdb3e4dbdabcb2a8e03d27b02e33f1cc96bb330b188bef49e4b75be686"),
    ("root_stage_48_argv_canonical_sha256",
     "792ada77045bc7e971dfe42aee1520b0ff233bea8e6971cf5a8a2d4a64a88fcb"),
    ("root_stage_48_argv_nul_sha256",
     "bba1a91a967a9ef76ac32034184010a3503a62ccc1288284cbd7ace68a506113"),
    ("prepare_argv_nul_sha256",
     "e65d53b5697de4d379abd7b77cfbebf51e3c02743a58a807bf80b147a8734696"),
    ("identity_helper_argv_nul_sha256",
     "6261bae1a8f92fde22d518f57ea754f26831ea18fd98f960a951a88719a035d0"),
    ("signal_helper_argv_nul_sha256",
     "50ab01a730795d53a877ac9da1dac996d59bda27d90329e5d3c3f8c1ad936d7d"),
)
IDENTITY_PIDFD_POLICY = (
    "complete-roster-open-before-target-proc;"
    "two-stable-proc-observations;pidfd-unready;"
    "explicit-exit-stat-or-pidfd-ready-exiting-never-absence;"
    "only-enoent-esrch-pidfd-open-can-prove-absence")
PROCESS_SIGNAL_PIDFD_POLICY = (
    IDENTITY_PIDFD_POLICY + ";validated-target-pidfd-held-through-"
    "pidfd-send-signal;canonical-action-receipt")
PROCESS_IDENTITY_KEYS = frozenset({
    "affinity", "argv", "cmdline_sha256", "executable", "pid", "ppid",
    "process_group", "processor", "session", "start_tick", "state", "uids",
})
STABLE_PROCESS_IDENTITY_KEYS = PROCESS_IDENTITY_KEYS - {"processor", "state"}
PROC_STAT_IDENTITY_KEYS = frozenset({
    "ppid", "process_group", "processor", "session", "start_tick", "state",
})
EXIT_PROCESS_STATES = frozenset({"X", "x", "Z"})
RENAME_NOREPLACE = 1
SHA256_RE = re.compile(r"[0-9a-f]{64}")
SEMANTIC_FD_STATUS_MASK = 0
for _flag_name in (
        "O_APPEND", "O_ASYNC", "O_DIRECT", "O_DSYNC", "O_NONBLOCK",
        "O_PATH", "O_SYNC"):
    SEMANTIC_FD_STATUS_MASK |= getattr(os, _flag_name, 0)
del _flag_name

CANDIDATE_SHA256 = (
    "3c291790b6169c07e53dc0924d383455a81addd6da858d8d4d231b9f1263f6a2")
P32_SHA256 = (
    "414806a687d184fc2ec2973fcd0ea4b8ab232d13d8e707404e926d394fefb821")
OLD_SOURCE_SHA256 = (
    "2b84efa91375a96a4a64e09ce5bfd7cba0b85b75028f5a93470cd4ae58aadb01")
OLD_BOOT_ID = "1788608a-7aa1-48de-8f7c-848485be3cc3"
OLD_CMDLINE_SHA256 = (
    "9626d8e11e0e49a6de6bf81c3f3d227f108540a8eb61f40f212c58f0da9fd069")
OLD_LAUNCHER_CMDLINE_SHA256 = (
    "b5cbaa94f51dfe477f0c600e27acdfc0e396b9ec46ab24098bf2a27e1df7c86f")

REPO_CANDIDATE = Path(
    "/tmp/wirehair-wh2-thermal-plausibility-v2/bench/"
    "wirehair_expo_thermal_sampler.py")
REPO_CONTROLLER = Path(
    "/tmp/wirehair-wh2-thermal-plausibility-v2/bench/"
    "wh2_thermal_sampler_transition.py")
REPO_P32 = Path(
    "/tmp/wirehair-wh2-thermal-plausibility-v2/bench/"
    "wh2_p32_dispatch_timing.py")
DRY_ROOT = Path("/dev/shm") / TRANSITION_ID
CODE_SEAL_PARENT = Path("/dev/shm") / (TRANSITION_ID + "-root-code")
CODE_SEAL_ROOT = CODE_SEAL_PARENT / "seal-v6"
ROOT_CODE_STAGE_RECEIPT = CODE_SEAL_PARENT / (
    "root-code-seal-stage-receipt.v6.json")
OLD_ROOT = Path("/dev/shm/wh2-rhs-fusion-batch16-formal-v3")
OLD_SOURCE = OLD_ROOT / "campaign/frozen/wirehair_expo_thermal_sampler.py"
OLD_THERMAL_DIR = OLD_ROOT / "thermal"
OLD_CSV = OLD_THERMAL_DIR / "thermal.csv"
OLD_PID_FILE = OLD_THERMAL_DIR / "sampler.pid"
OLD_ARCHIVE = OLD_THERMAL_DIR / (
    "thermal.pre-g8iv.6.22.1.3c291790.v6.csv")
OLD_UNCLEAN_ARCHIVE = OLD_THERMAL_DIR / (
    "thermal.pre-g8iv.6.22.1.3c291790.v6.unclean.csv")
OLD_STALE_PID_ARCHIVE = OLD_THERMAL_DIR / (
    "sampler.pre-g8iv.6.22.1.3c291790.v6.unclean.pid")
AUDIT_BINDING = OLD_THERMAL_DIR / (
    "future-audit-binding.g8iv.6.22.1.3c291790.v6.json")

OLD_ARGV = (
    "/usr/bin/python3.12", "-I", str(OLD_SOURCE),
    "--csv", str(OLD_CSV), "--pid-file", str(OLD_PID_FILE),
)
OLD_LAUNCH_ARGV = (
    "sudo", "-n", "-b", "/usr/bin/taskset", "--cpu-list", "127",
    *OLD_ARGV,
)
SEALED_CONTROLLER = CODE_SEAL_ROOT / "wh2_thermal_sampler_transition.py"
SEALED_CANDIDATE = CODE_SEAL_ROOT / "wirehair_expo_thermal_sampler.py"
SEALED_P32 = CODE_SEAL_ROOT / "wh2_p32_dispatch_timing.py"
SEALED_LEGACY = CODE_SEAL_ROOT / "legacy_wirehair_expo_thermal_sampler.py"

OLD_REPLACEMENT_ARGV = (
    "/usr/bin/python3.12", "-I", "-S", "-B", str(SEALED_LEGACY),
    "--csv", str(OLD_CSV), "--pid-file", str(OLD_PID_FILE),
)
OLD_REPLACEMENT_CMDLINE_SHA256 = (
    "c04bcd2331f6d580d4ad7f3d0d27c0ba2692046395a8c0b67b1084b57ff19c8f")

EXECUTE_FLAG_CONTRACT: Dict[str, int | bool] = {
    "bytes_warning": 0,
    "debug": 0,
    "dev_mode": False,
    "dont_write_bytecode": 1,
    "hash_randomization": 1,
    "ignore_environment": 1,
    "inspect": 0,
    "int_max_str_digits": 4300,
    "interactive": 0,
    "isolated": 1,
    "no_site": 1,
    "no_user_site": 1,
    "optimize": 0,
    "quiet": 0,
    "safe_path": True,
    "utf8_mode": 1,
    "verbose": 0,
    "warn_default_encoding": 0,
}
PREPARE_FLAG_CONTRACT: Dict[str, int | bool] = dict(EXECUTE_FLAG_CONTRACT)

# Exact host executables retained by the historical source/test contracts.
# Synthetic preparation seals their inode metadata, while every operational
# V6 entry point remains retired.
TOOL_CONTRACT: Tuple[Tuple[str, str, str], ...] = (
    ("bash", "/usr/bin/bash",
     "bc5945feb8bd26203ebfafea5ce1878bb2e32cb8fb50ab7ae395cfb1e1aaaef1"),
    ("chmod", "/usr/bin/chmod",
     "f7ba478ce7ff158ccd3ba01bac1c4231e68cbd490fe49b5562a6ea008a586696"),
    ("env", "/usr/bin/env",
     "0aefff8f912fb75716c5d4de3b6acde93edbe8fa280fc8ee895c1226d3e373ef"),
    ("fuser", "/usr/bin/fuser",
     "88aeee250ef0622fb638acd4694a4ffd96ec04d5e25dd2ca880739c4719197f5"),
    ("install", "/usr/bin/install",
     "0e328ae109217200da3207ece12514b867d44fb90b444958b4d64b6007736f33"),
    ("mkdir", "/usr/bin/mkdir",
     "af2cc5106b37532b2aa1b4981668685cf7aebdf5fea3ac7356d568cf98d7470b"),
    ("mpstat", "/usr/bin/mpstat",
     "74695dd7d010730cd1e19efa0664904c1e3b17d746d304c06243e183a5ac9f9c"),
    ("perl", "/usr/bin/perl",
     "56e5ea41974eb1eff0f7ea64677578b1938053d29818c2810bcb21e2ca68cafa"),
    ("python3", "/usr/bin/python3.12",
     "1643dacd9feaedc58f3cc581e4d22577dfe25c09b10282936186ccf0f2e61118"),
    ("sudo", "/usr/bin/sudo",
     "136f2e48b0295b9fc595b8259cf2411ac43f27ddbfe02b956649ddaa2e92b9fa"),
    ("sha256sum", "/usr/bin/sha256sum",
     "9992e1f1feb6f0f396bc8d6691ebc1adbfc269fd628bce84eda1d4ba5c3995c7"),
    ("stat", "/usr/bin/stat",
     "3b87d297111f11d30b3c51fd2663f131a161e09d8e130e1842adaefb74307efe"),
    ("taskset", "/usr/bin/taskset",
     "a9c851792e54e91fba7b827019380abee54e715b6817899c835e4f221354b260"),
    ("timeout", "/usr/bin/timeout",
     "4fccd5b0192653a2446b745d5385ea547b78e466150e07ade9e2caff2b7f4e08"),
)

ROOT_STAGE_TOOL_NAMES = (
    "bash", "chmod", "env", "install", "mkdir", "perl", "sha256sum",
    "stat", "sudo", "timeout")

# This opener is historical program data exercised only by the private,
# production-rejecting synthetic stage model.  The Linux O_NOFOLLOW value is
# deliberately literal for the host fixture.  It opens all four sources before
# validating any of them, retains the handles, clears close-on-exec, and only
# then execs the synthetic-test Bash authority.
ROOT_STAGE_PERL_OPENER = r'''$open_flags = 131072;
@held = ();
for ($i = 0; $i < 4; ++$i) {
    sysopen(my $fh, $ARGV[$i], $open_flags) or
        die "source open failed errno=" . (0 + $!) . "\n";
    push @held, $fh;
}
for ($i = 0; $i < 4; ++$i) {
    @s = stat($held[$i]);
    die "source descriptor is not regular nlink1\n"
        unless @s && (($s[2] & 0170000) == 0100000) && $s[3] == 1;
    die "cannot retain source descriptor across exec\n"
        unless defined fcntl($held[$i], 2, 0);
}
@cmd = @ARGV[4 .. $#ARGV];
push @cmd, map { fileno($_) } @held;
exec {$cmd[0]} @cmd;
die "Bash authority exec failed\n";
'''

# The synthetic root-authority fixture uses only pinned, absolute tools.  It receives source FDs
# from ROOT_STAGE_PERL_OPENER, never opens a source pathname, copies only from
# /proc/self/fd/N, and deliberately leaves all residue on failure.  Its stdout
# is exactly the canonical receipt after final replay.
ROOT_STAGE_BASH_PROGRAM = r'''set -Eeuo pipefail
umask 077
die() { printf '%s\n' "$1" >&2; exit 1; }
[[ $# -eq 25 ]] || die 'stage argument roster changed'
transition_id=$1
anchor=$2
parent=$3
root=$4
receipt=$5
shift 5
owner_uid=$1
owner_gid=$2
shift 2
names=(candidate controller legacy p32)
sources=()
destinations=()
expected_hashes=()
for ((i = 0; i < 4; ++i)); do
    sources+=("$1")
    destinations+=("$2")
    expected_hashes+=("$3")
    shift 3
done
perl_program_sha256=$1
bash_program_sha256=$2
shift 2
fds=("$@")
[[ ${#fds[@]} -eq 4 ]] || die 'source FD roster changed'
[[ $transition_id =~ ^[A-Za-z0-9._-]+$ ]] || die 'transition id is unsafe'
[[ $owner_uid =~ ^(0|[1-9][0-9]*)$ && $owner_uid -le 4294967294 && \
   $owner_gid =~ ^(0|[1-9][0-9]*)$ && $owner_gid -le 4294967294 ]] || \
    die 'stage owner identity is malformed'
for path in "$anchor" "$parent" "$root" "$receipt" \
        "${sources[@]}" "${destinations[@]}"; do
    [[ $path =~ ^/[A-Za-z0-9._/-]+$ ]] || die 'stage path is unsafe'
done
[[ $anchor == /dev/shm ]] || die 'stage anchor changed'
[[ $parent == "$anchor/$transition_id-root-code" ]] || \
    die 'stage parent ancestry changed'
[[ $root == "$parent/seal-v6" ]] || die 'stage root ancestry changed'
[[ $receipt == "$parent/root-code-seal-stage-receipt.v6.json" ]] || \
    die 'stage receipt path changed'

tool_names=(bash chmod env install mkdir perl sha256sum stat sudo timeout)
tool_paths=(/usr/bin/bash /usr/bin/chmod /usr/bin/env /usr/bin/install \
    /usr/bin/mkdir /usr/bin/perl /usr/bin/sha256sum /usr/bin/stat \
    /usr/bin/sudo /usr/bin/timeout)
tool_hashes=(
    bc5945feb8bd26203ebfafea5ce1878bb2e32cb8fb50ab7ae395cfb1e1aaaef1
    f7ba478ce7ff158ccd3ba01bac1c4231e68cbd490fe49b5562a6ea008a586696
    0aefff8f912fb75716c5d4de3b6acde93edbe8fa280fc8ee895c1226d3e373ef
    0e328ae109217200da3207ece12514b867d44fb90b444958b4d64b6007736f33
    af2cc5106b37532b2aa1b4981668685cf7aebdf5fea3ac7356d568cf98d7470b
    56e5ea41974eb1eff0f7ea64677578b1938053d29818c2810bcb21e2ca68cafa
    9992e1f1feb6f0f396bc8d6691ebc1adbfc269fd628bce84eda1d4ba5c3995c7
    3b87d297111f11d30b3c51fd2663f131a161e09d8e130e1842adaefb74307efe
    136f2e48b0295b9fc595b8259cf2411ac43f27ddbfe02b956649ddaa2e92b9fa
    4fccd5b0192653a2446b745d5385ea547b78e466150e07ade9e2caff2b7f4e08)
chmod_tool=${tool_paths[1]}
install_tool=${tool_paths[3]}
mkdir_tool=${tool_paths[4]}
sha_tool=${tool_paths[6]}
stat_tool=${tool_paths[7]}

tool_binding() {
    local path=$1 expected=$2 line kind uid gid mode nlink observed
    line=$("$stat_tool" -c '%F|%u|%g|%a|%h' -- "$path") || \
        die 'tool stat failed'
    IFS='|' read -r kind uid gid mode nlink <<<"$line"
    [[ $kind == 'regular file' && $uid == 0 && $gid == 0 && $nlink -ge 1 ]] || \
        die 'tool binding changed'
    (( (8#$mode & 8#22) == 0 && (8#$mode & 8#111) != 0 )) || \
        die 'tool mode changed'
    line=$("$sha_tool" -- "$path") || die 'tool hash failed'
    observed=${line%% *}
    [[ $observed == "$expected" ]] || die 'tool hash changed'
}
for ((i = 0; i < ${#tool_names[@]}; ++i)); do
    tool_binding "${tool_paths[$i]}" "${tool_hashes[$i]}"
done

directory_fields() {
    local path=$1 expected_mode=$2 expected_uid=$3 expected_gid=$4
    local line kind dev ino uid gid mode nlink
    line=$("$stat_tool" -c '%F|%d|%i|%u|%g|%a|%h' -- "$path") || \
        die 'directory stat failed'
    IFS='|' read -r kind dev ino uid gid mode nlink <<<"$line"
    [[ $kind == directory && $uid == "$expected_uid" && \
       $gid == "$expected_gid" && \
       $mode == "$expected_mode" && $nlink -ge 2 ]] || \
        die 'directory binding changed'
    printf '%s|%s|%s|%s|%s|%s' "$dev" "$ino" "$uid" "$gid" \
        "$((8#$mode))" "$nlink"
}
anchor_before_fields=$(directory_fields "$anchor" 1777 0 0)
IFS='|' read -r abdev abino abuid abgid abmode _abnlink \
    <<<"$anchor_before_fields"
[[ ! -e $parent && ! -L $parent && ! -e $root && ! -L $root && \
   ! -e $receipt && ! -L $receipt ]] || \
    die 'stage residue requires forensic review and forbids retry'

source_bindings=()
for ((i = 0; i < 4; ++i)); do
    fd="${fds[$i]}"
    [[ $fd =~ ^[0-9]+$ && $fd -ge 3 ]] || die 'source FD is malformed'
    fd_path=/proc/self/fd/$fd
    fd_line=$("$stat_tool" -L -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "$fd_path") || die 'source FD stat failed'
    path_line=$("$stat_tool" -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "${sources[$i]}") || die 'source path stat failed'
    [[ $fd_line == "$path_line" ]] || die 'source path does not bind held FD'
    IFS='|' read -r kind dev ino size uid gid mode nlink <<<"$fd_line"
    [[ $kind == 'regular file' && $nlink == 1 ]] || \
        die 'source is not exact regular nlink1'
    hash_line=$("$sha_tool" -- "$fd_path") || die 'source FD hash failed'
    observed_hash=${hash_line%% *}
    [[ $observed_hash == "${expected_hashes[$i]}" ]] || \
        die 'source FD hash changed'
    source_bindings+=("$dev|$ino|$size|$uid|$gid|$((8#$mode))|$nlink")
done

"$mkdir_tool" --mode=0700 -- "$parent" || \
    die 'parent creation failed; residue forbids retry'
directory_fields "$parent" 700 "$owner_uid" "$owner_gid" >/dev/null
"$mkdir_tool" --mode=0700 -- "$root" || \
    die 'root creation failed; residue forbids retry'
directory_fields "$root" 700 "$owner_uid" "$owner_gid" >/dev/null
anchor_fields=$(directory_fields "$anchor" 1777 0 0)
IFS='|' read -r adev aino auid agid amode _anlink <<<"$anchor_fields"
[[ $adev == "$abdev" && $aino == "$abino" && $auid == "$abuid" && \
   $agid == "$abgid" && $amode == "$abmode" ]] || \
    die 'shared anchor identity changed during creation'
shopt -s nullglob dotglob
entries=("$root"/*)
[[ ${#entries[@]} -eq 0 ]] || die 'new stage root is not empty'

file_records=()
source_records=()
destination_bindings=()
for ((i = 0; i < 4; ++i)); do
    fd_path="/proc/self/fd/${fds[$i]}"
    destination="${destinations[$i]}"
    [[ ${destination%/*} == "$root" ]] || die 'destination escaped root'
    [[ ${destination##*/} == \
       $(case "${names[$i]}" in
           candidate) printf wirehair_expo_thermal_sampler.py ;;
           controller) printf wh2_thermal_sampler_transition.py ;;
           legacy) printf legacy_wirehair_expo_thermal_sampler.py ;;
           p32) printf wh2_p32_dispatch_timing.py ;;
       esac) ]] || die 'destination roster changed'
    [[ ! -e $destination && ! -L $destination ]] || \
        die 'destination collision; residue forbids retry'
    "$install_tool" -o "$owner_uid" -g "$owner_gid" -m 0444 -T -- \
        "$fd_path" "$destination" || \
        die 'install failed; residue forbids retry'
    line=$("$stat_tool" -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "$destination") || die 'destination stat failed'
    IFS='|' read -r kind dev ino size uid gid mode nlink <<<"$line"
    [[ $kind == 'regular file' && $uid == "$owner_uid" && \
       $gid == "$owner_gid" && \
       $mode == 444 && $nlink == 1 ]] || die 'destination binding changed'
    hash_line=$("$sha_tool" -- "$destination") || \
        die 'destination hash failed'
    observed_hash=${hash_line%% *}
    [[ $observed_hash == "${expected_hashes[$i]}" ]] || \
        die 'destination hash changed'
    IFS='|' read -r sdev sino ssize suid sgid smode snlink \
        <<<"${source_bindings[$i]}"
    [[ $size == "$ssize" ]] || die 'destination size changed'
    printf -v "file_records[$i]" \
        '"%s":{"binding":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"sha256":"%s","size":%s,"uid":%s},"path":"%s"}' \
        "${names[$i]}" "$dev" "$gid" "$ino" "$((8#$mode))" "$nlink" \
        "$observed_hash" "$size" "$uid" "$destination"
    printf -v "source_records[$i]" \
        '"%s":{"binding":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"sha256":"%s","size":%s,"uid":%s},"fd":%s,"path":"%s","stability_observations":2}' \
        "${names[$i]}" "$sdev" "$sgid" "$sino" "$smode" "$snlink" \
        "${expected_hashes[$i]}" "$ssize" "$suid" "${fds[$i]}" \
        "${sources[$i]}"
    destination_bindings+=(
        "$dev|$ino|$size|$uid|$gid|$((8#$mode))|$nlink")
done
entries=("$root"/*)
[[ ${#entries[@]} -eq 4 ]] || die 'code root roster changed'
for destination in "${destinations[@]}"; do
    [[ -f $destination && ! -L $destination ]] || die 'code file vanished'
done

# Re-observe every held FD and pathname after all copies.  Loaded bytes derive
# from the held FDs even if an attacker can swap and restore a mutable path
# between these discrete observations.
for ((i = 0; i < 4; ++i)); do
    fd_path="/proc/self/fd/${fds[$i]}"
    fd_line=$("$stat_tool" -L -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "$fd_path") || die 'post-copy source FD stat failed'
    path_line=$("$stat_tool" -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "${sources[$i]}") || die 'post-copy source path stat failed'
    [[ $fd_line == "$path_line" ]] || die 'source path changed after copy'
    IFS='|' read -r kind dev ino size uid gid mode nlink <<<"$fd_line"
    IFS='|' read -r sdev sino ssize suid sgid smode snlink \
        <<<"${source_bindings[$i]}"
    [[ $kind == 'regular file' && $dev == "$sdev" && $ino == "$sino" && \
       $size == "$ssize" && $uid == "$suid" && $gid == "$sgid" && \
       $((8#$mode)) == "$smode" && $nlink == "$snlink" ]] || \
        die 'source FD binding changed after copy'
    hash_line=$("$sha_tool" -- "$fd_path") || die 'post-copy FD hash failed'
    [[ ${hash_line%% *} == "${expected_hashes[$i]}" ]] || \
        die 'source FD hash changed after copy'
done

"$chmod_tool" 0555 -- "$root" || die 'root seal chmod failed'
root_fields=$(directory_fields "$root" 555 "$owner_uid" "$owner_gid")
parent_fields=$(directory_fields "$parent" 700 "$owner_uid" "$owner_gid")
IFS='|' read -r adev aino auid agid amode anlink <<<"$anchor_fields"
IFS='|' read -r pdev pino puid pgid _pmode pnlink <<<"$parent_fields"
IFS='|' read -r rdev rino ruid rgid rmode rnlink <<<"$root_fields"
files_json=$(IFS=,; printf '%s' "${file_records[*]}")
sources_json=$(IFS=,; printf '%s' "${source_records[*]}")
printf -v unhashed \
    '{"authority":{"perl_opener_sha256":"%s","stage_program_sha256":"%s"},"directories":{"anchor":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"uid":%s},"parent":{"device":%s,"gid":%s,"inode":%s,"mode":365,"nlink":%s,"uid":%s},"root":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"uid":%s}},"files":{%s},"no_live_state_or_workload":true,"partial_or_residual_stage_policy":"hard-stop-no-blind-retry","schema":"wirehair.wh2.thermal_transition.root_code_stage_receipt.v6","sources":{%s},"transition_id":"%s"}' \
    "$perl_program_sha256" "$bash_program_sha256" \
    "$adev" "$agid" "$aino" "$amode" "$anlink" "$auid" \
    "$pdev" "$pgid" "$pino" "$pnlink" "$puid" \
    "$rdev" "$rgid" "$rino" "$rmode" "$rnlink" "$ruid" \
    "$files_json" "$sources_json" "$transition_id"
hash_line=$(printf '%s\n' "$unhashed" | "$sha_tool") || \
    die 'receipt self hash failed'
self_hash=${hash_line%% *}
printf -v final \
    '{"authority":{"perl_opener_sha256":"%s","stage_program_sha256":"%s"},"directories":{"anchor":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"uid":%s},"parent":{"device":%s,"gid":%s,"inode":%s,"mode":365,"nlink":%s,"uid":%s},"root":{"device":%s,"gid":%s,"inode":%s,"mode":%s,"nlink":%s,"uid":%s}},"files":{%s},"no_live_state_or_workload":true,"partial_or_residual_stage_policy":"hard-stop-no-blind-retry","schema":"wirehair.wh2.thermal_transition.root_code_stage_receipt.v6","self_sha256_excluding_field":"%s","sources":{%s},"transition_id":"%s"}' \
    "$perl_program_sha256" "$bash_program_sha256" \
    "$adev" "$agid" "$aino" "$amode" "$anlink" "$auid" \
    "$pdev" "$pgid" "$pino" "$pnlink" "$puid" \
    "$rdev" "$rgid" "$rino" "$rmode" "$rnlink" "$ruid" \
    "$files_json" "$self_hash" "$sources_json" "$transition_id"
(set -o noclobber; printf '%s\n' "$final" >"$receipt") || \
    die 'receipt creation failed; residue forbids retry'
"$chmod_tool" 0444 -- "$receipt" || die 'receipt chmod failed'
line=$("$stat_tool" -c '%F|%u|%g|%a|%h' -- "$receipt") || \
    die 'receipt stat failed'
[[ $line == "regular file|$owner_uid|$owner_gid|444|1" ]] || \
    die 'receipt binding changed'
"$chmod_tool" 0555 -- "$parent" || die 'parent seal chmod failed'
final_parent_fields=$(directory_fields \
    "$parent" 555 "$owner_uid" "$owner_gid")
[[ $final_parent_fields == \
   "$pdev|$pino|$puid|$pgid|365|$pnlink" ]] || \
    die 'parent binding changed while sealing'
final_anchor_fields=$(directory_fields "$anchor" 1777 0 0)
IFS='|' read -r fadev faino fauid fagid famode _fanlink \
    <<<"$final_anchor_fields"
[[ $(directory_fields "$root" 555 "$owner_uid" "$owner_gid") == \
   "$root_fields" && \
   $fadev == "$adev" && $faino == "$aino" && $fauid == "$auid" && \
   $fagid == "$agid" && $famode == "$amode" ]] || \
    die 'sealed ancestry changed during final replay'
entries=("$parent"/*)
[[ ${#entries[@]} -eq 2 && -d $root && -f $receipt && ! -L $receipt ]] || \
    die 'parent roster changed'
for entry in "${entries[@]}"; do
    [[ $entry == "$root" || $entry == "$receipt" ]] || \
        die 'parent roster gained an unknown entry'
done
for ((i = 0; i < 4; ++i)); do
    line=$("$stat_tool" -c '%F|%d|%i|%s|%u|%g|%a|%h' -- \
        "${destinations[$i]}") || die 'final destination stat failed'
    IFS='|' read -r kind dev ino size uid gid mode nlink <<<"$line"
    IFS='|' read -r ddev dino dsize duid dgid dmode dnlink \
        <<<"${destination_bindings[$i]}"
    [[ $kind == 'regular file' && $dev == "$ddev" && $ino == "$dino" && \
       $size == "$dsize" && $uid == "$duid" && $gid == "$dgid" && \
       $((8#$mode)) == "$dmode" && $nlink == "$dnlink" ]] || \
        die 'final destination binding changed'
    hash_line=$("$sha_tool" -- "${destinations[$i]}") || \
        die 'final destination replay failed'
    [[ ${hash_line%% *} == "${expected_hashes[$i]}" ]] || \
        die 'final destination hash changed'
done
expected_receipt_hash_line=$(printf '%s\n' "$final" | "$sha_tool") || \
    die 'final receipt expectation hash failed'
actual_receipt_hash_line=$("$sha_tool" -- "$receipt") || \
    die 'final receipt replay failed'
[[ ${actual_receipt_hash_line%% *} == \
   "${expected_receipt_hash_line%% *}" ]] || die 'final receipt hash changed'
printf '%s\n' "$final"
'''


class TransitionError(RuntimeError):
    """A substantive transition or evidence failure."""


class TransitionDeadline(TransitionError):
    """The bounded transition exhausted its normal or recovery time."""


class ProcessSnapshotRace(TransitionError):
    """A proc record entered teardown after an exact stat observation."""


def renameat2_function() -> Callable[..., int]:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise TransitionError("renameat2 is unavailable")
    renameat2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint)
    renameat2.restype = ctypes.c_int
    return renameat2


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_fd(fd: int) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    os.lseek(fd, 0, os.SEEK_SET)
    while True:
        block = os.read(fd, 1024 * 1024)
        if not block:
            break
        digest.update(block)
        size += len(block)
    os.lseek(fd, 0, os.SEEK_SET)
    return digest.hexdigest(), size


def stable_file_bytes(path: Path, attempts: int = 20) -> bytes:
    for _attempt in range(attempts):
        before = path.stat()
        raw = path.read_bytes()
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size,
                before.st_mtime_ns) == (
                after.st_dev, after.st_ino, after.st_size,
                after.st_mtime_ns) and len(raw) == after.st_size:
            return raw
        time.sleep(0.02)
    raise TransitionError("file did not become stable: %s" % path)


def sha256_file(path: Path) -> str:
    return sha256_bytes(stable_file_bytes(path))


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if not schema or "schema" in payload or \
            "self_sha256_excluding_field" in payload:
        raise TransitionError("invalid sealed record payload")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        raise TransitionError(name + " schema mismatch")
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA256_RE.fullmatch(claimed) is None:
        raise TransitionError(name + " self hash malformed")
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        raise TransitionError(name + " self hash mismatch")
    return value


def load_canonical(path: Path, schema: str, name: str) -> Dict[str, object]:
    raw = stable_file_bytes(path)
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TransitionError(name + " is not canonical JSON") from exc
    if not isinstance(value, dict) or canonical_json(value) != raw:
        raise TransitionError(name + " is not canonical JSON")
    return verify_sealed(value, schema, name)


def fsync_directory(path: Path) -> None:
    fd = os.open(
        str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def file_binding(path: Path, *, with_hash: bool) -> Dict[str, object]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | \
        getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(path), flags)
    try:
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode):
            raise TransitionError("bound path is not a regular file: %s" % path)
        value: Dict[str, object] = {
            "device": info.st_dev, "gid": info.st_gid,
            "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
            "nlink": info.st_nlink, "size": info.st_size,
            "uid": info.st_uid,
        }
        if with_hash:
            digest, size = sha256_fd(fd)
            if size != info.st_size:
                raise TransitionError("bound file size changed while hashing")
            value["sha256"] = digest
        post = os.stat(str(path), follow_symlinks=False)
        if (post.st_dev, post.st_ino, post.st_size,
                stat.S_IMODE(post.st_mode), post.st_uid, post.st_gid,
                post.st_nlink) != (
                info.st_dev, info.st_ino, info.st_size,
                stat.S_IMODE(info.st_mode), info.st_uid, info.st_gid,
                info.st_nlink):
            raise TransitionError("bound path changed while reading: %s" % path)
        return value
    finally:
        os.close(fd)


def capture_tool_records() -> Dict[str, Dict[str, object]]:
    """Bind every executable in the retained historical tool ledger."""
    records: Dict[str, Dict[str, object]] = {}
    for name, path_text, expected_sha256 in TOOL_CONTRACT:
        if name in records or SHA256_RE.fullmatch(expected_sha256) is None:
            raise TransitionError("tool contract is malformed")
        path = Path(path_text)
        if not path.is_absolute():
            raise TransitionError("tool path is not absolute: " + name)
        binding = file_binding(path, with_hash=True)
        if binding["sha256"] != expected_sha256 or binding["uid"] != 0 or \
                binding["gid"] != 0 or binding["mode"] & 0o022 or \
                not binding["mode"] & 0o111 or binding["nlink"] < 1:
            raise TransitionError("hardcoded tool identity changed: " + name)
        records[name] = {
            "binding": binding,
            "path": path_text,
            "sha256": expected_sha256,
        }
    if set(records) != {
            "bash", "chmod", "env", "fuser", "install", "mkdir", "mpstat",
            "perl", "python3", "sha256sum", "stat", "sudo", "taskset",
            "timeout"}:
        raise TransitionError("tool contract is incomplete")
    return records


def verify_tool_records(value: object) -> Dict[str, Dict[str, object]]:
    """Rehash and replay the exact preparation-time executable bindings."""
    if not isinstance(value, dict):
        raise TransitionError("sealed tool ledger is malformed")
    current = capture_tool_records()
    if set(value) != set(current):
        raise TransitionError("sealed tool ledger changed")
    for name, record in current.items():
        if value.get(name) != record:
            raise TransitionError("sealed tool identity changed: " + name)
    return current


def verify_running_interpreter(
        tool_records: Mapping[str, Mapping[str, object]], *,
        require_exact_path: bool = False) -> Dict[str, object]:
    """Prove the executing controller is the sealed Python inode."""
    try:
        expected = tool_records["python3"]["binding"]
        if not isinstance(expected, dict):
            raise TypeError("Python binding is not a mapping")
        info = os.stat("/proc/self/exe", follow_symlinks=True)
        observed = {"device": info.st_dev, "inode": info.st_ino}
        expected_path = str(tool_records["python3"]["path"])
        resolved_path = str(Path(sys.executable).resolve(strict=True))
        if observed != {
                "device": expected["device"], "inode": expected["inode"]} or \
                resolved_path != expected_path or \
                (require_exact_path and sys.executable != expected_path):
            raise TransitionError(
                "controller interpreter differs from the sealed Python tool")
        return {**observed, "argv_path": sys.executable,
                "resolved_path": resolved_path}
    except TransitionError:
        raise
    except (KeyError, OSError, RuntimeError, TypeError) as exc:
        raise TransitionError(
            "cannot bind the controller interpreter to the sealed Python tool") \
            from exc


def _directory_binding_fd(fd: int) -> Dict[str, int]:
    info = os.fstat(fd)
    if not stat.S_ISDIR(info.st_mode):
        raise TransitionError("bound parent is not a directory")
    return {
        "device": info.st_dev, "gid": info.st_gid,
        "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink, "uid": info.st_uid,
    }


def directory_binding(path: Path) -> Dict[str, int]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | os.O_DIRECTORY | \
        getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(path), flags)
    try:
        observed = _directory_binding_fd(fd)
        post = os.stat(str(path), follow_symlinks=False)
        if observed != {
                "device": post.st_dev, "gid": post.st_gid,
                "inode": post.st_ino, "mode": stat.S_IMODE(post.st_mode),
                "nlink": post.st_nlink, "uid": post.st_uid}:
            raise TransitionError("bound parent path changed")
        return observed
    finally:
        os.close(fd)


def write_new(path: Path, raw: bytes, *, owner_uid: int,
              mode: int = 0o444) -> Dict[str, object]:
    """Create, fsync, seal, and bind one controller-owned record."""
    parent_binding = directory_binding(path.parent)
    if parent_binding["uid"] != owner_uid or parent_binding["mode"] & 0o022:
        raise TransitionError("record parent is outside the controller trust boundary")
    dir_fd = os.open(
        str(path.parent), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    fd = -1
    try:
        if _directory_binding_fd(dir_fd) != parent_binding:
            raise TransitionError("record parent changed before creation")
        fd = os.open(
            path.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0), mode, dir_fd=dir_fd)
        view = memoryview(raw)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                raise OSError(errno.EIO, "short record write")
            view = view[written:]
        os.fsync(fd)
        os.fchmod(fd, mode)
        os.fsync(fd)
        info = os.fstat(fd)
        if info.st_uid != owner_uid or info.st_nlink != 1 or \
                stat.S_IMODE(info.st_mode) != mode:
            raise TransitionError("new record binding is unsafe")
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent_binding:
            raise TransitionError("record parent changed during creation")
    finally:
        if fd >= 0:
            os.close(fd)
        os.close(dir_fd)
    result = file_binding(path, with_hash=True)
    if directory_binding(path.parent) != parent_binding:
        raise TransitionError("record parent path changed after creation")
    if result["sha256"] != sha256_bytes(raw) or result["size"] != len(raw):
        raise TransitionError("new record content changed")
    return result


def copy_sealed(source: Path, destination: Path, expected_sha256: str,
                *, owner_uid: int) -> Dict[str, object]:
    if SHA256_RE.fullmatch(expected_sha256) is None:
        raise TransitionError("expected copy SHA256 is malformed")
    source_fd = os.open(
        str(source), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK |
        getattr(os, "O_NOFOLLOW", 0))
    destination_fd = -1
    dir_fd = -1
    try:
        source_info = os.fstat(source_fd)
        if not stat.S_ISREG(source_info.st_mode) or source_info.st_nlink != 1:
            raise TransitionError("copy source is not a single-link regular file")
        digest, _size = sha256_fd(source_fd)
        if digest != expected_sha256:
            raise TransitionError("copy source SHA256 changed")
        parent = directory_binding(destination.parent)
        if parent["uid"] != owner_uid or parent["mode"] & 0o022:
            raise TransitionError("copy destination parent is unsafe")
        dir_fd = os.open(
            str(destination.parent), os.O_RDONLY | os.O_CLOEXEC |
            os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0))
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("copy destination parent changed")
        destination_fd = os.open(
            destination.name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0), 0o444, dir_fd=dir_fd)
        os.lseek(source_fd, 0, os.SEEK_SET)
        while True:
            block = os.read(source_fd, 1024 * 1024)
            if not block:
                break
            view = memoryview(block)
            while view:
                written = os.write(destination_fd, view)
                if written <= 0:
                    raise OSError(errno.EIO, "short sealed-copy write")
                view = view[written:]
        os.fsync(destination_fd)
        os.fchmod(destination_fd, 0o444)
        os.fsync(destination_fd)
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("copy destination parent changed during copy")
    finally:
        if destination_fd >= 0:
            os.close(destination_fd)
        if dir_fd >= 0:
            os.close(dir_fd)
        os.close(source_fd)
    binding = file_binding(destination, with_hash=True)
    if directory_binding(destination.parent) != parent:
        raise TransitionError("copy destination parent path changed")
    if binding["sha256"] != expected_sha256 or binding["uid"] != owner_uid or \
            binding["mode"] != 0o444 or binding["nlink"] != 1:
        raise TransitionError("sealed copy binding changed")
    return binding


def rename_noreplace(source: Path, destination: Path,
                     expected_binding: Mapping[str, object], *,
                     parent_uid: int) -> Dict[str, object]:
    """Rename one exact inode with renameat2(RENAME_NOREPLACE)."""
    if source.parent != destination.parent or source.name in ("", ".", "..") or \
            destination.name in ("", ".", ".."):
        raise TransitionError("archive rename must stay within one parent")
    parent = directory_binding(source.parent)
    if parent["uid"] != parent_uid or parent["mode"] & 0o022:
        raise TransitionError("archive parent is outside its trust boundary")
    observed = file_binding(source, with_hash="sha256" in expected_binding)
    if observed != dict(expected_binding):
        raise TransitionError("archive source binding changed")
    dir_fd = os.open(
        str(source.parent), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    try:
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("archive parent changed before rename")
        try:
            os.stat(destination.name, dir_fd=dir_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise FileExistsError(errno.EEXIST, "archive destination exists",
                                  str(destination))
        renameat2 = renameat2_function()
        result = renameat2(
            dir_fd, os.fsencode(source.name), dir_fd,
            os.fsencode(destination.name), RENAME_NOREPLACE)
        if result != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), str(destination))
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("archive parent changed during rename")
        try:
            os.stat(source.name, dir_fd=dir_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise TransitionError("archive source name survived rename")
    finally:
        os.close(dir_fd)
    archived = file_binding(destination, with_hash="sha256" in expected_binding)
    if archived != dict(expected_binding):
        raise TransitionError("archive destination binding changed")
    if directory_binding(source.parent) != parent:
        raise TransitionError("archive parent binding changed")
    return archived


@dataclass(frozen=True)
class TransitionPlan:
    transition_id: str = TRANSITION_ID
    root: Path = DRY_ROOT
    controller_uid: int = 1000
    controller_gid: int = 1000
    code_owner_uid: int = 0
    code_owner_gid: int = 0
    code_anchor: Path = Path("/dev/shm")
    code_seal_parent: Path = CODE_SEAL_PARENT
    code_seal_root: Path = CODE_SEAL_ROOT
    deadline_s: float = 540.0
    recovery_reserve_s: float = 60.0
    emergency_recovery_s: float = 60.0
    # Supplied by an exact-byte fixture for deterministic source-test plan
    # construction and historical contract replay.
    controller_sha256: str = ""
    candidate_sha256: str = CANDIDATE_SHA256
    p32_sha256: str = P32_SHA256
    old_source_sha256: str = OLD_SOURCE_SHA256
    old_boot_id: str = OLD_BOOT_ID
    old_pid: int = 3320493
    old_start_tick: int = 160912119
    old_cmdline_sha256: str = OLD_CMDLINE_SHA256
    old_launcher_pid: int = 3320490
    old_launcher_start_tick: int = 160912119
    old_launcher_cmdline_sha256: str = OLD_LAUNCHER_CMDLINE_SHA256
    old_launcher_affinity: str = "16-59,81-123"
    old_launcher_uids: Tuple[int, int, int, int] = (1000, 0, 0, 0)
    old_process_group: int = 3320490
    old_session: int = 3320465
    old_csv_device: int = 28
    old_csv_inode: int = 10141
    old_pid_device: int = 28
    old_pid_inode: int = 10140
    old_cpu: int = 127
    candidate_cpu: int = 124
    candidate_sibling: int = 60
    stage_candidate_source: Path = REPO_CANDIDATE
    stage_controller_source: Path = REPO_CONTROLLER
    stage_legacy_source: Path = OLD_SOURCE
    stage_p32_source: Path = REPO_P32
    old_source: Path = OLD_SOURCE
    old_csv: Path = OLD_CSV
    old_pid_file: Path = OLD_PID_FILE
    old_archive: Path = OLD_ARCHIVE
    old_unclean_archive: Path = OLD_UNCLEAN_ARCHIVE
    old_stale_pid_archive: Path = OLD_STALE_PID_ARCHIVE
    audit_binding: Path = AUDIT_BINDING
    old_argv: Tuple[str, ...] = OLD_ARGV
    old_launch_argv: Tuple[str, ...] = OLD_LAUNCH_ARGV
    replacement_old_argv: Tuple[str, ...] = OLD_REPLACEMENT_ARGV
    replacement_old_cmdline_sha256: str = OLD_REPLACEMENT_CMDLINE_SHA256

    @property
    def receipts(self) -> Path:
        return self.root / "receipts"

    @property
    def code_stage_receipt(self) -> Path:
        return self.root / "root-code-seal-stage-receipt.v6.json"

    @property
    def root_code_stage_receipt(self) -> Path:
        return self.code_seal_parent / ROOT_CODE_STAGE_RECEIPT.name

    @property
    def sampler(self) -> Path:
        return self.code_seal_root / "wirehair_expo_thermal_sampler.py"

    @property
    def p32(self) -> Path:
        return self.code_seal_root / "wh2_p32_dispatch_timing.py"

    @property
    def controller(self) -> Path:
        return self.code_seal_root / "wh2_thermal_sampler_transition.py"

    @property
    def legacy_sampler(self) -> Path:
        return self.code_seal_root / "legacy_wirehair_expo_thermal_sampler.py"

    @property
    def plan_receipt(self) -> Path:
        return self.root / "transition-plan.json"


def _require_source_test_namespace(
        plan: TransitionPlan, operation: str) -> None:
    """Reject production authority from retained synthetic lower models."""
    production_roots = (DRY_ROOT, CODE_SEAL_PARENT)
    observed_paths = {
        plan.root, plan.code_seal_parent, plan.code_seal_root,
        plan.root_code_stage_receipt,
    }
    def within(path: Path, root: Path) -> bool:
        normalized = Path(os.path.abspath(os.fspath(path)))
        normalized_root = Path(os.path.abspath(os.fspath(root)))
        return normalized == normalized_root or \
            normalized_root in normalized.parents

    if re.fullmatch(
            re.escape(SOURCE_TEST_TRANSITION_PREFIX) + r"[A-Za-z0-9._-]+",
            plan.transition_id) is None or \
            plan.transition_id == TRANSITION_ID or \
            any(within(path, root) for path in observed_paths
                for root in production_roots) or \
            plan.controller_uid != os.geteuid() or \
            plan.controller_gid != os.getegid() or \
            plan.code_owner_uid != os.geteuid() or \
            plan.code_owner_gid != os.getegid():
        raise TransitionError(
            operation + " rejects the production V6 namespace")


class SplitCustodyRoot:
    """Route frozen-code lookups away from the writable evidence root.

    The reviewed P32 helper is byte-frozen and expresses its sampler lookup as
    ``root / "frozen/..."``.  This narrow path adapter preserves that API while
    routing only the frozen namespace to the separately root-owned code seal.
    Every other child is an ordinary ``Path`` beneath the UID-owned output
    root; returned children never retain implicit routing state.
    """

    def __init__(self, output_root: Path, code_seal_root: Path) -> None:
        self.output_root = Path(output_root)
        self.code_seal_root = Path(code_seal_root)

    def __fspath__(self) -> str:
        return os.fspath(self.output_root)

    def __str__(self) -> str:
        return str(self.output_root)

    def __truediv__(self, child: object) -> Path:
        text = os.fspath(child)
        relative = Path(text)
        if not text or relative.is_absolute() or ".." in relative.parts:
            raise TransitionError("split-custody child path is malformed")
        if text == "frozen":
            return self.code_seal_root
        prefix = "frozen/"
        if text.startswith(prefix):
            suffix = text[len(prefix):]
            if not suffix or suffix.startswith("/") or ".." in Path(suffix).parts:
                raise TransitionError("split-custody frozen path is malformed")
            return self.code_seal_root / suffix
        return self.output_root / text

    def resolve(self, strict: bool = False) -> Path:
        return self.output_root.resolve(strict=strict)


def _exact_directory(path: Path, *, uid: int, gid: int, mode: int,
                     description: str) -> Dict[str, int]:
    binding = directory_binding(path)
    if binding["uid"] != uid or binding["gid"] != gid or \
            binding["mode"] != mode or binding["nlink"] < 2:
        raise TransitionError(description + " binding is unsafe")
    return binding


def _exact_directory_roster(path: Path, expected: Sequence[str],
                            description: str) -> None:
    try:
        observed = sorted(entry.name for entry in path.iterdir())
    except OSError as exc:
        raise TransitionError(description + " roster is unreadable") from exc
    if observed != sorted(expected):
        raise TransitionError(description + " roster is incomplete or enlarged")


def _code_file_contract(plan: TransitionPlan) -> Tuple[Tuple[str, Path, str], ...]:
    contract = (
        ("candidate", plan.sampler, plan.candidate_sha256),
        ("controller", plan.controller, plan.controller_sha256),
        ("legacy", plan.legacy_sampler, plan.old_source_sha256),
        ("p32", plan.p32, plan.p32_sha256),
    )
    if len({path.name for _name, path, _digest in contract}) != len(contract) or \
            any(SHA256_RE.fullmatch(digest) is None
                for _name, _path, digest in contract):
        raise TransitionError("root code-seal contract is malformed")
    return contract


def _verify_code_anchor(plan: TransitionPlan) -> Dict[str, int]:
    if plan.code_seal_parent.parent != plan.code_anchor or \
            plan.code_seal_root.parent != plan.code_seal_parent or \
            plan.code_seal_root.name != "seal-v6":
        raise TransitionError("root code-seal ancestry changed")
    anchor = directory_binding(plan.code_anchor)
    if anchor["nlink"] < 2:
        raise TransitionError("root code-seal mount anchor is unsafe")
    if plan.code_anchor == Path("/dev/shm"):
        if anchor["uid"] != 0 or anchor["gid"] != 0 or \
                anchor["mode"] != 0o1777:
            raise TransitionError(
                "root code-seal /dev/shm anchor binding changed")
    elif anchor["uid"] != plan.code_owner_uid or \
            anchor["gid"] != plan.code_owner_gid or anchor["mode"] & 0o022:
        raise TransitionError("root code-seal mount anchor is writable")
    return anchor


def _root_stage_program_hashes() -> Dict[str, str]:
    return {
        "perl_opener_sha256": sha256_bytes(
            ROOT_STAGE_PERL_OPENER.encode("ascii")),
        "stage_program_sha256": sha256_bytes(
            ROOT_STAGE_BASH_PROGRAM.encode("ascii")),
    }


def _root_stage_tool_contract() -> Dict[str, Dict[str, str]]:
    by_name = {
        name: {"path": path, "sha256": digest}
        for name, path, digest in TOOL_CONTRACT
        if name in ROOT_STAGE_TOOL_NAMES
    }
    if tuple(sorted(by_name)) != ROOT_STAGE_TOOL_NAMES:
        raise TransitionError("root-stage tool contract is malformed")
    return by_name


def _literal_stage_path(path: Path, description: str) -> str:
    text = str(path)
    if not path.is_absolute() or ".." in path.parts or \
            re.fullmatch(r"/[A-Za-z0-9._/-]+", text) is None:
        raise TransitionError(description + " is not an exact safe path")
    return text


def _root_stage_source_contract(
        plan: TransitionPlan, *, controller_source: Path,
        candidate_source: Path, p32_source: Path,
        legacy_source: Optional[Path],
) -> Tuple[Tuple[str, Path, Path, str], ...]:
    legacy = plan.stage_legacy_source if legacy_source is None else \
        Path(legacy_source)
    sources = {
        "candidate": Path(candidate_source),
        "controller": Path(controller_source),
        "legacy": legacy,
        "p32": Path(p32_source),
    }
    expected_sources = {
        "candidate": plan.stage_candidate_source,
        "controller": plan.stage_controller_source,
        "legacy": plan.stage_legacy_source,
        "p32": plan.stage_p32_source,
    }
    if sources != expected_sources:
        raise TransitionError("root-stage source path contract changed")
    destinations = {
        name: path for name, path, _digest in _code_file_contract(plan)}
    expected = {
        name: digest for name, _path, digest in _code_file_contract(plan)}
    return tuple(
        (name, sources[name], destinations[name], expected[name])
        for name in ("candidate", "controller", "legacy", "p32"))


def _code_seal_stage_plan_value(
        plan: TransitionPlan, *, controller_source: Path,
        candidate_source: Path = REPO_CANDIDATE,
        p32_source: Path = REPO_P32,
        legacy_source: Optional[Path] = None,
) -> Dict[str, object]:
    """Build a literal authority for a synthetic source-test namespace.

    This function is intentionally pure with respect to source files, tools,
    subprocesses, and stage destinations.  It rejects the production V6
    namespace and exists only to exercise the historical authority grammar in
    isolated tests.  The mutable controller is never a privileged staging
    authority.
    """
    _require_source_test_namespace(plan, "root-stage source-test model")
    if type(plan.code_owner_uid) is not int or \
            type(plan.code_owner_gid) is not int or \
            not 0 <= plan.code_owner_uid <= 4294967294 or \
            not 0 <= plan.code_owner_gid <= 4294967294 or \
            plan.code_anchor != Path("/dev/shm") or \
            plan.code_seal_parent != \
                plan.code_anchor / (plan.transition_id + "-root-code") or \
            plan.code_seal_root != plan.code_seal_parent / "seal-v6" or \
            plan.root_code_stage_receipt != plan.code_seal_parent / \
                "root-code-seal-stage-receipt.v6.json":
        raise TransitionError("root-stage destination contract changed")
    contract = _root_stage_source_contract(
        plan, controller_source=Path(controller_source),
        candidate_source=Path(candidate_source), p32_source=Path(p32_source),
        legacy_source=legacy_source)
    for name, source, destination, digest in contract:
        _literal_stage_path(source, name + " stage source")
        _literal_stage_path(destination, name + " stage destination")
        if SHA256_RE.fullmatch(digest) is None:
            raise TransitionError(name + " stage digest is malformed")
    for path, description in (
            (plan.code_anchor, "stage anchor"),
            (plan.code_seal_parent, "stage parent"),
            (plan.code_seal_root, "stage root"),
            (plan.root_code_stage_receipt, "stage receipt")):
        _literal_stage_path(path, description)
    programs = _root_stage_program_hashes()
    tools = _root_stage_tool_contract()
    bash_argv = [
        tools["bash"]["path"], "--noprofile", "--norc", "-c",
        ROOT_STAGE_BASH_PROGRAM, "wirehair-root-code-stage",
        plan.transition_id, str(plan.code_anchor),
        str(plan.code_seal_parent), str(plan.code_seal_root),
        str(plan.root_code_stage_receipt),
        str(plan.code_owner_uid), str(plan.code_owner_gid),
    ]
    for _name, source, destination, digest in contract:
        bash_argv.extend((str(source), str(destination), digest))
    bash_argv.extend((
        programs["perl_opener_sha256"],
        programs["stage_program_sha256"]))
    authority_argv = [
        tools["sudo"]["path"], "-n", "--",
        tools["timeout"]["path"], "--signal=KILL",
        "--kill-after=1.000s", "30.000s",
        tools["env"]["path"], "-i", "HOME=/root", "LANG=C",
        "LC_ALL=C", "PATH=/usr/bin:/bin", "TZ=UTC",
        tools["perl"]["path"], "-e", ROOT_STAGE_PERL_OPENER,
        *(str(source) for _name, source, _destination, _digest in contract),
        *bash_argv,
    ]
    return sealed_record(CODE_SEAL_STAGE_SCHEMA, {
        "authority_argv": authority_argv,
        "authority_argv_canonical_sha256":
            sha256_bytes(canonical_json(authority_argv)),
        "authority_argv_nul_sha256": sha256_bytes(
            b"\0".join(value.encode("ascii") for value in authority_argv) +
            b"\0"),
        "destinations": {
            "parent": str(plan.code_seal_parent),
            "receipt": str(plan.root_code_stage_receipt),
            "root": str(plan.code_seal_root),
        },
        "no_live_state_or_workload": True,
        "partial_or_residual_stage_policy": "hard-stop-no-blind-retry",
        "programs": programs,
        "reviewer_authorization": (
            "synthetic-source-test-only; never production authorization"),
        "sources": {
            name: {"expected_sha256": digest, "path": str(source)}
            for name, source, _destination, digest in contract
        },
        "tools": tools,
        "transition_id": plan.transition_id,
    })


def code_seal_stage_plan(
        plan: TransitionPlan, *, controller_source: Path,
        candidate_source: Path = REPO_CANDIDATE,
        p32_source: Path = REPO_P32,
        legacy_source: Optional[Path] = None,
) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)


def _stable_anchor_binding(value: Mapping[str, object]) -> Dict[str, object]:
    """Return only the stable identity fields of the shared sticky anchor.

    Creating the separate UID-owned output root directly under ``/dev/shm``
    legitimately increments that directory's link count after the root-code
    stage.  Root-owned descendants remain exact, but shared-anchor nlink is an
    observation rather than immutable provenance.
    """
    keys = ("device", "gid", "inode", "mode", "uid")
    if any(type(value.get(key)) is not int or value[key] < 0 for key in keys):
        raise TransitionError("root code-seal anchor receipt is malformed")
    return {key: value[key] for key in keys}


def _verify_root_stage_receipt(
        plan: TransitionPlan, *, anchor: Mapping[str, object],
        parent: Mapping[str, object], root: Mapping[str, object],
        files: Mapping[str, object],
) -> Dict[str, object]:
    binding = file_binding(plan.root_code_stage_receipt, with_hash=True)
    if binding["uid"] != plan.code_owner_uid or \
            binding["gid"] != plan.code_owner_gid or \
            binding["mode"] != 0o444 or binding["nlink"] != 1:
        raise TransitionError("root code-stage receipt binding changed")
    receipt = load_canonical(
        plan.root_code_stage_receipt, ROOT_CODE_STAGE_RECEIPT_SCHEMA,
        "root code-stage receipt")
    directories = receipt.get("directories")
    receipt_anchor = directories.get("anchor") \
        if isinstance(directories, dict) else None
    if receipt.get("authority") != _root_stage_program_hashes() or \
            not isinstance(receipt_anchor, dict) or \
            _stable_anchor_binding(receipt_anchor) != \
                _stable_anchor_binding(anchor) or \
            receipt_anchor.get("nlink") is None or \
            type(receipt_anchor.get("nlink")) is not int or \
            receipt_anchor["nlink"] < 2 or \
            not isinstance(directories, dict) or \
            directories.get("parent") != dict(parent) or \
            directories.get("root") != dict(root) or \
            receipt.get("files") != dict(files) or \
            receipt.get("no_live_state_or_workload") is not True or \
            receipt.get("partial_or_residual_stage_policy") != \
                "hard-stop-no-blind-retry" or \
            receipt.get("transition_id") != plan.transition_id:
        raise TransitionError("root code-stage receipt contract changed")
    expected_sources = {
        "candidate": (plan.stage_candidate_source, plan.candidate_sha256),
        "controller": (plan.stage_controller_source, plan.controller_sha256),
        "legacy": (plan.stage_legacy_source, plan.old_source_sha256),
        "p32": (plan.stage_p32_source, plan.p32_sha256),
    }
    sources = receipt.get("sources")
    if not isinstance(sources, dict) or set(sources) != set(expected_sources):
        raise TransitionError("root code-stage source receipt roster changed")
    binding_keys = {
        "device", "gid", "inode", "mode", "nlink", "sha256", "size",
        "uid"}
    seen_fds = set()
    for name, (path, digest) in expected_sources.items():
        record = sources.get(name)
        source_binding = record.get("binding") \
            if isinstance(record, dict) else None
        fd = record.get("fd") if isinstance(record, dict) else None
        if not isinstance(source_binding, dict) or \
                set(source_binding) != binding_keys or \
                source_binding.get("sha256") != digest or \
                source_binding.get("nlink") != 1 or \
                any(type(source_binding.get(key)) is not int or
                    source_binding[key] < 0
                    for key in binding_keys - {"sha256"}) or \
                type(fd) is not int or fd < 3 or fd in seen_fds or \
                record.get("path") != str(path) or \
                record.get("stability_observations") != 2:
            raise TransitionError(
                "root code-stage source receipt changed: " + name)
        _literal_stage_path(
            Path(record["path"]), name + " root-stage receipt source")
        seen_fds.add(fd)
    return {"binding": binding, "path": str(plan.root_code_stage_receipt),
            "value": receipt}


def verify_code_seal(plan: TransitionPlan) -> Dict[str, object]:
    """Replay the separate root-owned, immutable executable-code seal."""
    anchor = _verify_code_anchor(plan)
    # /dev/shm is intentionally a sticky shared mount.  Every transition-owned
    # ancestor below that anchor is root-owned and non-writable.
    parent = _exact_directory(
        plan.code_seal_parent, uid=plan.code_owner_uid,
        gid=plan.code_owner_gid, mode=0o555,
        description="root code-seal parent")
    root = _exact_directory(
        plan.code_seal_root, uid=plan.code_owner_uid,
        gid=plan.code_owner_gid, mode=0o555,
        description="root code-seal directory")
    _exact_directory_roster(
        plan.code_seal_parent,
        (plan.code_seal_root.name, plan.root_code_stage_receipt.name),
        "root code-seal parent")
    contract = _code_file_contract(plan)
    _exact_directory_roster(
        plan.code_seal_root, tuple(path.name for _name, path, _digest in contract),
        "root code-seal")
    files: Dict[str, object] = {}
    for name, path, digest in contract:
        binding = file_binding(path, with_hash=True)
        if binding["sha256"] != digest or \
                binding["uid"] != plan.code_owner_uid or \
                binding["gid"] != plan.code_owner_gid or \
                binding["mode"] != 0o444 or binding["nlink"] != 1:
            raise TransitionError("root code-seal file changed: " + name)
        files[name] = {"binding": binding, "path": str(path)}
    stage_receipt = _verify_root_stage_receipt(
        plan, anchor=anchor, parent=parent, root=root, files=files)
    replay_anchor = _verify_code_anchor(plan)
    if directory_binding(plan.code_seal_parent) != parent or \
            directory_binding(plan.code_seal_root) != root or \
            _stable_anchor_binding(replay_anchor) != \
                _stable_anchor_binding(anchor):
        raise TransitionError("root code-seal ancestry changed during replay")
    return sealed_record(CODE_SEAL_RECEIPT_SCHEMA, {
        "anchor": {"binding": _stable_anchor_binding(anchor),
                   "path": str(plan.code_anchor),
                   "nlink_policy":
                       "volatile-shared-anchor-observed-minimum-two",
                   "shared_sticky_mount_exception":
                       plan.code_anchor == Path("/dev/shm")},
        "files": files,
        "parent": {"binding": parent, "path": str(plan.code_seal_parent)},
        "root": {"binding": root, "path": str(plan.code_seal_root)},
        "root_stage_receipt": stage_receipt,
        "transition_id": plan.transition_id,
    })


def execute_environment(plan: TransitionPlan) -> Dict[str, str]:
    """Historical exact environment, retained as a pure source/test fixture."""
    return {
        "HOME": str(plan.root / "runtime-home"),
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "TZ": "UTC",
    }


def prepare_environment(plan: TransitionPlan) -> Dict[str, str]:
    """Historical environment retained for synthetic preparation tests."""
    return execute_environment(plan)


def expected_prepare_orig_argv(plan: TransitionPlan) -> Tuple[str, ...]:
    if SHA256_RE.fullmatch(plan.controller_sha256) is None:
        raise TransitionError(
            "externally reviewed controller SHA256 is required")
    return (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--prepare-sealed-transition",
        "--expected-controller-sha256", plan.controller_sha256,
    )


def expected_prepare_command(plan: TransitionPlan) -> Tuple[str, ...]:
    environment = prepare_environment(plan)
    return (
        "/usr/bin/env", "-i",
        *("%s=%s" % (name, environment[name])
          for name in ("HOME", "PATH", "LANG", "LC_ALL", "TZ")),
        *expected_prepare_orig_argv(plan),
    )


def prepare_runtime_contract(
        plan: TransitionPlan,
        tool_records: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    try:
        python = tool_records["python3"]
        binding = python["binding"]
        path = python["path"]
        if not isinstance(binding, dict) or path != "/usr/bin/python3.12" or \
                type(binding.get("device")) is not int or \
                type(binding.get("inode")) is not int:
            raise TypeError("Python tool record is malformed")
    except (KeyError, TypeError) as exc:
        raise TransitionError(
            "prepare runtime lacks the sealed Python tool") from exc
    command = expected_prepare_command(plan)
    return {
        "command": list(command),
        "command_nul_sha256": sha256_bytes(
            b"\0".join(value.encode("ascii") for value in command) + b"\0"),
        "controller_interpreter": {
            "argv_path": path,
            "device": binding["device"],
            "inode": binding["inode"],
            "resolved_path": path,
        },
        "environment": prepare_environment(plan),
        "flags": dict(PREPARE_FLAG_CONTRACT),
        "sys_orig_argv": list(expected_prepare_orig_argv(plan)),
    }


def verify_prepare_runtime(
        plan: TransitionPlan,
        tool_records: Mapping[str, Mapping[str, object]], *,
        observed_environment: Optional[Mapping[str, str]] = None,
        observed_orig_argv: Optional[Sequence[str]] = None,
        observed_flags: Optional[Mapping[str, int | bool]] = None,
) -> Dict[str, object]:
    """Attest the isolated runtime of synthetic preparation tests."""
    environment = dict(os.environ if observed_environment is None
                       else observed_environment)
    orig_argv = tuple(sys.orig_argv if observed_orig_argv is None
                      else observed_orig_argv)
    flags = {
        name: getattr(sys.flags, name) for name in PREPARE_FLAG_CONTRACT
    } if observed_flags is None else dict(observed_flags)
    expected = prepare_runtime_contract(plan, tool_records)
    if environment != expected["environment"]:
        raise TransitionError(
            "prepare environment differs from exact env -i contract")
    if orig_argv != tuple(expected["sys_orig_argv"]):
        raise TransitionError(
            "prepare sys.orig_argv differs from frozen contract")
    if flags != PREPARE_FLAG_CONTRACT:
        raise TransitionError(
            "prepare sys.flags differ from isolated Python contract")
    interpreter = verify_running_interpreter(
        tool_records, require_exact_path=True)
    if interpreter != expected["controller_interpreter"]:
        raise TransitionError(
            "prepare interpreter differs from exact runtime contract")
    return expected


def expected_execute_orig_argv(plan: TransitionPlan) -> Tuple[str, ...]:
    """Build the retired V5-style execute argv without dispatching it."""
    if SHA256_RE.fullmatch(plan.controller_sha256) is None:
        raise TransitionError(
            "externally reviewed controller SHA256 is required")
    return (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--execute-sealed-transition", plan.transition_id,
        "--expected-controller-sha256", plan.controller_sha256,
        "--confirmation", EXECUTE_CONFIRMATION,
    )


def expected_execute_command(plan: TransitionPlan) -> Tuple[str, ...]:
    """Build the retired execute command as non-authoritative test data."""
    environment = execute_environment(plan)
    return (
        "/usr/bin/env", "-i",
        *("%s=%s" % (name, environment[name])
          for name in ("HOME", "PATH", "LANG", "LC_ALL", "TZ")),
        *expected_execute_orig_argv(plan),
    )


def python_isolation_contract(plan: TransitionPlan) -> Dict[str, object]:
    prefix = ["/usr/bin/python3.12", "-I", "-S", "-B"]
    return {
        "candidate_argv_prefix": prefix,
        "controller_orig_argv": list(expected_execute_orig_argv(plan)),
        "helper_argv_prefix": prefix,
        "replacement_old_argv": list(plan.replacement_old_argv),
    }


def verify_execute_runtime(
        plan: TransitionPlan, *,
        observed_environment: Optional[Mapping[str, str]] = None,
        observed_orig_argv: Optional[Sequence[str]] = None,
        observed_flags: Optional[Mapping[str, int | bool]] = None,
) -> Dict[str, object]:
    """Validate the retired runtime contract as pure historical evidence."""
    environment = dict(os.environ if observed_environment is None
                       else observed_environment)
    orig_argv = tuple(sys.orig_argv if observed_orig_argv is None
                      else observed_orig_argv)
    if observed_flags is None:
        flags = {
            name: getattr(sys.flags, name) for name in EXECUTE_FLAG_CONTRACT}
    else:
        flags = dict(observed_flags)
    expected_environment = execute_environment(plan)
    expected_orig_argv = expected_execute_orig_argv(plan)
    if environment != expected_environment:
        raise TransitionError("execute environment differs from exact env -i contract")
    if orig_argv != expected_orig_argv:
        raise TransitionError("execute sys.orig_argv differs from frozen contract")
    if flags != EXECUTE_FLAG_CONTRACT:
        raise TransitionError("execute sys.flags differ from isolated Python contract")
    return {
        "command": list(expected_execute_command(plan)),
        "environment": expected_environment,
        "flags": dict(EXECUTE_FLAG_CONTRACT),
        "sys_orig_argv": list(expected_orig_argv),
    }


IDENTITY_ABSENCE_POLICIES = {
    "none": frozenset(),
    "launcher": frozenset({"launcher"}),
    "child": frozenset({"child"}),
    "both": frozenset({"child", "launcher"}),
}


def replacement_launch_command(
        plan: TransitionPlan,
        tools: Mapping[str, Mapping[str, object]],
) -> Tuple[str, ...]:
    """Return the historical root-sealed legacy restart argv as test data."""
    try:
        sudo = tools["sudo"]["path"]
        env = tools["env"]["path"]
        taskset = tools["taskset"]["path"]
    except (KeyError, TypeError) as exc:
        raise TransitionError(
            "replacement launcher lacks its sealed tool paths") from exc
    environment = execute_environment(plan)
    return (
        str(sudo), "-n", "-b", str(env), "-i",
        *("%s=%s" % (name, environment[name])
          for name in ("HOME", "PATH", "LANG", "LC_ALL", "TZ")),
        str(taskset), "--cpu-list", str(plan.old_cpu),
        *plan.replacement_old_argv,
    )


def identity_inspection_environment(plan: TransitionPlan) -> Dict[str, str]:
    return {**execute_environment(plan), "PYTHONDONTWRITEBYTECODE": "1"}


def identity_inspection_request(
        plan: TransitionPlan, *, profile: str, child_pid: int,
        child_start_tick: int, launcher_pid: int,
        launcher_start_tick: int, controller_pid: int,
        allowed_absence: str,
) -> Dict[str, object]:
    """Validate and canonicalize one two-target identity request."""
    values = (child_pid, launcher_pid, controller_pid)
    ticks = (child_start_tick, launcher_start_tick)
    if profile not in ("original", "replacement") or \
            allowed_absence not in IDENTITY_ABSENCE_POLICIES or \
            any(type(value) is not int or value <= 1 for value in values) or \
            len(set(values)) != len(values) or \
            any(type(value) is not int or value < 0 for value in ticks):
        raise TransitionError("identity inspection request is malformed")
    if profile == "original":
        if child_pid != plan.old_pid or \
                child_start_tick != plan.old_start_tick or \
                launcher_pid != plan.old_launcher_pid or \
                launcher_start_tick != plan.old_launcher_start_tick:
            raise TransitionError(
                "original identity inspection request changed")
    elif allowed_absence not in ("none", "launcher"):
        raise TransitionError(
            "replacement identity absence policy is unsafe")
    return {
        "allowed_absence": allowed_absence,
        "child_pid": child_pid,
        "child_start_tick": child_start_tick,
        "controller_pid": controller_pid,
        "launcher_pid": launcher_pid,
        "launcher_start_tick": launcher_start_tick,
        "profile": profile,
    }


def _identity_executable_binding(
        tools: Mapping[str, Mapping[str, object]], name: str,
) -> Dict[str, int]:
    try:
        binding = tools[name]["binding"]
        device = binding["device"]
        inode = binding["inode"]
    except (KeyError, TypeError) as exc:
        raise TransitionError(
            "identity inspection tool binding is malformed: " + name) from exc
    if type(device) is not int or device < 0 or \
            type(inode) is not int or inode <= 0:
        raise TransitionError(
            "identity inspection tool binding is malformed: " + name)
    return {"device": device, "inode": inode}


def identity_target_contracts(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, Dict[str, object]]:
    """Derive every admitted process field from frozen plan/tool data."""
    canonical = identity_inspection_request(
        plan, profile=str(request.get("profile")),
        child_pid=request.get("child_pid"),
        child_start_tick=request.get("child_start_tick"),
        launcher_pid=request.get("launcher_pid"),
        launcher_start_tick=request.get("launcher_start_tick"),
        controller_pid=request.get("controller_pid"),
        allowed_absence=str(request.get("allowed_absence")))
    if dict(request) != canonical:
        raise TransitionError("identity inspection request is noncanonical")
    child_pid = int(canonical["child_pid"])
    child_tick = int(canonical["child_start_tick"])
    launcher_pid = int(canonical["launcher_pid"])
    launcher_tick = int(canonical["launcher_start_tick"])
    controller_pid = int(canonical["controller_pid"])
    python_executable = _identity_executable_binding(tools, "python3")
    sudo_executable = _identity_executable_binding(tools, "sudo")
    if canonical["profile"] == "original":
        child_argv = plan.old_argv
        child_cmdline_sha256 = plan.old_cmdline_sha256
        child_group = plan.old_process_group
        child_session = plan.old_session
        launcher_argv = plan.old_launch_argv
        launcher_cmdline_sha256 = plan.old_launcher_cmdline_sha256
        launcher_group = plan.old_process_group
        launcher_session = plan.old_session
        launcher_uids = list(plan.old_launcher_uids)
        launcher_affinity: Optional[str] = plan.old_launcher_affinity
        launcher_ppids = [1]
        child_capture_tick = False
    else:
        child_argv = plan.replacement_old_argv
        child_cmdline_sha256 = plan.replacement_old_cmdline_sha256
        child_group = launcher_pid
        child_session = launcher_pid
        launcher_argv = replacement_launch_command(plan, tools)
        launcher_cmdline_sha256 = sha256_bytes(
            b"\0".join(value.encode("ascii") for value in launcher_argv) +
            b"\0")
        launcher_group = launcher_pid
        launcher_session = launcher_pid
        launcher_uids = [plan.controller_uid, 0, 0, 0]
        # The sudo launcher inherits the controller's host affinity.  Its
        # first exact receipt is retained and later receipts compare the full
        # observed identity; the helper still requires a stable nonempty value.
        launcher_affinity = None
        launcher_ppids = [controller_pid, 1]
        child_capture_tick = child_tick == 0
    return {
        "child": {
            "affinity": str(plan.old_cpu),
            "allowed_ppids": [launcher_pid, 1],
            "argv": list(child_argv),
            "capture_start_tick": child_capture_tick,
            "cmdline_sha256": child_cmdline_sha256,
            "executable": python_executable,
            "pid": child_pid,
            "process_group": child_group,
            "session": child_session,
            "start_tick": child_tick,
            "uids": [0, 0, 0, 0],
        },
        "launcher": {
            "affinity": launcher_affinity,
            "allowed_ppids": launcher_ppids,
            "argv": list(launcher_argv),
            "capture_start_tick": False,
            "cmdline_sha256": launcher_cmdline_sha256,
            "executable": sudo_executable,
            "pid": launcher_pid,
            "process_group": launcher_group,
            "session": launcher_session,
            "start_tick": launcher_tick,
            "uids": launcher_uids,
        },
    }


def expected_identity_inspection_orig_argv(
        plan: TransitionPlan, request: Mapping[str, object],
) -> Tuple[str, ...]:
    canonical = identity_inspection_request(
        plan, profile=str(request.get("profile")),
        child_pid=request.get("child_pid"),
        child_start_tick=request.get("child_start_tick"),
        launcher_pid=request.get("launcher_pid"),
        launcher_start_tick=request.get("launcher_start_tick"),
        controller_pid=request.get("controller_pid"),
        allowed_absence=str(request.get("allowed_absence")))
    return (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--inspect-sealed-process-identities", str(canonical["profile"]),
        "--identity-child-pid", str(canonical["child_pid"]),
        "--identity-child-start-tick", str(canonical["child_start_tick"]),
        "--identity-launcher-pid", str(canonical["launcher_pid"]),
        "--identity-launcher-start-tick",
        str(canonical["launcher_start_tick"]),
        "--identity-controller-pid", str(canonical["controller_pid"]),
        "--identity-allowed-absence", str(canonical["allowed_absence"]),
        "--expected-controller-sha256", plan.controller_sha256,
    )


def identity_inspection_command(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Tuple[str, ...]:
    environment = identity_inspection_environment(plan)
    try:
        env_path = str(tools["env"]["path"])
        python_path = str(tools["python3"]["path"])
    except (KeyError, TypeError) as exc:
        raise TransitionError(
            "identity inspection lacks sealed runtime tools") from exc
    orig_argv = expected_identity_inspection_orig_argv(plan, request)
    if orig_argv[0] != python_path:
        raise TransitionError("identity inspection Python path changed")
    return (
        env_path, "-i",
        "HOME=" + environment["HOME"],
        "PATH=" + environment["PATH"],
        "LC_ALL=C", "LANG=C", "TZ=UTC",
        "PYTHONDONTWRITEBYTECODE=1",
        *orig_argv,
    )


def identity_inspection_runtime_contract(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    command = identity_inspection_command(plan, request, tools)
    python_binding = _identity_executable_binding(tools, "python3")
    return {
        "command": list(command),
        "command_nul_sha256": sha256_bytes(
            b"\0".join(value.encode("ascii") for value in command) + b"\0"),
        "controller_interpreter": {
            "argv_path": "/usr/bin/python3.12",
            **python_binding,
            "resolved_path": "/usr/bin/python3.12",
        },
        "environment": identity_inspection_environment(plan),
        "flags": dict(EXECUTE_FLAG_CONTRACT),
        "sys_orig_argv": list(
            expected_identity_inspection_orig_argv(plan, request)),
    }


def _runtime_with_code_seal(
        runtime: Mapping[str, object],
        code_seal: Mapping[str, object],
) -> Dict[str, object]:
    """Bind the complete replayed seal, not merely its controller digest."""
    verified = verify_sealed(
        dict(code_seal), CODE_SEAL_RECEIPT_SCHEMA,
        "identity helper code-seal receipt")
    return {**dict(runtime), "code_seal": verified}


def process_signal_request(
        plan: TransitionPlan, identity_request: Mapping[str, object], *,
        target: str, signum: int,
) -> Dict[str, object]:
    canonical = identity_inspection_request(
        plan, profile=str(identity_request.get("profile")),
        child_pid=identity_request.get("child_pid"),
        child_start_tick=identity_request.get("child_start_tick"),
        launcher_pid=identity_request.get("launcher_pid"),
        launcher_start_tick=identity_request.get("launcher_start_tick"),
        controller_pid=identity_request.get("controller_pid"),
        allowed_absence=str(identity_request.get("allowed_absence")))
    required_absence = {"child": "launcher", "launcher": "child"}
    if isinstance(signum, bool) or not isinstance(signum, int):
        raise TransitionError("process signal number is outside its contract")
    canonical_signum = int(signum)
    if canonical["profile"] != "original" or target not in required_absence or \
            canonical["allowed_absence"] != required_absence[target] or \
            canonical_signum not in (
                int(signal.SIGTERM), int(signal.SIGKILL)):
        raise TransitionError("process signal request is outside its contract")
    return {
        "identity_request": canonical,
        "signal": canonical_signum,
        "signal_name": signal.Signals(canonical_signum).name,
        "target": target,
    }


def expected_process_signal_orig_argv(
        plan: TransitionPlan, request: Mapping[str, object],
) -> Tuple[str, ...]:
    signal_request = process_signal_request(
        plan, request["identity_request"], target=str(request.get("target")),
        signum=request.get("signal"))
    identity = signal_request["identity_request"]
    return (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--signal-sealed-process-identities", str(identity["profile"]),
        "--identity-child-pid", str(identity["child_pid"]),
        "--identity-child-start-tick", str(identity["child_start_tick"]),
        "--identity-launcher-pid", str(identity["launcher_pid"]),
        "--identity-launcher-start-tick",
        str(identity["launcher_start_tick"]),
        "--identity-controller-pid", str(identity["controller_pid"]),
        "--identity-allowed-absence", str(identity["allowed_absence"]),
        "--signal-target", str(signal_request["target"]),
        "--signal-number", str(signal_request["signal"]),
        "--expected-controller-sha256", plan.controller_sha256,
    )


def process_signal_command(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Tuple[str, ...]:
    environment = identity_inspection_environment(plan)
    try:
        env_path = str(tools["env"]["path"])
        python_path = str(tools["python3"]["path"])
    except (KeyError, TypeError) as exc:
        raise TransitionError("process signal lacks sealed runtime tools") \
            from exc
    orig_argv = expected_process_signal_orig_argv(plan, request)
    if orig_argv[0] != python_path:
        raise TransitionError("process signal Python path changed")
    return (
        env_path, "-i", "HOME=" + environment["HOME"],
        "PATH=" + environment["PATH"], "LC_ALL=C", "LANG=C", "TZ=UTC",
        "PYTHONDONTWRITEBYTECODE=1", *orig_argv,
    )


def process_signal_runtime_contract(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    command = process_signal_command(plan, request, tools)
    return {
        "command": list(command),
        "command_nul_sha256": sha256_bytes(
            b"\0".join(value.encode("ascii") for value in command) + b"\0"),
        "controller_interpreter": {
            "argv_path": "/usr/bin/python3.12",
            **_identity_executable_binding(tools, "python3"),
            "resolved_path": "/usr/bin/python3.12",
        },
        "environment": identity_inspection_environment(plan),
        "flags": dict(EXECUTE_FLAG_CONTRACT),
        "sys_orig_argv": list(expected_process_signal_orig_argv(plan, request)),
    }


def identity_inspection_plan_contract(
        plan: TransitionPlan,
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    """Return historical helper-contract data without dispatch authority."""
    environment = identity_inspection_environment(plan)
    python_binding = _identity_executable_binding(tools, "python3")
    try:
        env_path = str(tools["env"]["path"])
        python_path = str(tools["python3"]["path"])
        sudo_path = str(tools["sudo"]["path"])
        timeout_path = str(tools["timeout"]["path"])
    except (KeyError, TypeError) as exc:
        raise TransitionError(
            "identity inspection plan lacks sealed runtime tools") from exc
    if python_path != "/usr/bin/python3.12":
        raise TransitionError("identity inspection plan Python path changed")
    return {
        "absence_policies": {
            "original": ["both", "child", "launcher", "none"],
            "replacement": ["launcher", "none"],
        },
        "argv_contract": {
            "environment_prefix": [
                env_path, "-i", "HOME=" + environment["HOME"],
                "PATH=" + environment["PATH"], "LC_ALL=C", "LANG=C",
                "TZ=UTC", "PYTHONDONTWRITEBYTECODE=1",
            ],
            "expected_controller_sha256": plan.controller_sha256,
            "ordered_request_options": [
                "--inspect-sealed-process-identities",
                "--identity-child-pid", "--identity-child-start-tick",
                "--identity-launcher-pid", "--identity-launcher-start-tick",
                "--identity-controller-pid", "--identity-allowed-absence",
                "--expected-controller-sha256",
            ],
            "python_prefix": [
                python_path, "-I", "-S", "-B", str(plan.controller),
            ],
        },
        "controller_interpreter": {
            "argv_path": python_path,
            **python_binding,
            "resolved_path": python_path,
        },
        "code_seal": {
            "full_receipt_in_helper_runtime_and_receipt": True,
            "replayed_before_target_proc": True,
            "schema": CODE_SEAL_RECEIPT_SCHEMA,
        },
        "environment": environment,
        "flags": dict(EXECUTE_FLAG_CONTRACT),
        "pidfd_policy": IDENTITY_PIDFD_POLICY,
        "process_signal": {
            "argv_contract": {
                "environment_prefix": [
                    env_path, "-i", "HOME=" + environment["HOME"],
                    "PATH=" + environment["PATH"], "LC_ALL=C", "LANG=C",
                    "TZ=UTC", "PYTHONDONTWRITEBYTECODE=1",
                ],
                "ordered_request_options": [
                    "--signal-sealed-process-identities",
                    "--identity-child-pid", "--identity-child-start-tick",
                    "--identity-launcher-pid",
                    "--identity-launcher-start-tick",
                    "--identity-controller-pid",
                    "--identity-allowed-absence", "--signal-target",
                    "--signal-number", "--expected-controller-sha256",
                ],
                "python_prefix": [
                    python_path, "-I", "-S", "-B", str(plan.controller),
                ],
            },
            "pidfd_policy": PROCESS_SIGNAL_PIDFD_POLICY,
            "receipt_schema": PROCESS_SIGNAL_SCHEMA,
            "signals": [int(signal.SIGTERM), int(signal.SIGKILL)],
            "targets": ["child", "launcher"],
            "recovery_action_evidence_schema":
                RECOVERY_ACTION_EVIDENCE_SCHEMA,
        },
        "privileged_wrapper_tools": {
            "sudo": sudo_path,
            "timeout": timeout_path,
        },
        "receipt_schema": IDENTITY_INSPECTION_SCHEMA,
    }


def _validate_identity_against_contract(
        identity: Mapping[str, object], contract: Mapping[str, object],
        description: str,
) -> None:
    if set(identity) != PROCESS_IDENTITY_KEYS:
        raise TransitionError(
            description + " process identity field roster changed")
    if re.fullmatch(r"[RSDTtKWPI]", str(identity.get("state"))) is None or \
            not _processor_is_allowed(
                identity.get("processor"), str(identity.get("affinity"))):
        raise TransitionError(
            description + " process live stat fields are malformed")
    exact_keys = (
        "argv", "cmdline_sha256", "executable", "pid", "process_group",
        "session", "uids")
    for key in exact_keys:
        if identity.get(key) != contract.get(key):
            raise TransitionError(
                description + " process identity changed: " + key)
    start_tick = identity.get("start_tick")
    if contract.get("capture_start_tick") is True:
        if type(start_tick) is not int or start_tick <= 0:
            raise TransitionError(
                description + " process start tick is invalid")
    elif start_tick != contract.get("start_tick"):
        raise TransitionError(
            description + " process identity changed: start_tick")
    if identity.get("ppid") not in contract.get("allowed_ppids", ()):
        raise TransitionError(
            description + " process identity changed: ppid")
    expected_affinity = contract.get("affinity")
    if expected_affinity is None:
        if not isinstance(identity.get("affinity"), str) or \
                not identity["affinity"]:
            raise TransitionError(
                description + " process affinity is malformed")
    elif identity.get("affinity") != expected_affinity:
        raise TransitionError(
            description + " process identity changed: affinity")


def stable_process_identity(
        identity: Mapping[str, object],
) -> Dict[str, object]:
    """Project out scheduler state that is evidence, not process identity."""
    if set(identity) != PROCESS_IDENTITY_KEYS:
        raise TransitionError("process identity field roster changed")
    return {key: identity[key] for key in sorted(STABLE_PROCESS_IDENTITY_KEYS)}


def _processor_is_allowed(processor: int, affinity: str) -> bool:
    if type(processor) is not int or processor < 0 or not isinstance(
            affinity, str) or not affinity:
        return False
    permitted = False
    for part in affinity.split(","):
        match = re.fullmatch(r"(0|[1-9][0-9]*)(?:-(0|[1-9][0-9]*))?", part)
        if match is None:
            return False
        first = int(match.group(1))
        last = first if match.group(2) is None else int(match.group(2))
        if last < first or last > 1048575:
            return False
        permitted = permitted or first <= processor <= last
    return permitted


def _validate_proc_stat_against_contract(
        value: Mapping[str, object], contract: Mapping[str, object],
        description: str,
) -> None:
    """Bind the stat fields that survive executable/cmdline teardown."""
    if set(value) != PROC_STAT_IDENTITY_KEYS or \
            not isinstance(value.get("state"), str) or \
            len(value["state"]) != 1 or \
            value.get("process_group") != contract.get("process_group") or \
            value.get("session") != contract.get("session") or \
            value.get("ppid") not in contract.get("allowed_ppids", ()):
        raise TransitionError(
            description + " process stat identity changed")
    start_tick = value.get("start_tick")
    if contract.get("capture_start_tick") is True:
        if type(start_tick) is not int or start_tick <= 0:
            raise TransitionError(
                description + " process stat start tick is invalid")
    elif start_tick != contract.get("start_tick"):
        raise TransitionError(
            description + " process stat start tick changed")
    affinity = contract.get("affinity")
    if affinity is not None and not _processor_is_allowed(
            value.get("processor"), str(affinity)):
        raise TransitionError(
            description + " process stat processor changed")


def _identity_proc_stat(identity: Mapping[str, object]) -> Dict[str, object]:
    if set(identity) != PROCESS_IDENTITY_KEYS:
        raise TransitionError("process identity field roster changed")
    return {key: identity[key] for key in sorted(PROC_STAT_IDENTITY_KEYS)}


def _read_proc_stat_at(proc_root: Path, pid: int) -> Dict[str, object]:
    return _parse_proc_stat((proc_root / str(pid) / "stat").read_bytes())


def _identity_matches_proc_stat(
        identity: Mapping[str, object], value: Mapping[str, object],
) -> bool:
    return all(identity.get(key) == value.get(key) for key in (
        "ppid", "process_group", "session", "start_tick"))


def _exiting_identity_target(
        contract: Mapping[str, object], *, snapshot_kind: str,
        proc_stat: Optional[Mapping[str, object]],
        identity: Optional[Mapping[str, object]], observations: int,
        pidfd_ready: bool = True,
) -> Dict[str, object]:
    if snapshot_kind not in {
            "exit-stat-before-pidfd-ready",
            "pidfd-ready-after-identity", "pidfd-ready-proc-disappeared",
            "pidfd-ready-with-stat"} or \
            type(observations) is not int or not 0 <= observations <= 2 or \
            type(pidfd_ready) is not bool:
        raise TransitionError("exiting process receipt is malformed")
    stat_value = dict(proc_stat) if proc_stat is not None else None
    identity_value = dict(identity) if identity is not None else None
    if stat_value is not None:
        _validate_proc_stat_against_contract(
            stat_value, contract, "exiting")
    if identity_value is not None:
        _validate_identity_against_contract(
            identity_value, contract, "exiting")
        if stat_value is None or not _identity_matches_proc_stat(
                identity_value, stat_value):
            raise TransitionError(
                "exiting process identity/stat evidence changed")
    if snapshot_kind == "exit-stat-before-pidfd-ready":
        malformed_kind = pidfd_ready or stat_value is None or \
            identity_value is not None or observations != 0 or \
            stat_value.get("state") not in EXIT_PROCESS_STATES
    elif not pidfd_ready:
        malformed_kind = True
    elif snapshot_kind == "pidfd-ready-proc-disappeared":
        malformed_kind = stat_value is not None or identity_value is not None \
            or observations != 0
    elif snapshot_kind == "pidfd-ready-after-identity":
        malformed_kind = stat_value is None or identity_value is None or \
            observations < 1
    else:
        malformed_kind = stat_value is None or identity_value is not None
    if malformed_kind:
        raise TransitionError("exiting process evidence kind changed")
    return {
        "expected_pid": contract["pid"],
        "expected_start_tick": contract["start_tick"],
        "identity": identity_value,
        "pidfd_opened_before_target_proc": True,
        "pidfd_ready_before_receipt": pidfd_ready,
        "proc_identity_observations": observations,
        "proc_stat": stat_value,
        "snapshot_kind": snapshot_kind,
        "start_tick_proven": stat_value is not None,
        "state": "exiting",
    }


def _pidfd_is_ready(pidfd: int) -> bool:
    poller = select.poll()
    poller.register(pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
    return bool(poller.poll(0))


def _pidfd_send_signal(pidfd: int, signum: int) -> None:
    signal.pidfd_send_signal(pidfd, signum, None, 0)


def capture_process_identity_targets(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]], *,
        proc_root: Path = Path("/proc"),
        boot_id_path: Path = Path("/proc/sys/kernel/random/boot_id"),
        pidfd_open: Optional[Callable[[int, int], int]] = None,
        pidfd_ready: Callable[[int], bool] = _pidfd_is_ready,
        pidfd_close: Callable[[int], None] = os.close,
        proc_stat_reader: Optional[
            Callable[[Path, int], Mapping[str, object]]] = None,
        signal_target: Optional[str] = None,
        signum: Optional[int] = None,
        pidfd_send_signal: Callable[[int, int], None] = _pidfd_send_signal,
) -> Dict[str, object]:
    """Open the complete PID roster first, then double-capture exact state.

    When a signal is requested, the validated target descriptor stays open
    from before the first target-proc read through ``pidfd_send_signal``.  An
    exit snapshot is evidence about that exact pidfd, never proof of absence.
    """
    canonical = identity_inspection_request(
        plan, profile=str(request.get("profile")),
        child_pid=request.get("child_pid"),
        child_start_tick=request.get("child_start_tick"),
        launcher_pid=request.get("launcher_pid"),
        launcher_start_tick=request.get("launcher_start_tick"),
        controller_pid=request.get("controller_pid"),
        allowed_absence=str(request.get("allowed_absence")))
    contracts = identity_target_contracts(plan, canonical, tools)
    allowed = IDENTITY_ABSENCE_POLICIES[str(canonical["allowed_absence"])]
    opener = os.pidfd_open if pidfd_open is None else pidfd_open
    canonical_signum = None if signum is None else \
        int(signum) if isinstance(signum, int) and not isinstance(signum, bool) \
        else -1
    if (signal_target is None) != (signum is None) or \
            signal_target is not None and (
                signal_target not in ("child", "launcher") or
                canonical_signum not in (
                    int(signal.SIGTERM), int(signal.SIGKILL)) or
                signal_target in allowed):
        raise TransitionError("pidfd signal capture request is malformed")
    opened: Dict[str, int] = {}
    absent = set()
    captured: Dict[str, Dict[str, object]] = {}
    signal_delivered = False
    stat_reader = _read_proc_stat_at if proc_stat_reader is None else \
        proc_stat_reader

    def absent_error(exc: OSError) -> bool:
        return exc.errno in (errno.ENOENT, errno.ESRCH)

    def capture_one(
            name: str,
            prior: Optional[Mapping[str, object]] = None,
    ) -> Dict[str, object]:
        contract = contracts[name]
        pid = int(contract["pid"])
        try:
            stat_value = dict(stat_reader(proc_root, pid))
        except OSError as exc:
            if not absent_error(exc):
                raise TransitionError(
                    name + " process stat became unavailable") from exc
            if not pidfd_ready(opened[name]):
                raise TransitionError(
                    name + " process stat vanished with a live pidfd") \
                    from exc
            if prior is not None:
                return _exiting_identity_target(
                    contract, snapshot_kind="pidfd-ready-after-identity",
                    proc_stat=prior["proc_stat"], identity=prior["identity"],
                    observations=int(prior["observations"]))
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-proc-disappeared",
                proc_stat=None, identity=None, observations=0)
        _validate_proc_stat_against_contract(stat_value, contract, name)
        ready_after_stat = pidfd_ready(opened[name])
        if stat_value["state"] in EXIT_PROCESS_STATES:
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-with-stat"
                if ready_after_stat else "exit-stat-before-pidfd-ready",
                proc_stat=stat_value, identity=None, observations=0,
                pidfd_ready=ready_after_stat)
        if ready_after_stat:
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-with-stat",
                proc_stat=stat_value, identity=None, observations=0)
        try:
            identity = _proc_identity_at(proc_root, pid)
        except ProcessSnapshotRace as exc:
            if not pidfd_ready(opened[name]):
                raise TransitionError(
                    name + " process snapshot tore with a live pidfd") from exc
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-with-stat",
                proc_stat=stat_value, identity=None, observations=0)
        except OSError as exc:
            if not absent_error(exc) or not pidfd_ready(opened[name]):
                raise TransitionError(
                    name + " process identity became unavailable") from exc
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-with-stat",
                proc_stat=stat_value, identity=None, observations=0)
        _validate_identity_against_contract(identity, contract, name)
        if not _identity_matches_proc_stat(identity, stat_value):
            raise TransitionError(
                name + " process stat changed during identity capture")
        observations = 1
        if prior is not None:
            if stable_process_identity(identity) != \
                    stable_process_identity(prior["identity"]):
                raise TransitionError(
                    name + " process identity changed during inspection")
            observations = int(prior["observations"]) + 1
        if pidfd_ready(opened[name]):
            return _exiting_identity_target(
                contract, snapshot_kind="pidfd-ready-after-identity",
                proc_stat=stat_value, identity=identity,
                observations=observations)
        if observations not in (1, 2):
            raise TransitionError(
                name + " process identity observation count changed")
        return {
            "identity": identity,
            "observations": observations,
            "proc_stat": stat_value,
            "state": "present-internal",
        }

    try:
        # No target proc path is touched until both pidfd_open attempts have
        # completed.  A permission or semantic failure is never absence.
        for name in ("child", "launcher"):
            pid = int(contracts[name]["pid"])
            try:
                opened[name] = opener(pid, 0)
            except OSError as exc:
                if name not in allowed or not absent_error(exc):
                    raise TransitionError(
                        name + " pidfd open failed without allowed absence") \
                        from exc
                absent.add(name)
        boot_before = boot_id_path.read_text(encoding="ascii").strip()
        if boot_before != plan.old_boot_id:
            raise TransitionError("identity inspection boot ID changed")
        for name in ("child", "launcher"):
            if name in absent:
                continue
            captured[name] = capture_one(name)
        for name in ("child", "launcher"):
            if name in absent or captured[name]["state"] == "exiting":
                continue
            second = capture_one(name, captured[name])
            if second["state"] == "exiting":
                captured[name] = second
                continue
            captured[name] = second
        boot_after = boot_id_path.read_text(encoding="ascii").strip()
        if boot_after != boot_before:
            raise TransitionError(
                "identity inspection boot ID changed during replay")
        targets: Dict[str, object] = {}
        for name in ("child", "launcher"):
            contract = contracts[name]
            if name not in absent and captured[name]["state"] == "exiting":
                targets[name] = captured[name]
                continue
            if name not in absent:
                targets[name] = {
                    "identity": captured[name]["identity"],
                    "pidfd_opened_before_target_proc": True,
                    "pidfd_unready_after_replay": True,
                    "proc_identity_observations":
                        captured[name]["observations"],
                    "state": "present",
                }
                continue
            pid = int(contract["pid"])
            try:
                probe = opener(pid, 0)
            except OSError as exc:
                if not absent_error(exc):
                    raise TransitionError(
                        name + " absence recheck failed") from exc
            else:
                pidfd_close(probe)
                raise TransitionError(
                    name + " PID appeared during absence inspection")
            targets[name] = {
                "expected_pid": pid,
                "expected_start_tick": int(contract["start_tick"]),
                "pidfd_absence_observations": 2,
                "state": "absent",
            }
        # This is deliberately the last observation before constructing the
        # return value.  Checking a child earlier and then inspecting its
        # launcher would allow the child to exit while still being described
        # as unready at receipt time.
        for name in ("child", "launcher"):
            if name not in opened:
                continue
            ready = pidfd_ready(opened[name])
            if targets[name]["state"] == "exiting":
                was_ready = targets[name]["pidfd_ready_before_receipt"]
                if was_ready and not ready:
                    raise TransitionError(
                        name + " exiting pidfd readiness regressed")
                if ready and not was_ready:
                    targets[name] = _exiting_identity_target(
                        contracts[name], snapshot_kind="pidfd-ready-with-stat",
                        proc_stat=targets[name]["proc_stat"], identity=None,
                        observations=0)
                continue
            if ready:
                targets[name] = _exiting_identity_target(
                    contracts[name],
                    snapshot_kind="pidfd-ready-after-identity",
                    proc_stat=captured[name]["proc_stat"],
                    identity=captured[name]["identity"],
                    observations=captured[name]["observations"])
        result: Dict[str, object] = {
            "boot_id": boot_before, "targets": targets}
        if signal_target is not None:
            target = targets[signal_target]
            if target.get("state") != "present":
                raise TransitionError(
                    signal_target + " exited before exact pidfd signal")
            target_fd = opened[signal_target]
            boot_at_signal = boot_id_path.read_text(encoding="ascii").strip()
            if boot_at_signal != boot_before:
                raise TransitionError(
                    "identity inspection boot ID changed before signal")
            try:
                pidfd_send_signal(target_fd, canonical_signum)
            except OSError as exc:
                raise TransitionError(
                    "pidfd_send_signal failed for exact " + signal_target) \
                    from exc
            signal_delivered = True
            result["signal_action"] = {
                "boot_id": boot_at_signal,
                "identity_target_sha256": sha256_bytes(
                    canonical_json(target)),
                "pid": contracts[signal_target]["pid"],
                "pidfd": target_fd,
                "pidfd_held_from_before_target_proc_through_signal": True,
                "pidfd_send_signal_result": 0,
                "signal": canonical_signum,
                "signal_name": signal.Signals(canonical_signum).name,
                "syscall": "pidfd_send_signal",
                "target": signal_target,
            }
        return result
    finally:
        close_error: Optional[OSError] = None
        for pidfd in opened.values():
            try:
                pidfd_close(pidfd)
            except OSError as exc:
                if close_error is None:
                    close_error = exc
        # After a successful signal, a descriptor-close failure in this
        # short-lived helper must not erase the canonical action receipt.  The
        # process exits immediately after stdout and the kernel closes it.
        if close_error is not None and not signal_delivered:
            raise TransitionError("identity inspection pidfd close failed") \
                from close_error


def capture_process_identity_targets_under_seal(
        plan: TransitionPlan, request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]], **capture_options: object,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    """Replay the complete code seal before any target process is touched."""
    code_seal = verify_code_seal(plan)
    captured = capture_process_identity_targets(
        plan, request, tools, **capture_options)
    return code_seal, captured


def _verify_identity_inspection_receipt_with_code_seal(
        plan: TransitionPlan, value: object,
        request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
        code_seal: Mapping[str, object],
) -> Dict[str, object]:
    code_seal = verify_sealed(
        dict(code_seal), CODE_SEAL_RECEIPT_SCHEMA,
        "process identity inspection code seal")
    receipt = verify_sealed(
        value, IDENTITY_INSPECTION_SCHEMA, "process identity inspection receipt")
    required_keys = {
        "boot_id", "code_seal", "helper_runtime", "pidfd_policy", "request",
        "schema", "self_sha256_excluding_field", "targets", "tools",
        "transition_id"}
    canonical = identity_inspection_request(
        plan, profile=str(request.get("profile")),
        child_pid=request.get("child_pid"),
        child_start_tick=request.get("child_start_tick"),
        launcher_pid=request.get("launcher_pid"),
        launcher_start_tick=request.get("launcher_start_tick"),
        controller_pid=request.get("controller_pid"),
        allowed_absence=str(request.get("allowed_absence")))
    contracts = identity_target_contracts(plan, canonical, tools)
    if set(receipt) != required_keys or receipt.get("request") != canonical or \
            receipt.get("boot_id") != plan.old_boot_id or \
            receipt.get("transition_id") != plan.transition_id or \
            receipt.get("code_seal") != code_seal or \
            receipt.get("tools") != dict(tools) or \
            receipt.get("helper_runtime") != \
                _runtime_with_code_seal(
                    identity_inspection_runtime_contract(
                        plan, canonical, tools), code_seal) or \
            receipt.get("pidfd_policy") != IDENTITY_PIDFD_POLICY:
        raise TransitionError(
            "process identity inspection receipt contract changed")
    targets = receipt.get("targets")
    allowed = IDENTITY_ABSENCE_POLICIES[str(canonical["allowed_absence"])]
    if not isinstance(targets, dict) or set(targets) != {"child", "launcher"}:
        raise TransitionError(
            "process identity inspection target roster changed")
    for name in ("child", "launcher"):
        target = targets[name]
        if not isinstance(target, dict):
            raise TransitionError(
                "process identity inspection target is malformed: " + name)
        if target.get("state") == "present":
            if set(target) != {
                    "identity", "pidfd_opened_before_target_proc",
                    "pidfd_unready_after_replay", "proc_identity_observations",
                    "state"} or \
                    target.get("pidfd_opened_before_target_proc") is not True or \
                    target.get("pidfd_unready_after_replay") is not True or \
                    target.get("proc_identity_observations") != 2 or \
                    not isinstance(target.get("identity"), dict):
                raise TransitionError(
                    "present identity target receipt is malformed: " + name)
            _validate_identity_against_contract(
                target["identity"], contracts[name], name)
        elif target.get("state") == "exiting":
            if set(target) != {
                    "expected_pid", "expected_start_tick", "identity",
                    "pidfd_opened_before_target_proc",
                    "pidfd_ready_before_receipt",
                    "proc_identity_observations", "proc_stat",
                    "snapshot_kind", "start_tick_proven", "state"}:
                raise TransitionError(
                    "exiting identity target receipt is malformed: " + name)
            expected = _exiting_identity_target(
                contracts[name], snapshot_kind=str(target["snapshot_kind"]),
                proc_stat=target["proc_stat"], identity=target["identity"],
                observations=target["proc_identity_observations"],
                pidfd_ready=target["pidfd_ready_before_receipt"])
            if target != expected:
                raise TransitionError(
                    "exiting identity target receipt changed: " + name)
        elif target.get("state") == "absent":
            if name not in allowed or set(target) != {
                    "expected_pid", "expected_start_tick",
                    "pidfd_absence_observations", "state"} or \
                    target.get("expected_pid") != contracts[name]["pid"] or \
                    target.get("expected_start_tick") != \
                        contracts[name]["start_tick"] or \
                    target.get("pidfd_absence_observations") != 2:
                raise TransitionError(
                    "absent identity target receipt is unsafe: " + name)
        else:
            raise TransitionError(
                "process identity inspection target state changed: " + name)
    return receipt


def verify_identity_inspection_receipt(
        plan: TransitionPlan, value: object,
        request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    return _verify_identity_inspection_receipt_with_code_seal(
        plan, value, request, tools, verify_code_seal(plan))


def _verify_process_signal_receipt_with_code_seal(
        plan: TransitionPlan, value: object,
        request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
        code_seal: Mapping[str, object],
) -> Dict[str, object]:
    code_seal = verify_sealed(
        dict(code_seal), CODE_SEAL_RECEIPT_SCHEMA,
        "exact process signal code seal")
    canonical = process_signal_request(
        plan, request["identity_request"], target=str(request.get("target")),
        signum=request.get("signal"))
    receipt = verify_sealed(
        value, PROCESS_SIGNAL_SCHEMA, "exact process signal receipt")
    if set(receipt) != {
            "boot_id", "code_seal", "helper_runtime", "pidfd_policy",
            "request", "schema", "self_sha256_excluding_field",
            "signal_action", "targets", "tools", "transition_id"} or \
            receipt.get("boot_id") != plan.old_boot_id or \
            receipt.get("code_seal") != code_seal or \
            receipt.get("helper_runtime") != _runtime_with_code_seal(
                process_signal_runtime_contract(
                    plan, canonical, tools), code_seal) or \
            receipt.get("pidfd_policy") != PROCESS_SIGNAL_PIDFD_POLICY or \
            receipt.get("request") != canonical or \
            receipt.get("tools") != dict(tools) or \
            receipt.get("transition_id") != plan.transition_id:
        raise TransitionError("exact process signal receipt contract changed")
    # Reuse the exact target verifier through a canonical inspection
    # projection.  This does not claim that a separate inspection occurred;
    # both projections describe the pidfds held by the signal helper itself.
    identity_projection = sealed_record(IDENTITY_INSPECTION_SCHEMA, {
        "boot_id": receipt["boot_id"],
        "code_seal": code_seal,
        "helper_runtime": _runtime_with_code_seal(
            identity_inspection_runtime_contract(
                plan, canonical["identity_request"], tools), code_seal),
        "pidfd_policy": IDENTITY_PIDFD_POLICY,
        "request": canonical["identity_request"],
        "targets": receipt["targets"],
        "tools": dict(tools),
        "transition_id": plan.transition_id,
    })
    _verify_identity_inspection_receipt_with_code_seal(
        plan, identity_projection, canonical["identity_request"], tools,
        code_seal)
    action = receipt.get("signal_action")
    target = receipt["targets"][canonical["target"]]
    contracts = identity_target_contracts(
        plan, canonical["identity_request"], tools)
    if not isinstance(action, dict) or set(action) != {
            "boot_id", "identity_target_sha256", "pid", "pidfd",
            "pidfd_held_from_before_target_proc_through_signal",
            "pidfd_send_signal_result", "signal", "signal_name", "syscall",
            "target"} or target.get("state") != "present" or \
            action.get("boot_id") != plan.old_boot_id or \
            action.get("identity_target_sha256") != sha256_bytes(
                canonical_json(target)) or \
            action.get("pid") != contracts[canonical["target"]]["pid"] or \
            type(action.get("pidfd")) is not int or action["pidfd"] < 0 or \
            action.get(
                "pidfd_held_from_before_target_proc_through_signal") is not True or \
            action.get("pidfd_send_signal_result") != 0 or \
            action.get("signal") != canonical["signal"] or \
            action.get("signal_name") != canonical["signal_name"] or \
            action.get("syscall") != "pidfd_send_signal" or \
            action.get("target") != canonical["target"]:
        raise TransitionError("exact process signal action receipt changed")
    return receipt


def verify_process_signal_receipt(
        plan: TransitionPlan, value: object,
        request: Mapping[str, object],
        tools: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    return _verify_process_signal_receipt_with_code_seal(
        plan, value, request, tools, verify_code_seal(plan))


class Deadline:
    def __init__(self, duration_s: float, recovery_reserve_s: float,
                 *, clock: Callable[[], float] = time.monotonic) -> None:
        if not math.isfinite(duration_s) or not math.isfinite(recovery_reserve_s) or \
                duration_s <= recovery_reserve_s or recovery_reserve_s <= 0:
            raise ValueError("transition deadline is malformed")
        self._clock = clock
        self.started = clock()
        self.absolute = self.started + duration_s
        self.normal = self.absolute - recovery_reserve_s

    def remaining(self) -> float:
        return max(0.0, self.absolute - self._clock())

    def now(self) -> float:
        return self._clock()

    def exhausted(self) -> bool:
        return self._clock() >= self.absolute

    def require_normal(self, phase: str) -> None:
        if self._clock() >= self.normal:
            raise TransitionDeadline(
                "normal deadline reached before phase " + phase)

    def start_emergency_recovery(
            self, duration_s: float) -> "EmergencyRecoveryBudget":
        return EmergencyRecoveryBudget(duration_s, clock=self._clock)


class EmergencyRecoveryBudget:
    """Fresh post-stop budget; it bounds waits but never vetoes safe actions."""

    def __init__(self, duration_s: float, *,
                 clock: Callable[[], float] = time.monotonic) -> None:
        if not math.isfinite(duration_s) or duration_s <= 0:
            raise ValueError("emergency recovery budget is malformed")
        self._clock = clock
        self.started = clock()
        self.absolute = self.started + duration_s
        self.safety_extension_count = 0
        self.safety_extension_deadline_max: Optional[float] = None

    def now(self) -> float:
        return self._clock()

    def remaining(self) -> float:
        return max(0.0, self.absolute - self._clock())

    def exhausted(self) -> bool:
        return self._clock() >= self.absolute

    def wait_deadline(self, maximum_wait_s: float, *,
                      minimum_safety_wait_s: float = 0.05) -> float:
        if not math.isfinite(maximum_wait_s) or maximum_wait_s <= 0 or \
                not math.isfinite(minimum_safety_wait_s) or \
                minimum_safety_wait_s <= 0 or \
                minimum_safety_wait_s > maximum_wait_s:
            raise ValueError("recovery wait bound is malformed")
        now = self._clock()
        deadline = max(
            min(self.absolute, now + maximum_wait_s),
            now + minimum_safety_wait_s)
        if deadline > self.absolute:
            self.safety_extension_count += 1
            self.safety_extension_deadline_max = max(
                deadline, self.safety_extension_deadline_max or deadline)
        return deadline

    def receipt(self) -> Dict[str, object]:
        observed = self._clock()
        return {
            "absolute_monotonic_s": self.absolute,
            "exhausted": observed >= self.absolute,
            "observed_monotonic_s": observed,
            "remaining_s": max(0.0, self.absolute - observed),
            "safety_extension_count": self.safety_extension_count,
            "safety_extension_deadline_max":
                self.safety_extension_deadline_max,
            "started_monotonic_s": self.started,
        }


class ReceiptJournal:
    def __init__(self, directory: Path, transition_id: str, owner_uid: int,
                 deadline: Deadline) -> None:
        self.directory = directory
        self.transition_id = transition_id
        self.owner_uid = owner_uid
        self.deadline = deadline
        self.sequence = 0
        self.previous_sha256: Optional[str] = None

    def record(self, phase: str, status_value: str,
               payload: Mapping[str, object]) -> Dict[str, object]:
        if not re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", phase) or \
                status_value not in {"started", "completed", "failed"}:
            raise TransitionError("phase receipt identity is malformed")
        sequence = self.sequence
        record = sealed_record(PHASE_SCHEMA, {
            "absolute_deadline_monotonic_s": self.deadline.absolute,
            "created_utc": utc_now(),
            "payload": dict(payload),
            "phase": phase,
            "previous_receipt_sha256": self.previous_sha256,
            "remaining_s": self.deadline.remaining(),
            "sequence": sequence,
            "status": status_value,
            "transition_id": self.transition_id,
        })
        raw = canonical_json(record)
        path = self.directory / (
            "%04d-%s-%s.json" % (sequence, phase, status_value))
        binding = write_new(path, raw, owner_uid=self.owner_uid)
        self.sequence += 1
        self.previous_sha256 = str(binding["sha256"])
        return record

    def replay(self) -> Dict[str, object]:
        """Re-read the complete durable prefix and verify its hash chain."""
        paths = sorted(self.directory.glob("*.json"))
        if len(paths) != self.sequence:
            raise TransitionError("phase receipt roster is incomplete or enlarged")
        previous: Optional[str] = None
        roster = []
        for sequence, path in enumerate(paths):
            match = re.fullmatch(
                r"([0-9]{4})-([a-z][a-z0-9_-]{0,63})-"
                r"(started|completed|failed)\.json", path.name)
            if match is None or int(match.group(1)) != sequence:
                raise TransitionError("phase receipt filename is noncanonical")
            # Reject FIFOs/devices/symlinks without first opening them through
            # the path-following stable-file reader.  The private parent is the
            # namespace trust boundary; this ordering also makes a malformed
            # receipt roster fail promptly instead of blocking recovery.
            binding = file_binding(path, with_hash=True)
            value = load_canonical(path, PHASE_SCHEMA, "phase receipt")
            if binding["uid"] != self.owner_uid or \
                    binding["mode"] != 0o444 or binding["nlink"] != 1 or \
                    value.get("sequence") != sequence or \
                    value.get("phase") != match.group(2) or \
                    value.get("status") != match.group(3) or \
                    value.get("transition_id") != self.transition_id or \
                    value.get("absolute_deadline_monotonic_s") != \
                        self.deadline.absolute or \
                    value.get("previous_receipt_sha256") != previous:
                raise TransitionError("phase receipt chain changed")
            roster.append({"binding": binding, "path": str(path),
                           "phase": match.group(2),
                           "sequence": sequence, "status": match.group(3)})
            previous = str(binding["sha256"])
        if previous != self.previous_sha256:
            raise TransitionError("phase receipt chain head changed")
        return {"count": len(roster), "head_sha256": previous,
                "roster": roster}


class DeferredSignals:
    def __init__(self) -> None:
        self.previous: Dict[int, object] = {}
        self.requested: list[int] = []

    def __enter__(self) -> "DeferredSignals":
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            try:
                self.previous[int(signum)] = signal.getsignal(signum)
                signal.signal(signum, self._record)
            except (OSError, ValueError) as exc:
                for installed, previous in self.previous.items():
                    signal.signal(installed, previous)
                raise TransitionError(
                    "cannot install transition signal guard") from exc
        return self

    def __exit__(self, _kind: object, _value: object, _trace: object) -> None:
        for signum, previous in self.previous.items():
            signal.signal(signum, previous)
        self.previous.clear()
        if _kind is None:
            self.raise_if_requested()

    def _record(self, signum: int, _frame: object) -> None:
        if not self.requested:
            self.requested.append(int(signum))

    def raise_if_requested(self) -> None:
        if self.requested:
            raise TransitionError(
                "deferred controller signal: " +
                signal.Signals(self.requested[0]).name)


def error_record(exc: BaseException) -> Dict[str, str]:
    return {"message": str(exc), "type": type(exc).__name__}


class TransitionController:
    """Retired state machine retained for exhaustive source-only tests."""

    def __init__(self, plan: TransitionPlan, backend: object,
                 journal: ReceiptJournal, deadline: Deadline,
                 signal_guard: Optional[DeferredSignals] = None) -> None:
        self.plan = plan
        self.backend = backend
        self.journal = journal
        self.deadline = deadline
        self.signal_guard = signal_guard
        self.recovery_armed = False
        self.dry_accepted = False
        self.primary_error: Optional[BaseException] = None
        self.recovery_errors: list[BaseException] = []
        self.old_stop: Optional[Dict[str, object]] = None
        self.archive: Optional[Dict[str, object]] = None
        self.candidate_accept: Optional[Dict[str, object]] = None
        self.restored: Optional[Dict[str, object]] = None
        self.audit: Optional[Dict[str, object]] = None
        self.final_replay: Optional[Dict[str, object]] = None
        self.receipt_chain: Optional[Dict[str, object]] = None
        self.recovery_budget: Optional[EmergencyRecoveryBudget] = None
        self.emergency_recovery_completion: Optional[Dict[str, object]] = None
        self.recovery_action_evidence: Optional[Dict[str, object]] = None
        self.scoring_deadline_error: Optional[TransitionDeadline] = None
        self.scoring_checkpoint_completed_monotonic_s: Optional[float] = None
        self.scoring_evidence_completed_monotonic_s: Optional[float] = None

    def _check_signal(self) -> None:
        if self.signal_guard is not None:
            self.signal_guard.raise_if_requested()

    def _phase(self, name: str, action: Callable[[], Mapping[str, object]]) \
            -> Dict[str, object]:
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " action")
            result = dict(action())
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _best_effort_recovery_receipt(
            self, name: str, status_value: str,
            payload: Mapping[str, object]) -> None:
        try:
            self.journal.record(name, status_value, payload)
        except BaseException as exc:
            self.recovery_errors.append(exc)

    def _capture_recovery_action_evidence(
            self,
    ) -> Optional[Dict[str, object]]:
        capture = getattr(self.backend, "recovery_action_evidence", None)
        if capture is None:
            return self.recovery_action_evidence
        value = dict(capture())
        self.recovery_action_evidence = value
        return value

    def _safe_recovery_action_evidence(
            self,
    ) -> Optional[Dict[str, object]]:
        try:
            return self._capture_recovery_action_evidence()
        except BaseException as exc:
            self.recovery_errors.append(exc)
            return self.recovery_action_evidence

    def _record_scoring_deadline_exhaustion(
            self, observed_monotonic_s: Optional[float] = None) -> None:
        observed = self.deadline.now() if observed_monotonic_s is None \
            else observed_monotonic_s
        if not math.isfinite(observed) or observed < self.deadline.absolute or \
                self.scoring_deadline_error is not None:
            return
        error = TransitionDeadline(
            "scoring absolute deadline exhausted; safe recovery continued")
        self.scoring_deadline_error = error
        self.recovery_errors.append(error)
        self._best_effort_recovery_receipt(
            "scoring_deadline", "failed", {
                "error": error_record(error),
                "observed_monotonic_s": observed,
                "recovery_actions_continue": True,
            })

    def _stop_old_phase(self) -> Dict[str, object]:
        """Arm in memory only after the durable stop-start receipt exists."""
        name = "old_stop"
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " signal")
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        # From this assignment onward, every exception and deferred signal must
        # enter recovery.  Before it, the live sampler is provably untouched.
        self.recovery_armed = True
        try:
            result = dict(self.backend.stop_old())
        except BaseException as exc:
            action_evidence = self._safe_recovery_action_evidence()
            try:
                self.journal.record(name, "failed", {
                    "error": error_record(exc),
                    "process_signal_actions": action_evidence,
                })
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _candidate_accept_phase(self) -> Dict[str, object]:
        """Retain acceptance evidence even if its completion receipt fails."""
        name = "candidate_accept"
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " action")
            result = dict(self.backend.accept_candidate())
            self.candidate_accept = result
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _recover(self) -> None:
        budget = self.deadline.start_emergency_recovery(
            self.plan.emergency_recovery_s)
        self.recovery_budget = budget
        emergency_setup_error: Optional[BaseException] = None
        emergency_setup: Dict[str, object] = {}
        try:
            emergency_setup = dict(
                self.backend.begin_emergency_recovery(budget))
        except BaseException as exc:
            emergency_setup_error = exc
            self.recovery_errors.append(exc)
        self._best_effort_recovery_receipt(
            "emergency_recovery", "started", {
                **budget.receipt(), "backend": emergency_setup})
        initial_actions = self._safe_recovery_action_evidence()
        self._best_effort_recovery_receipt(
            "process_signal_actions", "completed", {
                "checkpoint": "recovery-started",
                "evidence": initial_actions,
            })
        self._record_scoring_deadline_exhaustion()

        self._best_effort_recovery_receipt(
            "candidate_cleanup", "started", {})
        try:
            cleanup = dict(self.backend.cleanup_candidate())
        except BaseException as exc:
            self.recovery_errors.append(exc)
            self._best_effort_recovery_receipt(
                "candidate_cleanup", "failed", {"error": error_record(exc)})
        else:
            self._best_effort_recovery_receipt(
                "candidate_cleanup", "completed", cleanup)

        self._best_effort_recovery_receipt("old_restore", "started", {})
        try:
            self.restored = dict(self.backend.restore_old(self.archive))
        except BaseException as exc:
            self.recovery_errors.append(exc)
            action_evidence = self._safe_recovery_action_evidence()
            self._best_effort_recovery_receipt(
                "old_restore", "failed", {
                    "error": error_record(exc),
                    "process_signal_actions": action_evidence,
                })
        else:
            action_evidence = self._safe_recovery_action_evidence()
            if self.archive is None and isinstance(
                    self.restored.get("archive_record"), dict):
                self.archive = dict(self.restored["archive_record"])
            self._best_effort_recovery_receipt(
                "old_restore", "completed", {
                    **self.restored,
                    "process_signal_actions": action_evidence,
                })

            self._best_effort_recovery_receipt(
                "audit_binding", "started", {})
            try:
                pre_audit_receipt_chain = self.journal.replay()
                self.audit = dict(self.backend.publish_audit_binding(
                    self.archive, self.restored, self.dry_accepted,
                    self.candidate_accept, pre_audit_receipt_chain))
            except BaseException as exc:
                self.recovery_errors.append(exc)
                self._best_effort_recovery_receipt(
                    "audit_binding", "failed", {"error": error_record(exc)})
            else:
                self._best_effort_recovery_receipt(
                    "audit_binding", "completed", self.audit)

                self._best_effort_recovery_receipt(
                    "final_replay", "started", {})
                try:
                    self.final_replay = dict(self.backend.final_replay(
                        self.candidate_accept, self.archive, self.restored,
                        self.audit))
                except BaseException as exc:
                    self.recovery_errors.append(exc)
                    self._best_effort_recovery_receipt(
                        "final_replay", "failed", {"error": error_record(exc)})
                else:
                    self._best_effort_recovery_receipt(
                        "final_replay", "completed", self.final_replay)

        self._record_scoring_deadline_exhaustion()
        final_actions = self._safe_recovery_action_evidence()
        self._best_effort_recovery_receipt(
            "process_signal_actions", "completed", {
                "checkpoint": "recovery-finished",
                "evidence": final_actions,
            })
        completion = budget.receipt()
        self.emergency_recovery_completion = completion
        emergency_error = emergency_setup_error
        if completion["exhausted"]:
            emergency_error = TransitionDeadline(
                "emergency recovery budget exhausted after safe actions were attempted")
            self.recovery_errors.append(emergency_error)
        self._best_effort_recovery_receipt(
            "emergency_recovery",
            "failed" if emergency_error is not None else "completed",
            {**completion,
             **({"error": error_record(emergency_error)}
                if emergency_error is not None else {})})

    def run(self) -> Dict[str, object]:
        try:
            self._phase("hard_preflight", self.backend.hard_preflight)
            self._phase("recovery_arm", self.backend.arm_recovery)
            self.old_stop = self._stop_old_phase()
            if self.old_stop.get("forced") is True:
                raise TransitionError(
                    "forced old stop vetoes candidate dry run")
            self.archive = self._phase("old_archive", self.backend.archive_old)
            if self.archive.get("forced_stop") is True:
                raise TransitionError(
                    "forced old stop vetoes candidate dry run")
            self._phase("candidate_start", self.backend.start_candidate)
            self._phase("candidate_exercise", self.backend.exercise_candidate)
            self._phase("candidate_stop", self.backend.stop_candidate)
            self._candidate_accept_phase()
            self.dry_accepted = True
        except BaseException as exc:
            self.primary_error = exc
        finally:
            if self.recovery_armed:
                try:
                    self._recover()
                except BaseException as exc:
                    self.recovery_errors.append(exc)

        # Signals received during recovery are deliberately held until cleanup,
        # old restoration, and audit publication have all had their chance.  Do
        # not silently turn such a request into a successful transition.
        if self.signal_guard is not None:
            try:
                self.signal_guard.raise_if_requested()
            except BaseException as exc:
                if self.primary_error is None:
                    self.primary_error = exc
                else:
                    self.recovery_errors.append(exc)

        # Include a scoring overrun that happened during the final recovery
        # receipt/signal replay.  This can only change acceptance; recovery has
        # already had every safe action attempted.
        if self.recovery_armed:
            self._record_scoring_deadline_exhaustion()
            # This independent backend snapshot survives every restore path,
            # including quiesce, archive, or restart failure.  It is not
            # conditional on a populated ``self.restored`` object.
            self._safe_recovery_action_evidence()

        # This is the only terminal journal write.  It deliberately contains
        # the complete evidence candidate but no outcome claim.  Only after the
        # durable write returns can the absolute scoring boundary be evaluated
        # without leaving a serialized success that a late write contradicts.
        try:
            self.journal.record("terminal", "started", {
                "archived_pre_dry": self.archive,
                "candidate_accept": self.candidate_accept,
                "dry_accepted": self.dry_accepted,
                "future_audit_binding": self.audit,
                "final_replay": self.final_replay,
                "old_sampler_restored": self.restored,
                "process_signal_actions": self.recovery_action_evidence,
                "primary_error": error_record(self.primary_error)
                    if self.primary_error is not None else None,
                "outcome_pending": True,
                "recovery_errors_before_checkpoint": [
                    error_record(exc) for exc in self.recovery_errors],
                "recovery_attempts_finished": self.recovery_armed,
                "transition_id": self.plan.transition_id,
                "transition_success": None,
            })
            self.scoring_checkpoint_completed_monotonic_s = self.deadline.now()
        except BaseException as exc:
            self.recovery_errors.append(exc)
        self._record_scoring_deadline_exhaustion(
            self.scoring_checkpoint_completed_monotonic_s)
        try:
            self.receipt_chain = self.journal.replay()
        except BaseException as exc:
            self.recovery_errors.append(exc)
        # Capture the final scoring observation exactly once.  A deliberately
        # slow replay is scored, while later terminal-object serialization can
        # neither contradict this observation nor reopen the recovery path.
        self.scoring_evidence_completed_monotonic_s = self.deadline.now()
        if self.scoring_evidence_completed_monotonic_s >= \
                self.deadline.absolute and \
                self.scoring_deadline_error is None:
            self._record_scoring_deadline_exhaustion(
                self.scoring_evidence_completed_monotonic_s)
            try:
                self.receipt_chain = self.journal.replay()
            except BaseException as exc:
                self.recovery_errors.append(exc)

        success = self.dry_accepted and self.restored is not None and \
            self.audit is not None and self.final_replay is not None and \
            self.receipt_chain is not None and \
            self.primary_error is None and not self.recovery_errors
        terminal = sealed_record(TRANSITION_SCHEMA, {
            "archived_pre_dry": self.archive,
            "candidate_accept": self.candidate_accept,
            "dry_accepted": self.dry_accepted,
            "finished_utc": utc_now(),
            "future_audit_binding": self.audit,
            "final_replay": self.final_replay,
            "old_sampler_restored": self.restored,
            "process_signal_actions": self.recovery_action_evidence,
            "primary_error": error_record(self.primary_error)
                if self.primary_error is not None else None,
            "recovery_armed": self.recovery_armed,
            "emergency_recovery": self.emergency_recovery_completion,
            "recovery_errors": [error_record(exc)
                                for exc in self.recovery_errors],
            "receipt_chain_prefix": self.receipt_chain,
            "scoring_deadline": {
                "absolute_monotonic_s": self.deadline.absolute,
                "checkpoint_completed_monotonic_s":
                    self.scoring_checkpoint_completed_monotonic_s,
                "evidence_completed_monotonic_s":
                    self.scoring_evidence_completed_monotonic_s,
                "evidence_replay_defines_scoring_boundary": True,
                "exhausted": self.scoring_evidence_completed_monotonic_s >=
                    self.deadline.absolute,
                "remaining_s": max(
                    0.0, self.deadline.absolute -
                    self.scoring_evidence_completed_monotonic_s),
                "terminal_checkpoint_is_outcome_neutral": True,
            },
            "success": success,
            "transition_id": self.plan.transition_id,
        })
        if not success:
            messages = []
            if self.primary_error is not None:
                messages.append("primary=" + str(self.primary_error))
            messages.extend("recovery=" + str(exc)
                            for exc in self.recovery_errors)
            raise TransitionError(
                "thermal transition failed: " + "; ".join(messages))
        return terminal


def _parse_proc_stat(raw: bytes) -> Dict[str, int | str]:
    right = raw.rfind(b")")
    fields = raw[right + 2:].split() if right >= 0 else []
    if len(fields) < 37:
        raise TransitionError("process stat is truncated")
    try:
        return {
            "state": fields[0].decode("ascii"),
            "ppid": int(fields[1]), "process_group": int(fields[2]),
            "session": int(fields[3]), "start_tick": int(fields[19]),
            "processor": int(fields[36]),
        }
    except (UnicodeDecodeError, ValueError) as exc:
        raise TransitionError("process stat is malformed") from exc


def capture_owned_session_leader(
        pid: int, *, proc_root: Path = Path("/proc"),
) -> int:
    if type(pid) is not int or pid <= 1:
        raise TransitionError("owned launcher PID is outside its domain")
    try:
        launcher = _parse_proc_stat(
            (proc_root / str(pid) / "stat").read_bytes())
    except FileNotFoundError:
        # sudo -b may have exited after leaving its exact background session.
        # The session cleanup helper treats zero as valid only while the leader
        # PID remains absent, and refuses a live/reused leader at this identity.
        return 0
    if launcher["process_group"] != pid or launcher["session"] != pid:
        raise TransitionError(
            "old restart launcher did not create an exact session")
    return int(launcher["start_tick"])


def _status_field(status: str, name: str) -> str:
    match = re.search(r"^" + re.escape(name) + r":\s*(\S(?:.*\S)?)\s*$",
                      status, re.MULTILINE)
    if match is None:
        raise TransitionError("process status lacks " + name)
    return match.group(1)


def parse_mpstat_idle_receipt(
        raw: bytes, cpu_pair: Sequence[int], *, expected_intervals: int = 3,
) -> Dict[int, float]:
    """Validate the exact finite per-CPU idle receipt requested from mpstat."""
    pair = set(cpu_pair)
    if len(pair) != 2 or any(type(cpu) is not int or cpu < 0 for cpu in pair) or \
            type(expected_intervals) is not int or expected_intervals <= 0:
        raise TransitionError("mpstat receipt contract is malformed")

    def reject_constant(value: str) -> object:
        raise ValueError("non-finite JSON constant: " + value)

    try:
        document = json.loads(
            raw.decode("ascii"), parse_constant=reject_constant)
        hosts = document["sysstat"]["hosts"]
        if not isinstance(hosts, list) or len(hosts) != 1:
            raise TransitionError("mpstat JSON has an unexpected host count")
        statistics = hosts[0]["statistics"]
        if not isinstance(statistics, list) or \
                len(statistics) != expected_intervals:
            raise TransitionError("mpstat JSON has an unexpected interval count")
        idle_by_cpu = {cpu: [] for cpu in pair}
        for interval in statistics:
            loads = interval["cpu-load"]
            if not isinstance(loads, list):
                raise TransitionError("mpstat cpu-load receipt is malformed")
            observed: Dict[int, float] = {}
            for load in loads:
                if not isinstance(load, dict):
                    raise TransitionError("mpstat CPU receipt is malformed")
                cpu_text = str(load["cpu"])
                if not cpu_text.isdecimal() or int(cpu_text) not in pair:
                    continue
                cpu = int(cpu_text)
                if cpu in observed:
                    raise TransitionError("mpstat CPU receipt is duplicated")
                idle_value = load["idle"]
                if type(idle_value) not in (int, float):
                    raise TransitionError("mpstat idle receipt is not numeric")
                idle = float(idle_value)
                if not math.isfinite(idle) or not 0.0 <= idle <= 100.0:
                    raise TransitionError("mpstat idle receipt is outside its domain")
                observed[cpu] = idle
            if set(observed) != pair:
                raise TransitionError("mpstat interval lacks the exact CPU pair")
            for cpu, idle in observed.items():
                idle_by_cpu[cpu].append(idle)
    except TransitionError:
        raise
    except (KeyError, TypeError, ValueError, UnicodeDecodeError,
            json.JSONDecodeError) as exc:
        raise TransitionError("mpstat JSON is malformed") from exc
    if any(len(values) != expected_intervals or
           sum(values) / expected_intervals < 99.0
           for values in idle_by_cpu.values()):
        raise TransitionError("candidate CPU pair is not at least 99% idle")
    return {
        cpu: sum(values) / expected_intervals
        for cpu, values in idle_by_cpu.items()
    }


def _proc_identity(pid: int) -> Dict[str, object]:
    if pid <= 1:
        raise TransitionError("process PID is outside its domain")
    root = Path("/proc") / str(pid)
    stat_raw = (root / "stat").read_bytes()
    cmdline_raw = (root / "cmdline").read_bytes()
    status = (root / "status").read_text(encoding="ascii")
    executable_info = os.stat(str(root / "exe"), follow_symlinks=True)
    if not cmdline_raw.endswith(b"\0"):
        raise TransitionError("process cmdline is not NUL terminated")
    try:
        argv = tuple(
            item.decode("ascii") for item in cmdline_raw[:-1].split(b"\0"))
    except UnicodeDecodeError as exc:
        raise TransitionError("process cmdline is not ASCII") from exc
    parsed = _parse_proc_stat(stat_raw)
    uid_text = _status_field(status, "Uid")
    try:
        uids = tuple(int(value) for value in uid_text.split())
    except ValueError as exc:
        raise TransitionError("process UID receipt is malformed") from exc
    if len(uids) != 4:
        raise TransitionError("process UID receipt is malformed")
    return {
        **parsed, "affinity": _status_field(status, "Cpus_allowed_list"),
        "argv": list(argv), "cmdline_sha256": sha256_bytes(cmdline_raw),
        "executable": {"device": executable_info.st_dev,
                       "inode": executable_info.st_ino},
        "pid": pid, "uids": list(uids),
    }


def _mkdir_exact(path: Path, mode: int, owner_uid: int,
                 owner_gid: Optional[int] = None) -> None:
    path.mkdir(mode=mode)
    os.chmod(path, mode)
    observed = directory_binding(path)
    if observed["uid"] != owner_uid or \
            (owner_gid is not None and observed["gid"] != owner_gid) or \
            observed["mode"] != mode:
        raise TransitionError("new directory binding is unsafe: %s" % path)
    fsync_directory(path.parent)


def _prepare_transition_source_test_model(
        plan: TransitionPlan, *, controller_source: Path,
        candidate_source: Path = REPO_CANDIDATE,
        p32_source: Path = REPO_P32) -> Dict[str, object]:
    """Exercise historical preparation in a synthetic, disjoint namespace.

    ``candidate_source`` and ``p32_source`` remain keyword-compatible with the
    v1 test harness, but any path other than the verified seal is rejected.
    The production V6 namespace is rejected before any filesystem operation.
    """
    _require_source_test_namespace(plan, "prepare source-test model")
    if os.geteuid() != plan.controller_uid:
        raise TransitionError("transition preparation controller UID changed")
    if Path(controller_source).resolve() != plan.controller.resolve() or \
            Path(candidate_source).resolve() != plan.sampler.resolve() or \
            Path(p32_source).resolve() != plan.p32.resolve() or \
            SHA256_RE.fullmatch(plan.controller_sha256) is None or \
            sha256_file(controller_source) != plan.controller_sha256:
        raise TransitionError(
            "preparation must execute entirely from the reviewed root code seal")
    try:
        Deadline(plan.deadline_s, plan.recovery_reserve_s)
        EmergencyRecoveryBudget(plan.emergency_recovery_s)
    except ValueError as exc:
        raise TransitionError("transition deadline plan is malformed") from exc
    code_seal = verify_code_seal(plan)
    tool_records = capture_tool_records()
    prepare_runtime = verify_prepare_runtime(plan, tool_records)
    if plan.root.exists() or plan.root.is_symlink():
        raise TransitionError("transition root already exists")
    _mkdir_exact(
        plan.root, 0o700, plan.controller_uid, plan.controller_gid)
    try:
        for directory in (
                plan.root / "segments", plan.root / "interrupted",
                plan.root / "runtime-home", plan.receipts):
            _mkdir_exact(
                directory, 0o700, plan.controller_uid, plan.controller_gid)
        code_stage_binding = write_new(
            plan.code_stage_receipt, canonical_json(code_seal),
            owner_uid=plan.controller_uid)
        code_files = code_seal["files"]
        plan_value = sealed_record(PLAN_SCHEMA, {
            "absolute_deadline_s": plan.deadline_s,
            "candidate": code_files["candidate"]["binding"],
            "candidate_cpu": plan.candidate_cpu,
            "candidate_sha256": plan.candidate_sha256,
            "candidate_sibling": plan.candidate_sibling,
            "code_owner_gid": plan.code_owner_gid,
            "code_owner_uid": plan.code_owner_uid,
            "code_seal": code_seal,
            "code_stage_receipt": {
                "binding": code_stage_binding,
                "path": str(plan.code_stage_receipt),
            },
            "controller": code_files["controller"]["binding"],
            "controller_gid": plan.controller_gid,
            "controller_sha256": plan.controller_sha256,
            "controller_uid": plan.controller_uid,
            "created_utc": utc_now(),
            "destinations": {
                "audit_binding": str(plan.audit_binding),
                "old_archive": str(plan.old_archive),
                "old_stale_pid_archive": str(plan.old_stale_pid_archive),
                "old_unclean_archive": str(plan.old_unclean_archive),
            },
            "emergency_recovery_s": plan.emergency_recovery_s,
            "old": {
                "argv": list(plan.old_argv),
                "boot_id": plan.old_boot_id,
                "cmdline_sha256": plan.old_cmdline_sha256,
                "cpu": plan.old_cpu, "pid": plan.old_pid,
                "csv": {
                    "device": plan.old_csv_device,
                    "inode": plan.old_csv_inode,
                    "path": str(plan.old_csv),
                },
                "launcher": {
                    "affinity": plan.old_launcher_affinity,
                    "cmdline_sha256": plan.old_launcher_cmdline_sha256,
                    "pid": plan.old_launcher_pid,
                    "start_tick": plan.old_launcher_start_tick,
                    "uids": list(plan.old_launcher_uids),
                },
                "process_group": plan.old_process_group,
                "pid_file": {
                    "device": plan.old_pid_device,
                    "inode": plan.old_pid_inode,
                    "path": str(plan.old_pid_file),
                },
                "session": plan.old_session,
                "source_path": str(plan.old_source),
                "source_sha256": plan.old_source_sha256,
                "start_tick": plan.old_start_tick,
                "replacement_argv": list(plan.replacement_old_argv),
                "replacement_cmdline_sha256":
                    plan.replacement_old_cmdline_sha256,
                "replacement_source":
                    code_files["legacy"]["binding"],
            },
            "p32": code_files["p32"]["binding"],
            "p32_split_custody": {
                "frozen_root": str(plan.code_seal_root),
                "only_routed_namespace": "frozen",
                "output_root": str(plan.root),
            },
            "p32_sha256": plan.p32_sha256,
            "prepare_runtime": prepare_runtime,
            "recovery_reserve_s": plan.recovery_reserve_s,
            "root": str(plan.root),
            "execution": {
                "enabled": False,
                "retirement": SOURCE_ONLY_RETIREMENT,
            },
            "tools": tool_records,
            "transition_id": plan.transition_id,
        })
        plan_binding = write_new(
            plan.plan_receipt, canonical_json(plan_value),
            owner_uid=plan.controller_uid)
        return {"plan": plan_binding, "value": plan_value}
    except BaseException:
        # Preparation happens before any live inspection or stop.  Preserve a
        # partial root for forensic inspection rather than recursively remove.
        raise


def prepare_transition(plan: TransitionPlan, *, controller_source: Path,
                       candidate_source: Path = REPO_CANDIDATE,
                       p32_source: Path = REPO_P32) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)


def verify_transition_plan(plan: TransitionPlan) -> Dict[str, object]:
    if os.geteuid() != plan.controller_uid:
        raise TransitionError("controller UID differs from the frozen plan")
    value = load_canonical(plan.plan_receipt, PLAN_SCHEMA, "transition plan")
    if value.get("root") != str(plan.root) or \
            value.get("transition_id") != plan.transition_id or \
            value.get("controller_sha256") != plan.controller_sha256 or \
            value.get("controller_uid") != plan.controller_uid or \
            value.get("controller_gid") != plan.controller_gid or \
            value.get("code_owner_uid") != plan.code_owner_uid or \
            value.get("code_owner_gid") != plan.code_owner_gid or \
            value.get("candidate_sha256") != plan.candidate_sha256 or \
            value.get("p32_sha256") != plan.p32_sha256 or \
            value.get("candidate_cpu") != plan.candidate_cpu or \
            value.get("candidate_sibling") != plan.candidate_sibling or \
            value.get("absolute_deadline_s") != plan.deadline_s or \
            value.get("recovery_reserve_s") != plan.recovery_reserve_s or \
            value.get("emergency_recovery_s") != plan.emergency_recovery_s:
        raise TransitionError("transition plan contract changed")
    try:
        Deadline(plan.deadline_s, plan.recovery_reserve_s)
        EmergencyRecoveryBudget(plan.emergency_recovery_s)
    except ValueError as exc:
        raise TransitionError("transition deadline plan is malformed") from exc
    code_seal = verify_code_seal(plan)
    if value.get("code_seal") != code_seal:
        raise TransitionError("root code-seal receipt changed")
    expected_old = {
        "argv": list(plan.old_argv),
        "boot_id": plan.old_boot_id,
        "cmdline_sha256": plan.old_cmdline_sha256,
        "cpu": plan.old_cpu, "pid": plan.old_pid,
        "csv": {
            "device": plan.old_csv_device,
            "inode": plan.old_csv_inode,
            "path": str(plan.old_csv),
        },
        "launcher": {
            "affinity": plan.old_launcher_affinity,
            "cmdline_sha256": plan.old_launcher_cmdline_sha256,
            "pid": plan.old_launcher_pid,
            "start_tick": plan.old_launcher_start_tick,
            "uids": list(plan.old_launcher_uids),
        },
        "process_group": plan.old_process_group,
        "pid_file": {
            "device": plan.old_pid_device,
            "inode": plan.old_pid_inode,
            "path": str(plan.old_pid_file),
        },
        "session": plan.old_session,
        "source_path": str(plan.old_source),
        "source_sha256": plan.old_source_sha256,
        "start_tick": plan.old_start_tick,
        "replacement_argv": list(plan.replacement_old_argv),
        "replacement_cmdline_sha256": plan.replacement_old_cmdline_sha256,
        "replacement_source":
            code_seal["files"]["legacy"]["binding"],
    }
    if value.get("old") != expected_old:
        raise TransitionError("frozen old-sampler plan changed")
    if value.get("destinations") != {
            "audit_binding": str(plan.audit_binding),
            "old_archive": str(plan.old_archive),
            "old_stale_pid_archive": str(plan.old_stale_pid_archive),
            "old_unclean_archive": str(plan.old_unclean_archive)}:
        raise TransitionError("transition destination plan changed")
    if value.get("p32_split_custody") != {
            "frozen_root": str(plan.code_seal_root),
            "only_routed_namespace": "frozen",
            "output_root": str(plan.root)}:
        raise TransitionError("P32 split-custody routing plan changed")
    for path, key, digest in (
            (plan.sampler, "candidate", plan.candidate_sha256),
            (plan.p32, "p32", plan.p32_sha256),
            (plan.controller, "controller", plan.controller_sha256)):
        binding = file_binding(path, with_hash=True)
        if binding != value.get(key) or \
                (digest is not None and binding["sha256"] != digest) or \
                binding["uid"] != plan.code_owner_uid or \
                binding["gid"] != plan.code_owner_gid or \
                binding["mode"] != 0o444 or binding["nlink"] != 1:
            raise TransitionError("root-sealed transition file changed: " + key)
    code_stage = load_canonical(
        plan.code_stage_receipt, CODE_SEAL_RECEIPT_SCHEMA,
        "root code-seal stage receipt")
    stage_binding = file_binding(plan.code_stage_receipt, with_hash=True)
    if code_stage != code_seal or value.get("code_stage_receipt") != {
            "binding": stage_binding, "path": str(plan.code_stage_receipt)} or \
            stage_binding["uid"] != plan.controller_uid or \
            stage_binding["gid"] != plan.controller_gid or \
            stage_binding["mode"] != 0o444 or stage_binding["nlink"] != 1:
        raise TransitionError("root code-seal stage receipt changed")
    for directory in (
            plan.root, plan.root / "segments", plan.root / "interrupted",
            plan.root / "runtime-home", plan.receipts):
        binding = directory_binding(directory)
        if binding["uid"] != plan.controller_uid or \
                binding["gid"] != plan.controller_gid or \
                binding["mode"] != 0o700:
            raise TransitionError("transition output parent changed")
    tools = verify_tool_records(value.get("tools"))
    if value.get("prepare_runtime") != \
            prepare_runtime_contract(plan, tools):
        raise TransitionError("frozen prepare runtime contract changed")
    if value.get("execution") != {
            "enabled": False, "retirement": SOURCE_ONLY_RETIREMENT} or \
            any(key in value for key in (
                "execute_contract", "identity_inspection",
                "python_isolation")):
        raise TransitionError("source-only execution retirement changed")
    if not plan.old_argv or \
            plan.old_argv[0] != tools["python3"]["path"] or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.old_argv) + b"\0") != \
                plan.old_cmdline_sha256 or \
            tuple(plan.replacement_old_argv[:4]) != (
                tools["python3"]["path"], "-I", "-S", "-B") or \
            tuple(plan.replacement_old_argv[4:5]) != \
                (str(plan.legacy_sampler),) or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.replacement_old_argv) + b"\0") != \
                plan.replacement_old_cmdline_sha256 or \
            len(plan.old_launch_argv) < 4 or \
            tuple(plan.old_launch_argv[:3]) != ("sudo", "-n", "-b") or \
            plan.old_launch_argv[3] != tools["taskset"]["path"] or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.old_launch_argv) + b"\0") != \
                plan.old_launcher_cmdline_sha256:
        raise TransitionError(
            "old sampler argv differs from the sealed executable ledger")
    return value


def load_verified_p32(path: Path, expected_sha256: str) -> ModuleType:
    binding = file_binding(path, with_hash=True)
    if binding["sha256"] != expected_sha256 or binding["mode"] != 0o444 or \
            binding["nlink"] != 1:
        raise TransitionError("P32 helper binding changed")
    name = "_wh2_transition_frozen_p32"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise TransitionError("cannot load frozen P32 helper")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    if sha256_file(path) != expected_sha256:
        raise TransitionError("P32 helper changed while loading")
    return module


def _special_file_record(path: Path) -> Dict[str, object]:
    try:
        before = path.lstat()
        if stat.S_ISLNK(before.st_mode):
            raise TransitionError("special path is a symlink: %s" % path)
        after = path.lstat()
    except OSError as exc:
        raise TransitionError("special path is unavailable: %s" % path) from exc
    fields = (before.st_dev, before.st_ino, before.st_mode, before.st_rdev,
              before.st_uid, before.st_gid, before.st_nlink)
    if fields != (after.st_dev, after.st_ino, after.st_mode, after.st_rdev,
                  after.st_uid, after.st_gid, after.st_nlink):
        raise TransitionError("special path changed while binding: %s" % path)
    kind = "char" if stat.S_ISCHR(before.st_mode) else \
        "regular" if stat.S_ISREG(before.st_mode) else "other"
    return {
        "device": before.st_dev, "gid": before.st_gid,
        "inode": before.st_ino, "kind": kind,
        "major": os.major(before.st_rdev) if stat.S_ISCHR(before.st_mode) else 0,
        "minor": os.minor(before.st_rdev) if stat.S_ISCHR(before.st_mode) else 0,
        "mode": stat.S_IMODE(before.st_mode), "nlink": before.st_nlink,
        "rdev": before.st_rdev, "size": before.st_size,
        "uid": before.st_uid,
    }


def _fd_flags(fdinfo: bytes) -> int:
    try:
        text = fdinfo.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TransitionError("target fdinfo is not ASCII") from exc
    matches = re.findall(r"^flags:\s*([0-7]+)\s*$", text, re.MULTILINE)
    if len(matches) != 1:
        raise TransitionError("target fdinfo lacks one canonical flags field")
    try:
        return int(matches[0], 8)
    except ValueError as exc:
        raise TransitionError("target fd flags are malformed") from exc


def _capture_fd_record(fd_path: Path, fd: int) -> Dict[str, object]:
    try:
        flags_before = _fd_flags(
            fd_path.parent.parent.joinpath("fdinfo", str(fd)).read_bytes())
        info_before = fd_path.stat()
        link = os.readlink(fd_path)
        info_after = fd_path.stat()
        flags_after = _fd_flags(
            fd_path.parent.parent.joinpath("fdinfo", str(fd)).read_bytes())
    except OSError as exc:
        raise TransitionError("target descriptor changed during inspection") from exc
    identity_before = (
        info_before.st_dev, info_before.st_ino, info_before.st_mode,
        info_before.st_rdev, info_before.st_uid, info_before.st_gid,
        info_before.st_nlink)
    identity_after = (
        info_after.st_dev, info_after.st_ino, info_after.st_mode,
        info_after.st_rdev, info_after.st_uid, info_after.st_gid,
        info_after.st_nlink)
    if identity_before != identity_after or flags_before != flags_after:
        raise TransitionError("target descriptor was reused during inspection")
    kind = "char" if stat.S_ISCHR(info_before.st_mode) else \
        "regular" if stat.S_ISREG(info_before.st_mode) else "other"
    return {
        "access_mode": flags_before & os.O_ACCMODE,
        "device": info_before.st_dev, "fd": fd, "flags": flags_before,
        "gid": info_before.st_gid, "inode": info_before.st_ino,
        "kind": kind, "link_path": link,
        "major": os.major(info_before.st_rdev)
            if stat.S_ISCHR(info_before.st_mode) else 0,
        "minor": os.minor(info_before.st_rdev)
            if stat.S_ISCHR(info_before.st_mode) else 0,
        "mode": stat.S_IMODE(info_before.st_mode),
        "nlink": info_before.st_nlink, "rdev": info_before.st_rdev,
        "size": info_before.st_size, "uid": info_before.st_uid,
    }


def _fd_roster(fd_dir: Path) -> Tuple[Path, ...]:
    try:
        names = tuple(sorted(fd_dir.iterdir(), key=lambda path: int(path.name)))
    except (OSError, ValueError) as exc:
        raise TransitionError("target descriptor roster is unreadable") from exc
    if any(not path.name.isdecimal() or str(int(path.name)) != path.name
           for path in names):
        raise TransitionError("target descriptor name is noncanonical")
    return names


def capture_sampler_fd_snapshot(
        pid: int, *, source_path: Path, csv_path: Path,
        proc_root: Path = Path("/proc"),
        boot_id_path: Path = Path("/proc/sys/kernel/random/boot_id"),
        i2c_paths: Sequence[Path] = (Path("/dev/i2c-1"),
                                    Path("/dev/i2c-2")),
) -> Dict[str, object]:
    """Capture one race-checked target process/FD observation."""
    if type(pid) is not int or pid <= 1:
        raise TransitionError("FD inspection PID is outside its domain")
    if len(i2c_paths) != 2 or len({str(path) for path in i2c_paths}) != 2:
        raise TransitionError("FD inspection I2C path contract is malformed")
    pidfd = -1
    try:
        pidfd_open = getattr(os, "pidfd_open", None)
        if pidfd_open is None:
            raise TransitionError("pidfd_open is unavailable")
        pidfd = pidfd_open(pid, 0)
        boot_before = boot_id_path.read_text(encoding="ascii").strip()
        identity_before = _proc_identity(pid) if proc_root == Path("/proc") else \
            _proc_identity_at(proc_root, pid)
        fd_dir = proc_root / str(pid) / "fd"
        names = _fd_roster(fd_dir)
        records = [_capture_fd_record(path, int(path.name)) for path in names]
        csv_record = {"binding": file_binding(csv_path, with_hash=False),
                      "path": str(csv_path)}
        i2c_records = {
            str(path): _special_file_record(path) for path in i2c_paths}
        source_record = {
            "binding": file_binding(source_path, with_hash=True),
            "path": str(source_path)}
        names_after = _fd_roster(fd_dir)
        identity_after = _proc_identity(pid) if proc_root == Path("/proc") else \
            _proc_identity_at(proc_root, pid)
        boot_after = boot_id_path.read_text(encoding="ascii").strip()
        if tuple(path.name for path in names) != \
                tuple(path.name for path in names_after):
            raise TransitionError(
                "target descriptor roster changed during inspection")
        if identity_before != identity_after or boot_before != boot_after:
            raise TransitionError("target process changed during FD inspection")
        return {
            "boot_id": boot_before,
            "csv": csv_record,
            "fds": records,
            "i2c": i2c_records,
            "identity": identity_before,
            "pid": pid,
            "source": source_record,
        }
    finally:
        if pidfd >= 0:
            os.close(pidfd)


def _proc_identity_at(proc_root: Path, pid: int) -> Dict[str, object]:
    """Synthetic-proc variant used only by deterministic unit tests."""
    root = proc_root / str(pid)
    stat_raw = (root / "stat").read_bytes()
    cmdline_raw = (root / "cmdline").read_bytes()
    if not cmdline_raw:
        raise ProcessSnapshotRace("process cmdline entered teardown")
    if not cmdline_raw.endswith(b"\0"):
        raise TransitionError("process cmdline is not NUL terminated")
    status = (root / "status").read_text(encoding="ascii")
    executable_info = (root / "exe").stat()
    try:
        argv = tuple(
            item.decode("ascii") for item in cmdline_raw[:-1].split(b"\0"))
    except UnicodeDecodeError as exc:
        raise TransitionError("process cmdline is not ASCII") from exc
    parsed = _parse_proc_stat(stat_raw)
    try:
        uids = tuple(int(value) for value in _status_field(status, "Uid").split())
    except ValueError as exc:
        raise TransitionError("process UID receipt is malformed") from exc
    if len(uids) != 4:
        raise TransitionError("process UID receipt is malformed")
    return {
        **parsed, "affinity": _status_field(status, "Cpus_allowed_list"),
        "argv": list(argv), "cmdline_sha256": sha256_bytes(cmdline_raw),
        "executable": {"device": executable_info.st_dev,
                       "inode": executable_info.st_ino},
        "pid": pid, "uids": list(uids),
    }


def _fd_flags_match(
        record: Mapping[str, object], expected_access: int, *,
        required_status: int = 0, allowed_status: int = 0,
) -> bool:
    flags = record.get("flags")
    return type(flags) is int and flags >= 0 and \
        type(record.get("access_mode")) is int and \
        record.get("access_mode") == expected_access and \
        flags & os.O_ACCMODE == expected_access and \
        flags & required_status == required_status and \
        flags & (SEMANTIC_FD_STATUS_MASK & ~allowed_status) == 0


def validate_sampler_fd_snapshot(
        snapshot: Mapping[str, object], *, kind: str,
        expected_pid: int, expected_start_tick: int, expected_boot_id: str,
        expected_argv: Sequence[str], expected_python: Mapping[str, object],
        expected_source_sha256: str,
) -> Dict[str, object]:
    """Validate candidate or legacy source/CSV/I2C descriptor custody."""
    if kind not in ("candidate", "legacy"):
        raise TransitionError("FD inspection kind is malformed")
    try:
        identity = snapshot["identity"]
        source = snapshot["source"]
        csv_value = snapshot["csv"]
        i2c = snapshot["i2c"]
        fds = snapshot["fds"]
    except (KeyError, TypeError) as exc:
        raise TransitionError("FD inspection snapshot is incomplete") from exc
    if not isinstance(identity, dict) or not isinstance(source, dict) or \
            not isinstance(csv_value, dict) or not isinstance(i2c, dict) or \
            not isinstance(fds, list):
        raise TransitionError("FD inspection snapshot is malformed")
    if len(expected_argv) < 5 or snapshot.get("pid") != expected_pid or \
            identity.get("pid") != expected_pid or \
            identity.get("start_tick") != expected_start_tick or \
            snapshot.get("boot_id") != expected_boot_id or \
            identity.get("argv") != list(expected_argv) or \
            identity.get("cmdline_sha256") != sha256_bytes(
                b"\0".join(value.encode("ascii") for value in expected_argv) +
                b"\0") or identity.get("uids") != [0, 0, 0, 0] or \
            identity.get("executable") != {
                "device": expected_python.get("device"),
                "inode": expected_python.get("inode")}:
        raise TransitionError("FD inspection process identity changed")
    source_binding = source.get("binding")
    if source.get("path") != str(expected_argv[4]) or \
            not isinstance(source_binding, dict) or \
            source_binding.get("sha256") != expected_source_sha256 or \
            source_binding.get("uid") != 0 or source_binding.get("gid") != 0 or \
            source_binding.get("mode") != 0o444 or \
            source_binding.get("nlink") != 1:
        raise TransitionError("FD inspection root-sealed source changed")
    csv_binding = csv_value.get("binding")
    if not isinstance(csv_binding, dict) or \
            csv_binding.get("uid") != 0 or csv_binding.get("mode") != 0o444 or \
            csv_binding.get("nlink") != 1 or csv_binding.get("size", 0) <= 0:
        raise TransitionError("FD inspection CSV path binding is unsafe")
    if len({record.get("fd") for record in fds
            if isinstance(record, dict)}) != len(fds):
        raise TransitionError("FD inspection descriptor roster is duplicated")
    source_matches = [record for record in fds if isinstance(record, dict) and
                      record.get("kind") == "regular" and
                      (record.get("device"), record.get("inode")) ==
                      (source_binding.get("device"),
                       source_binding.get("inode"))]
    regular_identity_keys = (
        "device", "gid", "inode", "mode", "nlink", "size", "uid")
    if kind == "candidate":
        if len(source_matches) != 1 or \
                source_matches[0].get("link_path") != source.get("path") or \
                not _fd_flags_match(
                    source_matches[0], os.O_RDONLY,
                    required_status=os.O_NONBLOCK,
                    allowed_status=os.O_NONBLOCK) or \
                any(source_matches[0].get(key) != source_binding.get(key)
                    for key in regular_identity_keys):
            raise TransitionError(
                "candidate did not retain one exact self-open source FD")
        source_claim = {
            "basis": "retained-self-open-root-seal-source-fd",
            "fd": source_matches[0], "retained_source_fd": True,
        }
    else:
        if source_matches:
            raise TransitionError(
                "legacy retained an unclaimed root-seal source FD")
        source_claim = {
            "basis": "root-custody-normal-path-and-exact-argv",
            "retained_source_fd": False,
        }
    csv_matches = [record for record in fds if isinstance(record, dict) and
                   record.get("kind") == "regular" and
                   (record.get("device"), record.get("inode")) ==
                   (csv_binding.get("device"), csv_binding.get("inode"))]
    expected_csv_access = os.O_RDWR if kind == "candidate" else os.O_WRONLY
    if len(csv_matches) != 1 or \
            csv_matches[0].get("link_path") != csv_value.get("path") or \
            not _fd_flags_match(csv_matches[0], expected_csv_access) or \
            any(csv_matches[0].get(key) != csv_binding.get(key)
                for key in regular_identity_keys):
        raise TransitionError("sampler did not retain one exact growing CSV FD")
    expected_i2c_paths = {"/dev/i2c-1", "/dev/i2c-2"}
    if set(i2c) != expected_i2c_paths:
        raise TransitionError("FD inspection I2C path roster changed")
    for path in expected_i2c_paths:
        binding = i2c[path]
        if not isinstance(binding, dict) or binding.get("kind") != "char" or \
                binding.get("major") != 89 or \
                binding.get("minor") != int(path.rsplit("-", 1)[1]):
            raise TransitionError("FD inspection I2C device binding changed")
    i2c_fds = [record for record in fds if isinstance(record, dict) and
               record.get("kind") == "char" and record.get("major") == 89]
    if len(i2c_fds) != 2:
        raise TransitionError("sampler has a missing, duplicate, or additional I2C FD")
    by_path: Dict[str, Dict[str, object]] = {}
    for record in i2c_fds:
        path = record.get("link_path")
        if path not in expected_i2c_paths or path in by_path or \
                not _fd_flags_match(record, os.O_RDWR):
            raise TransitionError("sampler I2C FD path or flags changed")
        binding = i2c[path]
        if any(record.get(key) != binding.get(key) for key in (
                "device", "gid", "inode", "major", "minor", "mode",
                "nlink", "rdev", "size", "uid")):
            raise TransitionError("sampler I2C FD/path device binding changed")
        by_path[path] = record
    if set(by_path) != expected_i2c_paths or \
            len({record["fd"] for record in by_path.values()}) != 2 or \
            len({record["rdev"] for record in by_path.values()}) != 2:
        raise TransitionError("sampler I2C descriptor ownership is ambiguous")
    return sealed_record(FD_INSPECTION_SCHEMA, {
        "boot_id": expected_boot_id,
        "csv": {**csv_value, "fd": csv_matches[0]},
        "i2c": {path: {"device": i2c[path], "fd": by_path[path]}
                for path in sorted(expected_i2c_paths)},
        "kind": kind,
        "loaded_source": {**source, **source_claim},
        "ownership_scope":
            "target-retained-fds-and-per-device-observations-only;"
            "no-continuous-host-wide-exclusivity-claim",
        "process_identity": identity,
    })


def inspect_sampler_fd_provenance(
        plan: TransitionPlan, *, kind: str, pid: int, start_tick: int,
        boot_id: str, argv: Sequence[str], csv_path: Path,
        snapshot: Optional[Mapping[str, object]] = None,
) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)
    # Historical root-only FD-inspection entry model; unreachable in V6.
    if os.geteuid() != 0:
        raise TransitionError("FD inspection helper is not root")
    verify_code_seal(plan)
    source_path = plan.sampler if kind == "candidate" else plan.legacy_sampler
    expected_digest = plan.candidate_sha256 if kind == "candidate" else \
        plan.old_source_sha256
    captured = capture_sampler_fd_snapshot(
        pid, source_path=source_path, csv_path=csv_path) \
        if snapshot is None else dict(snapshot)
    receipt = validate_sampler_fd_snapshot(
        captured, kind=kind, expected_pid=pid,
        expected_start_tick=start_tick, expected_boot_id=boot_id,
        expected_argv=argv,
        expected_python=capture_tool_records()["python3"]["binding"],
        expected_source_sha256=expected_digest)
    # ``verify_code_seal`` above is a precondition rather than an appended
    # field, so the returned record retains its own canonical self hash.
    return receipt


class LiveBackend:
    """Historical operations model; V6 has no path that instantiates it."""

    def __init__(self, plan: TransitionPlan, p32: ModuleType,
                 deadline: Deadline,
                 tool_records: Mapping[str, Mapping[str, object]],
                 execute_runtime: Mapping[str, object]) -> None:
        self.plan = plan
        self.p32 = p32
        self.deadline = deadline
        self.old_preflight: Optional[Dict[str, object]] = None
        self.stop_record: Optional[Dict[str, object]] = None
        self.archive_record: Optional[Dict[str, object]] = None
        self.candidate_owner: Optional[object] = None
        self.candidate_timing_start: Optional[float] = None
        self.candidate_benchmark_end: Optional[float] = None
        self.candidate_terminal: Optional[Dict[str, object]] = None
        self.candidate_fd_bootstrap: Optional[Dict[str, object]] = None
        self.candidate_fd_replay: Optional[Dict[str, object]] = None
        self.legacy_fd_bootstrap: Optional[Dict[str, object]] = None
        self.legacy_fd_replay: Optional[Dict[str, object]] = None
        self.original_identity_preflight: Optional[Dict[str, object]] = None
        self.original_identity_arm: Optional[Dict[str, object]] = None
        self.replacement_identity_bootstrap: Optional[Dict[str, object]] = None
        self.original_signal_receipts: list[Dict[str, object]] = []
        self.original_signal_failures: Dict[str, Dict[str, object]] = {}
        self.candidate_cleanup_complete = False
        self.recovery_budget: Optional[EmergencyRecoveryBudget] = None
        self.tools = verify_tool_records(tool_records)
        self.execute_runtime = dict(execute_runtime)
        self.p32_root = SplitCustodyRoot(plan.root, plan.code_seal_root)
        self.design: Dict[str, object] = {
            "controller_uid": plan.controller_uid,
            "immutable_files": {
                "frozen/wirehair_expo_thermal_sampler.py":
                    plan.candidate_sha256,
            },
            "thermal_core": plan.candidate_cpu,
            "tools": self.tools,
            "python_isolation": python_isolation_contract(plan),
            "transition_execute_runtime": self.execute_runtime,
        }

    def begin_emergency_recovery(
            self, budget: EmergencyRecoveryBudget) -> Mapping[str, object]:
        if self.recovery_budget is not None or \
                not isinstance(budget, EmergencyRecoveryBudget):
            raise TransitionError("emergency recovery budget was installed twice")
        self.recovery_budget = budget
        tools = self._verify_tools()
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed before emergency recovery")
        return {"controller_interpreter": verify_running_interpreter(tools),
                "execute_runtime": runtime,
                "tools": tools}

    def _recovery_wait_deadline(
            self, maximum_wait_s: float, *,
            minimum_safety_wait_s: float = 0.05) -> float:
        if self.recovery_budget is None:
            raise TransitionError("recovery action lacks its emergency budget")
        return self.recovery_budget.wait_deadline(
            maximum_wait_s,
            minimum_safety_wait_s=minimum_safety_wait_s)

    def _recovery_now(self) -> float:
        if self.recovery_budget is None:
            raise TransitionError("recovery action lacks its emergency budget")
        return self.recovery_budget.now()

    def _verify_tools(self) -> Dict[str, Dict[str, object]]:
        return verify_tool_records(self.tools)

    def _environment(self) -> Dict[str, str]:
        return self.p32.sanitized_environment(
            self.plan.root / "runtime-home", allocator=False)

    def _i2c_readers(self) -> Tuple[int, ...]:
        return self.p32.sole_i2c_readers(
            Path(str(self.tools["fuser"]["path"])),
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])))

    def _fuser(self, path: Path) -> Tuple[int, ...]:
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])),
            (str(self.tools["fuser"]["path"]), str(path)),
            self._environment())
        if rc == 1 and not stdout and not stderr:
            return ()
        if rc != 0 or re.fullmatch(rb"(?: +[1-9][0-9]*)+", stdout) is None:
            raise TransitionError("fuser result is not canonical for %s" % path)
        expected_label = re.escape(os.fsencode(str(path))) + rb":[ ]*\n"
        if re.fullmatch(expected_label, stderr) is None:
            raise TransitionError("fuser label is not canonical for %s" % path)
        return tuple(sorted(set(int(value) for value in stdout.split())))

    def _inspect_sampler_fds(
            self, kind: str, *, pid: int, start_tick: int, boot_id: str,
            argv: Sequence[str], csv_path: Path,
    ) -> Dict[str, object]:
        if kind not in ("candidate", "legacy") or pid <= 1 or start_tick <= 0:
            raise TransitionError("privileged FD inspection request is malformed")
        argv_hex = canonical_json(list(argv)).hex()
        helper_argv = (
            str(self.tools["env"]["path"]), "-i",
            "HOME=" + str(self.plan.root / "runtime-home"),
            "PATH=/usr/bin:/bin", "LC_ALL=C", "LANG=C", "TZ=UTC",
            "PYTHONDONTWRITEBYTECODE=1",
            str(self.tools["python3"]["path"]), "-I", "-S", "-B",
            str(self.plan.controller),
            "--inspect-sealed-sampler-fds", kind,
            "--target-pid", str(pid),
            "--target-start-tick", str(start_tick),
            "--target-boot-id", boot_id,
            "--target-csv", str(csv_path),
            "--target-argv-json-hex", argv_hex,
            "--expected-controller-sha256", self.plan.controller_sha256,
        )
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), helper_argv,
            self._environment(), stdout_limit=32768, stderr_limit=4096)
        if rc != 0 or stderr:
            raise TransitionError("privileged sampler FD inspection failed")
        try:
            value = json.loads(stdout.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise TransitionError(
                "privileged sampler FD receipt is not JSON") from exc
        if not isinstance(value, dict) or canonical_json(value) != stdout:
            raise TransitionError(
                "privileged sampler FD receipt is noncanonical")
        receipt = verify_sealed(
            value, FD_INSPECTION_SCHEMA, "privileged sampler FD receipt")
        process = receipt.get("process_identity")
        csv_value = receipt.get("csv")
        source = receipt.get("loaded_source")
        expected_source = self.plan.sampler if kind == "candidate" else \
            self.plan.legacy_sampler
        expected_digest = self.plan.candidate_sha256 if kind == "candidate" \
            else self.plan.old_source_sha256
        if receipt.get("kind") != kind or receipt.get("boot_id") != boot_id or \
                not isinstance(process, dict) or process.get("pid") != pid or \
                process.get("start_tick") != start_tick or \
                process.get("argv") != list(argv) or \
                not isinstance(csv_value, dict) or \
                csv_value.get("path") != str(csv_path) or \
                not isinstance(source, dict) or \
                source.get("path") != str(expected_source) or \
                not isinstance(source.get("binding"), dict) or \
                source["binding"].get("sha256") != expected_digest or \
                (kind == "candidate" and
                 source.get("retained_source_fd") is not True) or \
                (kind == "legacy" and
                 (source.get("retained_source_fd") is not False or
                  source.get("basis") !=
                    "root-custody-normal-path-and-exact-argv")):
            raise TransitionError(
                "privileged sampler FD receipt changed semantically")
        readers = self._i2c_readers()
        if readers != (pid,):
            raise TransitionError(
                "sampler was not the sole per-device I2C reader at FD receipt")
        return {
            "fd_receipt": receipt,
            "per_device_sole_reader_observation": list(readers),
            "scope": "recorded-observation-only; sampler-retained-own-fds; "
                     "no-continuous-host-wide-exclusivity-claim",
        }

    def _inspect_process_identities(
            self, *, profile: str, child_pid: int, child_start_tick: int,
            launcher_pid: int, launcher_start_tick: int,
            allowed_absence: str,
    ) -> Dict[str, object]:
        """Run one single-shot root helper and validate its exact receipt."""
        request = identity_inspection_request(
            self.plan, profile=profile, child_pid=child_pid,
            child_start_tick=child_start_tick, launcher_pid=launcher_pid,
            launcher_start_tick=launcher_start_tick,
            controller_pid=os.getpid(), allowed_absence=allowed_absence)
        helper_argv = identity_inspection_command(
            self.plan, request, self.tools)
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), helper_argv,
            self._environment(), stdout_limit=65536, stderr_limit=4096)
        if rc != 0 or stderr:
            raise TransitionError(
                "privileged process identity inspection failed")
        try:
            value = json.loads(stdout.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise TransitionError(
                "privileged process identity receipt is not JSON") from exc
        if not isinstance(value, dict) or canonical_json(value) != stdout:
            raise TransitionError(
                "privileged process identity receipt is noncanonical")
        # Deliberately no retry: any semantic mismatch binds this inspection
        # attempt to an unaccepted identity and must fail closed.
        return verify_identity_inspection_receipt(
            self.plan, value, request, self.tools)

    def _inspect_original_identities(
            self, allowed_absence: str = "none",
    ) -> Dict[str, object]:
        return self._inspect_process_identities(
            profile="original", child_pid=self.plan.old_pid,
            child_start_tick=self.plan.old_start_tick,
            launcher_pid=self.plan.old_launcher_pid,
            launcher_start_tick=self.plan.old_launcher_start_tick,
            allowed_absence=allowed_absence)

    def _inspect_replacement_identities(
            self, *, child_pid: int, child_start_tick: int,
            launcher_pid: int, launcher_start_tick: int,
            allowed_absence: str,
    ) -> Dict[str, object]:
        return self._inspect_process_identities(
            profile="replacement", child_pid=child_pid,
            child_start_tick=child_start_tick,
            launcher_pid=launcher_pid,
            launcher_start_tick=launcher_start_tick,
            allowed_absence=allowed_absence)

    def _signal_original_target(
            self, target: str, signum: int,
    ) -> Dict[str, object]:
        allowed_absence = {"child": "launcher", "launcher": "child"}.get(
            target)
        if allowed_absence is None:
            raise TransitionError("original signal target is malformed")
        failures = getattr(self, "original_signal_failures", None)
        if failures is None:
            failures = {}
            self.original_signal_failures = failures
        if not isinstance(failures, dict):
            raise TransitionError("original signal failure log changed")
        if target in failures:
            raise TransitionError(
                "prior %s signal helper failure forbids identity retry" %
                target)
        receipts = getattr(self, "original_signal_receipts", None)
        if receipts is None:
            receipts = []
            self.original_signal_receipts = receipts
        if type(receipts) is not list:
            raise TransitionError("original signal receipt log changed")
        identity_request = identity_inspection_request(
            self.plan, profile="original", child_pid=self.plan.old_pid,
            child_start_tick=self.plan.old_start_tick,
            launcher_pid=self.plan.old_launcher_pid,
            launcher_start_tick=self.plan.old_launcher_start_tick,
            controller_pid=os.getpid(), allowed_absence=allowed_absence)
        request = process_signal_request(
            self.plan, identity_request, target=target, signum=signum)
        command = process_signal_command(self.plan, request, self.tools)
        attempt_invoked = False
        try:
            attempt_invoked = True
            rc, stdout, stderr = self.p32.run_privileged_bounded(
                Path(str(self.tools["sudo"]["path"])),
                Path(str(self.tools["timeout"]["path"])), command,
                self._environment(), stdout_limit=65536, stderr_limit=4096)
            if rc != 0 or stderr:
                raise TransitionError(
                    "exact root-sealed pidfd signal helper failed")
            try:
                value = json.loads(stdout.decode("ascii"))
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                raise TransitionError(
                    "exact process signal receipt is not JSON") from exc
            if not isinstance(value, dict) or canonical_json(value) != stdout:
                raise TransitionError(
                    "exact process signal receipt is noncanonical")
            receipt = verify_process_signal_receipt(
                self.plan, value, request, self.tools)
        except BaseException as exc:
            if attempt_invoked:
                failures[target] = {
                    "error": error_record(exc),
                    "no_retry": True,
                    "signal": int(signum),
                    "target": target,
                }
            raise
        receipts.append(receipt)
        return receipt

    def recovery_action_evidence(self) -> Dict[str, object]:
        """Snapshot actions independently of successful old restoration."""
        receipts = getattr(self, "original_signal_receipts", None)
        failures = getattr(self, "original_signal_failures", None)
        if type(receipts) is not list or type(failures) is not dict:
            raise TransitionError("recovery process signal log is malformed")
        verified_receipts = []
        for value in receipts:
            if not isinstance(value, dict):
                raise TransitionError(
                    "recovery process signal log is malformed")
            try:
                request = value["request"]
                tools = value["tools"]
                code_seal = value["code_seal"]
            except (KeyError, TypeError) as exc:
                raise TransitionError(
                    "recovery process signal log is malformed") from exc
            verified_receipts.append(
                _verify_process_signal_receipt_with_code_seal(
                    self.plan, value, request, tools, code_seal))
        for target, record in failures.items():
            if target not in ("child", "launcher") or \
                    not isinstance(record, dict) or set(record) != {
                        "error", "no_retry", "signal", "target"} or \
                    record.get("target") != target or \
                    record.get("no_retry") is not True or \
                    record.get("signal") not in (
                        int(signal.SIGTERM), int(signal.SIGKILL)) or \
                    not isinstance(record.get("error"), dict) or \
                    set(record["error"]) != {"message", "type"} or \
                    not all(isinstance(value, str)
                            for value in record["error"].values()):
                raise TransitionError(
                    "recovery process signal failure log is malformed")
        # Canonical JSON round-tripping makes this a detached snapshot; later
        # appends cannot rewrite an earlier journal or terminal payload.
        detached_receipts = json.loads(
            canonical_json(verified_receipts).decode("ascii"))
        detached_failures = json.loads(canonical_json(failures).decode("ascii"))
        return sealed_record(RECOVERY_ACTION_EVIDENCE_SCHEMA, {
            "failed_attempts": detached_failures,
            "receipt_count": len(detached_receipts),
            "signal_receipts": detached_receipts,
            "transition_id": self.plan.transition_id,
        })

    @staticmethod
    def _replay_fd_growth(
            first: Mapping[str, object], second: Mapping[str, object],
            *, description: str) -> None:
        try:
            first_receipt = first["fd_receipt"]
            second_receipt = second["fd_receipt"]
            first_csv = first_receipt["csv"]
            second_csv = second_receipt["csv"]
            first_binding = first_csv["binding"]
            second_binding = second_csv["binding"]
            first_fd = first_csv["fd"]
            second_fd = second_csv["fd"]
        except (KeyError, TypeError) as exc:
            raise TransitionError(description + " FD replay is malformed") from exc
        identity_keys = ("device", "inode", "uid", "gid", "mode", "nlink")
        if any(first_binding.get(key) != second_binding.get(key)
               for key in identity_keys) or \
                any(first_fd.get(key) != second_fd.get(key)
                    for key in ("fd", "device", "inode", "access_mode",
                                "link_path")) or \
                type(first_binding.get("size")) is not int or \
                type(second_binding.get("size")) is not int or \
                second_binding["size"] <= first_binding["size"] or \
                second_fd.get("size") != second_binding["size"]:
            raise TransitionError(
                description + " CSV FD/path did not remain exact and grow")

    def _wait_for_legacy_csv_growth(
            self, bootstrap: Mapping[str, object]) -> None:
        try:
            initial_size = bootstrap["fd_receipt"]["csv"]["binding"]["size"]
        except (KeyError, TypeError) as exc:
            raise TransitionError(
                "legacy bootstrap CSV provenance is malformed") from exc
        if type(initial_size) is not int or initial_size <= 0:
            raise TransitionError(
                "legacy bootstrap CSV size is outside its domain")
        binding = file_binding(self.plan.old_csv, with_hash=False)
        if binding["size"] > initial_size:
            return
        deadline = self._recovery_wait_deadline(
            3.0, minimum_safety_wait_s=1.25)
        while self._recovery_now() < deadline:
            binding = file_binding(self.plan.old_csv, with_hash=False)
            if binding["size"] > initial_size:
                return
            time.sleep(0.05)
        raise TransitionError(
            "legacy sampler CSV did not grow before FD provenance replay")

    def _verify_topology_and_occupancy(self) -> Dict[str, object]:
        pair = {self.plan.candidate_cpu, self.plan.candidate_sibling}
        expected = "%d,%d" % tuple(sorted(pair))
        for cpu in sorted(pair):
            value = (Path("/sys/devices/system/cpu") / ("cpu%d" % cpu) /
                     "topology/thread_siblings_list").read_text(
                         encoding="ascii").strip()
            if value != expected:
                raise TransitionError("candidate CPU sibling topology changed")
        online = set(self.p32.parse_cpu_list(
            Path("/sys/devices/system/cpu/online").read_text(
                encoding="ascii").strip()))
        if not pair <= online:
            raise TransitionError("candidate CPU pair is offline")

        def scan() -> Tuple[list[Dict[str, object]], list[Dict[str, object]]]:
            pinned = []
            executing = []
            for task in Path("/proc").glob("[0-9]*/task/[0-9]*"):
                try:
                    pid = int(task.parent.parent.name)
                    tid = int(task.name)
                    cmdline = (Path("/proc") / str(pid) / "cmdline").read_bytes()
                    if not cmdline or pid == os.getpid():
                        continue
                    status = (task / "status").read_text(encoding="ascii")
                    allowed = set(self.p32.parse_cpu_list(
                        _status_field(status, "Cpus_allowed_list")))
                    proc_stat = _parse_proc_stat((task / "stat").read_bytes())
                except (OSError, TransitionError, ValueError) as exc:
                    if task.exists():
                        raise TransitionError(
                            "cannot classify a live task on the candidate CPU pair") \
                            from exc
                    continue
                record = {"pid": pid, "tid": tid,
                          "cmdline_sha256": sha256_bytes(cmdline)}
                if allowed and allowed <= pair:
                    pinned.append(record)
                if proc_stat["state"] == "R" and proc_stat["processor"] in pair:
                    executing.append(record)
            return pinned, executing

        pinned_before, executing_before = scan()
        command = [str(self.tools["mpstat"]["path"]), "-o", "JSON", "-P",
                   "%d,%d" % tuple(sorted(pair)), "1", "3"]
        rc, stdout, stderr = self.p32.run_bounded(
            command, self._environment(), 6.0, 128 * 1024, 4096)
        if rc != 0 or stderr:
            raise TransitionError("mpstat occupancy check failed")
        idle_average = parse_mpstat_idle_receipt(stdout, tuple(sorted(pair)))
        pinned_after, executing_after = scan()
        if pinned_before or executing_before or pinned_after or executing_after:
            raise TransitionError("candidate CPU pair has a user-space occupant")
        return {
            "cpu_pair": sorted(pair),
            "idle_average_pct": {
                str(cpu): average for cpu, average in idle_average.items()},
            "mpstat_sha256": sha256_bytes(stdout),
        }

    def _validate_old_child_identity(self) -> Dict[str, object]:
        receipt = self._inspect_original_identities("launcher")
        child = receipt["targets"]["child"]
        if child.get("state") != "present":
            raise TransitionError("old sampler identity is absent")
        return {"child": child["identity"], "boot_id": receipt["boot_id"],
                "identity_receipt": receipt}

    def _validate_old_identity(self) -> Dict[str, object]:
        receipt = self._inspect_original_identities("none")
        child = receipt["targets"]["child"]
        launcher = receipt["targets"]["launcher"]
        if child.get("state") != "present" or \
                launcher.get("state") != "present":
            raise TransitionError("old sampler/launcher identity is absent")
        if child["identity"].get("ppid") != self.plan.old_launcher_pid:
            raise TransitionError("old sampler launcher parent changed")
        return {"boot_id": receipt["boot_id"],
                "child": child["identity"],
                "launcher": launcher["identity"],
                "identity_receipt": receipt}

    def _require_exclusive_old_readers(self) -> Dict[str, object]:
        i2c_readers = self._i2c_readers()
        if i2c_readers != (self.plan.old_pid,):
            raise TransitionError("old sampler is not the sole I2C reader")
        csv_readers = self._fuser(self.plan.old_csv)
        if csv_readers != (self.plan.old_pid,):
            raise TransitionError(
                "fuzz/GFNI timing/auditor or unknown process is reading the live CSV")
        return {"csv_readers": list(csv_readers),
                "i2c_readers": list(i2c_readers)}

    def _validate_thermal_parent(
            self, expected: Optional[Mapping[str, object]] = None,
    ) -> Dict[str, int]:
        paths = (
            self.plan.old_csv, self.plan.old_pid_file, self.plan.old_archive,
            self.plan.old_unclean_archive, self.plan.old_stale_pid_archive,
            self.plan.audit_binding,
        )
        parent = self.plan.old_csv.parent
        if any(path.parent != parent for path in paths):
            raise TransitionError("thermal evidence destinations changed parent")
        observed = directory_binding(parent)
        if observed["uid"] != self.plan.controller_uid or \
                observed["mode"] != 0o700:
            raise TransitionError(
                "thermal evidence parent is outside the controller trust boundary")
        if expected is not None and observed != dict(expected):
            raise TransitionError("thermal evidence parent binding changed")
        return observed

    def hard_preflight(self) -> Mapping[str, object]:
        value = verify_transition_plan(self.plan)
        if value.get("tools") != self.tools:
            raise TransitionError("backend tool ledger differs from the frozen plan")
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed after entry")
        # Discover an unavailable libc/kernel interface before the old sampler
        # can be signalled.  The actual archive still uses RENAME_NOREPLACE and
        # verifies its exact post-rename inode/content binding.
        renameat2_function()
        thermal_parent = self._validate_thermal_parent()
        for path in (
                self.plan.old_archive, self.plan.old_unclean_archive,
                self.plan.old_stale_pid_archive, self.plan.audit_binding):
            if path.exists() or path.is_symlink():
                raise TransitionError("reserved transition destination exists: %s" % path)
        source = file_binding(self.plan.old_source, with_hash=True)
        if source["sha256"] != self.plan.old_source_sha256 or \
                source["uid"] != self.plan.controller_uid or \
                source["mode"] & 0o222 or source["nlink"] != 1:
            raise TransitionError("old sampler source binding changed")
        csv_binding = file_binding(self.plan.old_csv, with_hash=False)
        if (csv_binding["device"], csv_binding["inode"],
                csv_binding["uid"], csv_binding["mode"],
                csv_binding["nlink"]) != (
                self.plan.old_csv_device, self.plan.old_csv_inode,
                0, 0o444, 1):
            raise TransitionError("live old CSV binding changed")
        pid_binding = file_binding(self.plan.old_pid_file, with_hash=True)
        if (pid_binding["device"], pid_binding["inode"],
                pid_binding["uid"], pid_binding["mode"],
                pid_binding["nlink"], pid_binding["size"],
                pid_binding["sha256"]) != (
                self.plan.old_pid_device, self.plan.old_pid_inode,
                0, 0o444, 1, len(str(self.plan.old_pid)) + 1,
                sha256_bytes((str(self.plan.old_pid) + "\n").encode("ascii"))):
            raise TransitionError("live old PID binding changed")
        identity = self._validate_old_identity()
        self.original_identity_preflight = identity["identity_receipt"]
        topology = self._verify_topology_and_occupancy()
        reader_seal = self._require_exclusive_old_readers()
        csv_readers = tuple(reader_seal["csv_readers"])
        self.old_preflight = {
            "csv": csv_binding, "csv_readers": list(csv_readers),
            "identity": identity, "pid_file": pid_binding,
            "source": source, "thermal_parent": thermal_parent,
            "tools": tools, "controller_interpreter": interpreter,
            "execute_runtime": runtime,
            "topology": topology,
        }
        return self.old_preflight

    def arm_recovery(self) -> Mapping[str, object]:
        # The controller writes the completed recovery_arm phase receipt before
        # calling stop_old.  This method performs one last no-side-effect seal.
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed before stop arm")
        identity = self._validate_old_identity()
        self.original_identity_arm = identity["identity_receipt"]
        if self.old_preflight is None or not isinstance(
                self.old_preflight.get("thermal_parent"), dict):
            raise TransitionError("recovery arm lacks the thermal-parent seal")
        self._validate_thermal_parent(self.old_preflight["thermal_parent"])
        self._require_exclusive_old_readers()
        return {"armed_for_pid": self.plan.old_pid,
                "controller_interpreter": interpreter,
                "execute_runtime": runtime,
                "start_tick": identity["child"]["start_tick"],
                "identity_receipt": self.original_identity_arm,
                "tools": tools}

    def _old_pidfd_signal(self, signum: int) -> Dict[str, object]:
        return self._signal_original_target("child", signum)

    def _old_identity_lives(self) -> bool:
        receipt = self._inspect_original_identities("both")
        return receipt["targets"]["child"].get("state") == "present"

    def _old_launcher_identity_lives(self) -> bool:
        receipt = self._inspect_original_identities("both")
        return receipt["targets"]["launcher"].get("state") == "present"

    def _old_launcher_pidfd_signal(self, signum: int) -> Dict[str, object]:
        return self._signal_original_target("launcher", signum)

    def _identity_exit_retry_window(self) -> Tuple[Callable[[], float], float]:
        budget = getattr(self, "recovery_budget", None)
        if budget is not None:
            return self._recovery_now, self._recovery_wait_deadline(
                5.0, minimum_safety_wait_s=1.0)
        now = self.deadline.now
        return now, min(self.deadline.normal, now() + 5.0)

    def _require_original_target_absent(
            self, target: str,
    ) -> Dict[str, object]:
        if target not in ("child", "launcher"):
            raise TransitionError("old absence target is malformed")
        now, wait_deadline = self._identity_exit_retry_window()
        exiting_receipts = []
        while True:
            receipt = self._inspect_original_identities("both")
            state = receipt["targets"][target].get("state")
            if state == "absent":
                return {
                    "absence_receipt": receipt,
                    "exiting_receipts": exiting_receipts,
                    "target": target,
                }
            if state != "exiting":
                raise TransitionError(
                    "original %s remains with a non-exit identity" % target)
            exiting_receipts.append(receipt)
            if now() >= wait_deadline:
                raise TransitionError(
                    "original %s exit race did not reach pidfd-open absence" %
                    target)
            time.sleep(0.05)

    def _require_old_pid_absent(self) -> Dict[str, object]:
        return self._require_original_target_absent("child")

    def _require_old_launcher_absent(self) -> Dict[str, object]:
        return self._require_original_target_absent("launcher")

    def stop_old(self) -> Mapping[str, object]:
        if self.old_preflight is None:
            raise TransitionError("old stop lacks its preflight seal")
        tools_at_stop = self._verify_tools()
        interpreter_at_stop = verify_running_interpreter(tools_at_stop)
        self._require_exclusive_old_readers()
        signal_receipts = getattr(self, "original_signal_receipts", None)
        if signal_receipts is None:
            signal_receipts = []
            self.original_signal_receipts = signal_receipts
        signal_start = len(signal_receipts)
        self._old_pidfd_signal(signal.SIGTERM)
        forced = False
        graceful_deadline = min(self.deadline.absolute - 10.0,
                                time.monotonic() + 15.0)
        while self._old_identity_lives() and time.monotonic() < graceful_deadline:
            time.sleep(0.05)
        if self._old_identity_lives():
            forced = True
            self._old_pidfd_signal(signal.SIGKILL)
            kill_deadline = min(self.deadline.absolute - 5.0,
                                time.monotonic() + 5.0)
            while self._old_identity_lives() and time.monotonic() < kill_deadline:
                time.sleep(0.05)
        if self._old_identity_lives():
            raise TransitionError("exact old sampler survived bounded pidfd stop")
        child_absence = self._require_old_pid_absent()
        parent_deadline = min(self.deadline.absolute - 5.0,
                              time.monotonic() + 5.0)
        while self._old_launcher_identity_lives() and \
                time.monotonic() < parent_deadline:
            time.sleep(0.05)
        if self._old_launcher_identity_lives():
            raise TransitionError("old sudo launcher survived child shutdown")
        launcher_absence = self._require_old_launcher_absent()
        if self._i2c_readers():
            raise TransitionError("I2C reader remained after old sampler stop")
        if not forced and (self.plan.old_pid_file.exists() or
                           self.plan.old_pid_file.is_symlink()):
            raise TransitionError("graceful old stop retained its PID file")
        self.stop_record = {
            "absence_receipts": {
                "child": child_absence, "launcher": launcher_absence},
            "forced": forced, "i2c_readers_after": [],
            "controller_interpreter": interpreter_at_stop,
            "old_pid": self.plan.old_pid,
            "signal_receipts": signal_receipts[signal_start:],
            "tools": tools_at_stop,
        }
        return self.stop_record

    def archive_old(self) -> Mapping[str, object]:
        if self.stop_record is None:
            raise TransitionError("old archive lacks a completed stop")
        if self._i2c_readers() or self._fuser(self.plan.old_csv):
            raise TransitionError("old evidence remains open at archive")
        binding = file_binding(self.plan.old_csv, with_hash=True)
        if self._fuser(self.plan.old_csv):
            raise TransitionError("old CSV gained a reader while being archived")
        preflight_csv = self.old_preflight["csv"] if self.old_preflight else None
        if not isinstance(preflight_csv, dict) or any(
                binding[key] != preflight_csv[key]
                for key in ("device", "inode", "uid", "gid", "mode", "nlink")):
            raise TransitionError("old CSV inode changed before archive")
        destination = self.plan.old_unclean_archive \
            if self.stop_record["forced"] else self.plan.old_archive
        rename_noreplace(
            self.plan.old_csv, destination, binding,
            parent_uid=self.plan.controller_uid)
        stale_pid = None
        if self.plan.old_pid_file.exists() or self.plan.old_pid_file.is_symlink():
            if not self.stop_record["forced"]:
                raise TransitionError("graceful stop unexpectedly retained PID evidence")
            stale_binding = file_binding(self.plan.old_pid_file, with_hash=True)
            rename_noreplace(
                self.plan.old_pid_file, self.plan.old_stale_pid_archive,
                stale_binding, parent_uid=self.plan.controller_uid)
            stale_pid = {"binding": stale_binding,
                         "path": str(self.plan.old_stale_pid_archive)}
        self.archive_record = {
            "binding": binding, "forced_stop": self.stop_record["forced"],
            "path": str(destination), "stale_pid": stale_pid,
        }
        return self.archive_record

    def start_candidate(self) -> Mapping[str, object]:
        if self.stop_record is None or self.stop_record.get("forced") is not False:
            raise TransitionError("candidate start is vetoed after forced old stop")
        if self.archive_record is None or \
                self.archive_record.get("forced_stop") is not False or \
                self._i2c_readers():
            raise TransitionError("candidate start lacks an empty-I2C archive seal")
        launch_topology = self._verify_topology_and_occupancy()
        if self._i2c_readers():
            raise TransitionError("I2C ownership appeared during CPU recheck")
        self.candidate_owner = self.p32.start_sampler(
            self.p32_root, self.design, 0)
        self.candidate_timing_start = time.monotonic()
        self.candidate_fd_bootstrap = self._inspect_sampler_fds(
            "candidate", pid=self.candidate_owner.pid,
            start_tick=self.candidate_owner.identity["start_tick"],
            boot_id=self.candidate_owner.identity["boot_id"],
            argv=self.candidate_owner.identity["cmdline"],
            csv_path=self.candidate_owner.csv_part)
        return {
            "argv": list(self.candidate_owner.identity["cmdline"]),
            "pid": self.candidate_owner.pid,
            "start_tick": self.candidate_owner.identity["start_tick"],
            "cmdline_sha256": self.candidate_owner.identity["cmdline_sha256"],
            "expected_output_owner_uid": self.plan.controller_uid,
            "expected_source_sha256": self.plan.candidate_sha256,
            "fd_provenance": self.candidate_fd_bootstrap,
            "launch_topology": launch_topology,
        }

    def exercise_candidate(self) -> Mapping[str, object]:
        if self.candidate_owner is None or self.candidate_timing_start is None:
            raise TransitionError("candidate exercise lacks an owner")
        paths = self.p32._sampler_evidence_paths(
            self.plan.root, 0, final=False)
        last_count = 0
        exercise_deadline = min(self.deadline.normal, time.monotonic() + 12.0)
        while time.monotonic() < exercise_deadline:
            if not self.p32.process_identity_matches(
                    self.candidate_owner.identity, self.plan.candidate_cpu,
                    self.candidate_owner.csv_part):
                raise TransitionError("candidate identity changed during dry exercise")
            rows, samples = self.p32.stable_validated_thermal_rows(
                paths["csv"], paths["validation"], self.plan.candidate_sha256,
                self.plan.controller_uid, attempts=3)
            last_count = len(rows)
            if last_count >= 5:
                self.candidate_benchmark_end = time.monotonic()
                if self.candidate_fd_bootstrap is None:
                    raise TransitionError(
                        "candidate exercise lacks bootstrap FD provenance")
                self.candidate_fd_replay = self._inspect_sampler_fds(
                    "candidate", pid=self.candidate_owner.pid,
                    start_tick=self.candidate_owner.identity["start_tick"],
                    boot_id=self.candidate_owner.identity["boot_id"],
                    argv=self.candidate_owner.identity["cmdline"],
                    csv_path=self.candidate_owner.csv_part)
                self._replay_fd_growth(
                    self.candidate_fd_bootstrap, self.candidate_fd_replay,
                    description="candidate")
                return {"paired_rows_before_stop": last_count,
                        "last_decision": samples[-1]["decision"],
                        "fd_provenance_replay": self.candidate_fd_replay}
            time.sleep(0.05)
        raise TransitionDeadline(
            "candidate did not produce five paired rows: %d" % last_count)

    def stop_candidate(self) -> Mapping[str, object]:
        if self.candidate_owner is None or \
                self.candidate_timing_start is None or \
                self.candidate_benchmark_end is None:
            raise TransitionError("candidate stop lacks its exercise interval")
        self.candidate_terminal = self.p32.stop_sampler(
            self.candidate_owner, self.p32_root, self.design,
            self.candidate_timing_start, self.candidate_benchmark_end)
        self.candidate_cleanup_complete = True
        if self._i2c_readers():
            raise TransitionError("candidate stop retained an I2C reader")
        return self.candidate_terminal

    def accept_candidate(self) -> Mapping[str, object]:
        if self.candidate_terminal is None:
            raise TransitionError("candidate acceptance lacks terminal evidence")
        if self.candidate_fd_bootstrap is None or \
                self.candidate_fd_replay is None:
            raise TransitionError(
                "candidate acceptance lacks replayed FD provenance")
        evidence, validation_raw = self.p32.validate_sampler_terminal_evidence(
            self.p32_root, self.design, 0, final=True)
        paths = self.p32._sampler_evidence_paths(self.plan.root, 0, final=True)
        segment_entries = self._candidate_artifact_roster(paths)
        raw_csv = stable_file_bytes(paths["csv"])
        rows = self.p32._parse_thermal_csv(
            raw_csv, allow_raw_dimm_faults=True)
        samples = self.p32._parse_thermal_validation(
            validation_raw, rows, self.plan.candidate_sha256,
            self.plan.controller_uid)
        if self.candidate_timing_start is None or \
                self.candidate_benchmark_end is None:
            raise TransitionError(
                "candidate acceptance lacks its measured interval")
        thermal = self.p32.validate_thermal_interval(
            raw_csv, self.candidate_timing_start, self.candidate_benchmark_end,
            validation_raw=validation_raw,
            expected_source_sha256=self.plan.candidate_sha256,
            expected_output_owner_uid=self.plan.controller_uid)
        if thermal != self.candidate_terminal.get("thermal_summary"):
            raise TransitionError(
                "candidate terminal thermal summary did not replay")
        summary = evidence["sampler_summary"]
        if len(rows) < 6 or any(sample["decision"] != "continue"
                                for sample in samples) or \
                any(sensor["reason"] != "ok" or not sensor["valid"] or
                    sensor["hot"] or sensor["attempt_errors"] != 0
                    for sample in samples
                    for sensor in sample["sensors"].values()) or \
                summary["dimm_invalid_samples_total"] != 0 or \
                summary["dimm_read_error_samples_total"] != 0 or \
                summary["dimm_attempt_errors_total"] != 0 or \
                summary["edac_ce_delta"] != 0 or summary["edac_ue_delta"] != 0:
            raise TransitionError("candidate dry samples are not strictly clean")
        if thermal["cpu_max_c"] >= 85.0 or \
                max(thermal["dimm_max_c"].values()) >= 84.0 or \
                thermal["max_gap_s"] > 2.25:
            raise TransitionError("candidate dry environmental gate failed")
        artifact_bindings = {
            name: {"binding": file_binding(path, with_hash=True),
                   "path": str(path)}
            for name, path in paths.items()
        }
        expected_hashes = {
            "csv": self.candidate_terminal["thermal_csv_sha256"],
            "receipt": self.candidate_terminal["thermal_sampler_evidence"]
                ["terminal_receipt"]["sha256"],
            "validation": evidence["validation_jsonl"]["sha256"],
        }
        if any(record["binding"]["sha256"] != expected_hashes[name] or
               record["binding"]["uid"] != 0 or
               record["binding"]["mode"] != 0o444 or
               record["binding"]["nlink"] != 1
               for name, record in artifact_bindings.items()):
            raise TransitionError("candidate artifact binding changed at acceptance")
        source_binding = file_binding(self.plan.sampler, with_hash=True)
        if source_binding["sha256"] != self.plan.candidate_sha256 or \
                source_binding["uid"] != self.plan.code_owner_uid or \
                source_binding["gid"] != self.plan.code_owner_gid or \
                source_binding["mode"] != 0o444 or \
                source_binding["nlink"] != 1:
            raise TransitionError("candidate source binding changed at acceptance")
        code_seal = verify_code_seal(self.plan)
        return {
            "artifact_bindings": artifact_bindings,
            "artifact_roster": [str(path) for path in segment_entries],
            "candidate_receipt_sha256":
                self.candidate_terminal["thermal_sampler_evidence"]
                    ["terminal_receipt"]["sha256"],
            "raw_sha256": self.candidate_terminal["thermal_csv_sha256"],
            "sample_count": len(rows),
            "fd_provenance": {
                "bootstrap": self.candidate_fd_bootstrap,
                "replay": self.candidate_fd_replay,
            },
            "root_code_seal_sha256":
                code_seal["self_sha256_excluding_field"],
            "source": {"binding": source_binding,
                       "path": str(self.plan.sampler)},
            "validation_sha256": evidence["validation_jsonl"]["sha256"],
        }

    def _candidate_artifact_roster(
            self, paths: Mapping[str, Path]) -> list[Path]:
        if set(paths) != {"csv", "receipt", "validation"} or \
                any(not isinstance(path, Path) for path in paths.values()):
            raise TransitionError("candidate final artifact paths are malformed")
        try:
            segment_entries = sorted(self.plan.root.joinpath("segments").iterdir())
        except OSError as exc:
            raise TransitionError("candidate final artifact roster is unreadable") \
                from exc
        expected_paths = set(paths.values())
        if set(segment_entries) != expected_paths or any(
                path.is_symlink() or not path.is_file()
                for path in segment_entries):
            raise TransitionError("candidate final artifact roster changed")
        return segment_entries

    def _owned_session_members(
            self, session_id: int, *, proc_root: Path = Path("/proc"),
    ) -> Tuple[int, ...]:
        if type(session_id) is not int or session_id <= 1:
            raise TransitionError("candidate session identity is malformed")
        members = []
        for process_path in proc_root.glob("[0-9]*"):
            try:
                pid = int(process_path.name)
                identity = _parse_proc_stat(
                    (process_path / "stat").read_bytes())
            except FileNotFoundError:
                continue
            except (OSError, TransitionError, ValueError) as exc:
                if process_path.exists():
                    raise TransitionError(
                        "cannot classify a live process during candidate cleanup") \
                        from exc
                continue
            if identity["session"] == session_id or \
                    identity["process_group"] == session_id:
                members.append(pid)
        return tuple(sorted(set(members)))

    def _cleanup_proof_pause(self) -> None:
        # This separation is an evidence requirement, not a recovery wait.
        # Even an exhausted emergency budget must not collapse two probes into
        # the same instant and manufacture a vacuous stable-empty proof.
        time.sleep(0.05)

    @staticmethod
    def _candidate_pid_present(pid: int, *,
                               proc_root: Path = Path("/proc")) -> bool:
        if type(pid) is not int or pid <= 1:
            raise TransitionError("candidate PID identity is malformed")
        try:
            (proc_root / str(pid)).stat()
        except FileNotFoundError:
            return False
        except OSError as exc:
            raise TransitionError(
                "cannot prove candidate PID absence during cleanup") from exc
        return True

    def _candidate_cleanup_observation(
            self, owner: object) -> Dict[str, object]:
        returncode = getattr(owner.process, "poll")()
        identity_live = self.p32.process_identity_matches(
            owner.identity, self.plan.candidate_cpu, owner.csv_part)
        session_members = self._owned_session_members(owner.process.pid)
        i2c_readers = self._i2c_readers()
        final_paths = self.p32._sampler_evidence_paths(
            self.plan.root, owner.segment, final=True)
        final_csv = final_paths.get("csv")
        if not isinstance(final_csv, Path):
            raise TransitionError("candidate final CSV path is malformed")
        csv_readers = set()
        for path in {owner.csv_part, final_csv}:
            try:
                info = path.lstat()
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(info.st_mode):
                raise TransitionError(
                    "candidate cleanup CSV path is not a regular file")
            csv_readers.update(self._fuser(path))
        return {
            "candidate_pid_present": self._candidate_pid_present(owner.pid),
            "csv_readers": sorted(csv_readers),
            "exact_identity_live": identity_live,
            "i2c_readers": list(i2c_readers),
            "launcher_returncode": returncode,
            "session_members": list(session_members),
        }

    @staticmethod
    def _cleanup_observation_is_empty(
            observation: Mapping[str, object]) -> bool:
        return observation.get("launcher_returncode") is not None and \
            observation.get("candidate_pid_present") is False and \
            observation.get("exact_identity_live") is False and \
            observation.get("session_members") == [] and \
            observation.get("i2c_readers") == [] and \
            observation.get("csv_readers") == []

    def cleanup_candidate(self) -> Mapping[str, object]:
        owner = self.candidate_owner
        self.candidate_cleanup_complete = False
        kill_error: Optional[BaseException] = None
        if owner is not None:
            try:
                # The privileged launcher may already have exited while the
                # exact root-owned session remains alive.  Always direct the
                # identity-verifying helper at the sealed owned session.
                self.p32._kill_owned_sampler_session(
                    owner, self.plan.root, self.design)
            except BaseException as exc:
                # The helper can report a transport/reap error after it has
                # already cleared the exact session.  Do not infer either
                # success or failure: independently prove stable emptiness.
                kill_error = exc
            try:
                first = self._candidate_cleanup_observation(owner)
                self._cleanup_proof_pause()
                second = self._candidate_cleanup_observation(owner)
                if not self._cleanup_observation_is_empty(first) or \
                        not self._cleanup_observation_is_empty(second):
                    raise TransitionError(
                        "candidate cleanup could not prove stable empty ownership: "
                        "%r then %r" % (first, second))
            except BaseException as proof_error:
                if kill_error is not None:
                    raise TransitionError(
                        "candidate kill helper failed and independent empty-"
                        "ownership proof failed: helper=%s; proof=%s" %
                        (kill_error, proof_error)) from proof_error
                raise
            self.candidate_cleanup_complete = True
            if kill_error is not None:
                raise TransitionError(
                    "candidate kill helper failed after independently proven "
                    "empty ownership: %s" % kill_error) from kill_error
            return {
                "empty_ownership_observations": [first, second],
                "i2c_readers_after": [],
                "owner_created": True,
            }

        first_readers = self._i2c_readers()
        self._cleanup_proof_pause()
        second_readers = self._i2c_readers()
        allowed = ((), (self.plan.old_pid,))
        if first_readers not in allowed or second_readers not in allowed:
            raise TransitionError(
                "candidate cleanup found an unowned I2C reader: "
                "%r then %r" % (first_readers, second_readers))
        self.candidate_cleanup_complete = True
        return {"i2c_readers_after": list(second_readers),
                "owner_created": False,
                "stable_observations": 2}

    def _quiesce_old_for_recovery(self) -> Dict[str, object]:
        signal_receipts = getattr(self, "original_signal_receipts", None)
        if signal_receipts is None:
            signal_receipts = []
            self.original_signal_receipts = signal_receipts
        signal_count_before_recovery = len(signal_receipts)
        if self._old_identity_lives():
            # A stop helper can fail after delivering SIGTERM.  Complete the
            # same exact-identity stop before any restart rather than assume the
            # old process is healthy or launch a second reader.
            try:
                self._old_pidfd_signal(signal.SIGTERM)
            except BaseException:
                pass
            wait_deadline = self._recovery_wait_deadline(
                5.0, minimum_safety_wait_s=1.0)
            while self._old_identity_lives() and \
                    self._recovery_now() < wait_deadline:
                time.sleep(0.05)
            if self._old_identity_lives():
                self._old_pidfd_signal(signal.SIGKILL)
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._old_identity_lives() and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        if self._old_identity_lives():
            raise TransitionError("old sampler could not be quiesced for recovery")
        child_absence = self._require_old_pid_absent()
        if self._old_launcher_identity_lives():
            try:
                self._old_launcher_pidfd_signal(signal.SIGTERM)
            except BaseException:
                pass
            wait_deadline = self._recovery_wait_deadline(
                5.0, minimum_safety_wait_s=1.0)
            while self._old_launcher_identity_lives() and \
                    self._recovery_now() < wait_deadline:
                time.sleep(0.05)
            if self._old_launcher_identity_lives():
                self._old_launcher_pidfd_signal(signal.SIGKILL)
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._old_launcher_identity_lives() and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        launcher_absence = self._require_old_launcher_absent()
        readers = self._i2c_readers()
        if readers:
            raise TransitionError("I2C is not empty before old recovery")
        return {
            "absence_receipts": {
                "child": child_absence, "launcher": launcher_absence},
            "i2c_readers_after": [],
            # Keep the complete cumulative log.  A successful stop signal can
            # precede a later stop failure or failed phase-completion journal
            # write; recovery must carry that otherwise unjournaled action
            # into restored/final audit evidence.
            "signal_receipts": list(signal_receipts),
            "signal_receipts_before_recovery":
                signal_count_before_recovery,
            "signal_attempt_failures": dict(getattr(
                self, "original_signal_failures", {})),
        }

    def _ensure_recovery_archive(
            self, archive: Optional[Mapping[str, object]]) -> Dict[str, object]:
        def validate_csv_binding(binding: Mapping[str, object]) -> None:
            preflight = getattr(self, "old_preflight", None)
            expected = preflight.get("csv") \
                if isinstance(preflight, dict) else None
            if isinstance(expected, dict) and any(
                    binding.get(key) != expected.get(key)
                    for key in (
                        "device", "inode", "uid", "gid", "mode", "nlink")):
                raise TransitionError(
                    "recovery archive no longer binds the preflight CSV inode")

        def preserve_stale_pid() -> Optional[Dict[str, object]]:
            canonical_exists = self.plan.old_pid_file.exists() or \
                self.plan.old_pid_file.is_symlink()
            archived_exists = self.plan.old_stale_pid_archive.exists() or \
                self.plan.old_stale_pid_archive.is_symlink()
            if canonical_exists and archived_exists:
                raise TransitionError("stale PID archive collision during recovery")
            if canonical_exists:
                pid_binding = file_binding(
                    self.plan.old_pid_file, with_hash=True)
                preflight = getattr(self, "old_preflight", None)
                expected_pid = preflight.get("pid_file") \
                    if isinstance(preflight, dict) else None
                if isinstance(expected_pid, dict) and \
                        pid_binding != expected_pid:
                    raise TransitionError(
                        "recovery stale PID no longer binds the preflight inode")
                rename_noreplace(
                    self.plan.old_pid_file, self.plan.old_stale_pid_archive,
                    pid_binding, parent_uid=self.plan.controller_uid)
                return {"binding": pid_binding,
                        "path": str(self.plan.old_stale_pid_archive)}
            if archived_exists:
                return {
                    "binding": file_binding(
                        self.plan.old_stale_pid_archive, with_hash=True),
                    "path": str(self.plan.old_stale_pid_archive),
                }
            return None

        if archive is not None:
            path = Path(str(archive["path"]))
            forced = archive.get("forced_stop")
            expected_path = self.plan.old_unclean_archive \
                if forced is True else self.plan.old_archive
            if type(forced) is not bool or path != expected_path:
                raise TransitionError("pre-dry archive path changed during recovery")
            other_path = self.plan.old_archive \
                if forced else self.plan.old_unclean_archive
            if other_path.exists() or other_path.is_symlink():
                raise TransitionError(
                    "alternate pre-dry archive appeared during recovery")
            binding = file_binding(path, with_hash=True)
            if binding != archive.get("binding"):
                raise TransitionError("pre-dry archive changed during recovery")
            validate_csv_binding(binding)
            result = dict(archive)
            if not forced and (result.get("stale_pid") is not None or
                               self.plan.old_pid_file.exists() or
                               self.plan.old_pid_file.is_symlink() or
                               self.plan.old_stale_pid_archive.exists() or
                               self.plan.old_stale_pid_archive.is_symlink()):
                raise TransitionError(
                    "graceful archive unexpectedly has stale PID evidence")
            result["stale_pid"] = result.get("stale_pid") or \
                (preserve_stale_pid() if forced else None)
            return result
        discovered = [path for path in (
            self.plan.old_archive, self.plan.old_unclean_archive)
            if path.exists() or path.is_symlink()]
        if len(discovered) > 1:
            raise TransitionError("multiple pre-dry archives exist during recovery")
        if discovered:
            path = discovered[0]
            binding = file_binding(path, with_hash=True)
            validate_csv_binding(binding)
            forced = path == self.plan.old_unclean_archive
            if not forced and (self.plan.old_pid_file.exists() or
                               self.plan.old_pid_file.is_symlink() or
                               self.plan.old_stale_pid_archive.exists() or
                               self.plan.old_stale_pid_archive.is_symlink()):
                raise TransitionError(
                    "graceful archive unexpectedly has stale PID evidence")
            record = {
                "binding": binding,
                "forced_stop": forced,
                "path": str(path),
                "stale_pid": preserve_stale_pid() if forced else None,
            }
            self.archive_record = record
            return record
        if not (self.plan.old_csv.exists() or self.plan.old_csv.is_symlink()):
            raise TransitionError("recovery lacks both old CSV and archive")
        binding = file_binding(self.plan.old_csv, with_hash=True)
        validate_csv_binding(binding)
        if self._fuser(self.plan.old_csv):
            raise TransitionError("old CSV gained a reader during recovery archive")
        destination = self.plan.old_unclean_archive
        rename_noreplace(
            self.plan.old_csv, destination, binding,
            parent_uid=self.plan.controller_uid)
        stale_pid = preserve_stale_pid()
        record = {"binding": binding, "forced_stop": True,
                  "path": str(destination), "stale_pid": stale_pid}
        self.archive_record = record
        return record

    def _replacement_launcher_roster(
            self, restored: Mapping[str, object]) -> list[Dict[str, object]]:
        session_id = restored.get("launcher_session")
        child_pid = restored.get("pid")
        launcher_tick = restored.get("launcher_start_tick")
        expected_command = self._replacement_launch_command()
        expected_command_sha256 = sha256_bytes(
            b"\0".join(value.encode("ascii")
                        for value in expected_command) + b"\0")
        if type(session_id) is not int or session_id <= 1 or \
                type(child_pid) is not int or child_pid <= 1 or \
                type(launcher_tick) is not int or launcher_tick < 0 or \
                restored.get("launcher_command") != list(expected_command) or \
                restored.get("launcher_command_sha256") != \
                    expected_command_sha256:
            raise TransitionError("replacement launcher identity is malformed")
        members = self._owned_session_members(session_id)
        allowed = {session_id, child_pid}
        if child_pid not in members or any(pid not in allowed for pid in members):
            raise TransitionError(
                "replacement launcher session roster has an unknown member")
        launcher_present = session_id in members
        receipt = self._inspect_replacement_identities(
            child_pid=child_pid,
            child_start_tick=int(restored.get("start_tick", 0)),
            launcher_pid=session_id, launcher_start_tick=launcher_tick,
            allowed_absence="none" if launcher_present else "launcher")
        child = receipt["targets"]["child"]
        launcher = receipt["targets"]["launcher"]
        if child.get("state") != "present" or \
                (launcher.get("state") == "present") != launcher_present:
            raise TransitionError(
                "replacement launcher roster identity changed")
        roster = [child["identity"]]
        if launcher_present:
            launcher_identity = restored.get("launcher_identity")
            identity = launcher["identity"]
            if launcher_tick == 0 or not isinstance(
                    launcher_identity, dict) or \
                    stable_process_identity(identity) != \
                    stable_process_identity(launcher_identity):
                raise TransitionError(
                    "replacement sudo launcher changed after creation")
            self._validate_replacement_launcher_identity(
                identity, session_id, launcher_tick)
            roster.insert(0, identity)
        return roster

    def _replacement_launch_command(self) -> Tuple[str, ...]:
        return replacement_launch_command(self.plan, self.tools)

    def _validate_replacement_launcher_identity(
            self, identity: Mapping[str, object], session_id: int,
            start_tick: int) -> None:
        command = self._replacement_launch_command()
        expected_hash = sha256_bytes(
            b"\0".join(value.encode("ascii") for value in command) + b"\0")
        exact = {
            "argv": list(command),
            "cmdline_sha256": expected_hash,
            "pid": session_id,
            "process_group": session_id,
            "session": session_id,
            "start_tick": start_tick,
            "uids": [self.plan.controller_uid, 0, 0, 0],
        }
        sudo_binding = self.tools["sudo"]["binding"]
        exact["executable"] = {
            "device": sudo_binding["device"],
            "inode": sudo_binding["inode"],
        }
        if any(identity.get(key) != value for key, value in exact.items()) or \
                identity.get("ppid") not in (os.getpid(), 1) or \
                not isinstance(identity.get("affinity"), str) or \
                not identity["affinity"]:
            raise TransitionError(
                "replacement sudo launcher identity is not exact")

    def _launch_old(self) -> Dict[str, object]:
        self._require_old_launcher_absent()
        if self._i2c_readers() or any(
                path.exists() or path.is_symlink()
                for path in (self.plan.old_csv, self.plan.old_pid_file)):
            raise TransitionError("old restart preconditions are not empty")
        code_seal = verify_code_seal(self.plan)
        source_binding = file_binding(self.plan.legacy_sampler, with_hash=True)
        expected_source = code_seal["files"]["legacy"]["binding"]
        if source_binding != expected_source or \
                source_binding["sha256"] != self.plan.old_source_sha256:
            raise TransitionError("root-sealed legacy source changed before restart")
        environment = execute_environment(self.plan)
        command = self._replacement_launch_command()
        launcher_start_tick = 0
        launcher_error: Optional[BaseException] = None
        boot_id = Path("/proc/sys/kernel/random/boot_id").read_text(
            encoding="ascii").strip()
        process = subprocess.Popen(
            command, env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True, close_fds=True)
        last_error: Optional[BaseException] = None
        restored: Optional[Dict[str, object]] = None
        try:
            try:
                launcher_start_tick = capture_owned_session_leader(process.pid)
            except BaseException as exc:
                launcher_start_tick = 0
                launcher_error = exc
            if launcher_error is not None:
                raise launcher_error
            bootstrap_deadline = self._recovery_wait_deadline(
                20.0, minimum_safety_wait_s=5.0)
            while self._recovery_now() < bootstrap_deadline:
                try:
                    raw_pid = stable_file_bytes(self.plan.old_pid_file)
                    if re.fullmatch(rb"[1-9][0-9]*\n", raw_pid) is None:
                        raise TransitionError("replacement old PID is noncanonical")
                    pid = int(raw_pid)
                    raw = stable_file_bytes(self.plan.old_csv, attempts=3)
                except (OSError, TransitionError, ValueError) as exc:
                    last_error = exc
                    time.sleep(0.05)
                    continue
                # The frozen classifier labels every no-newline prefix as
                # ``incomplete`` before considering its bytes.  NUL and
                # non-ASCII can never be a valid CSV prefix, so reject them
                # terminally instead of retrying malformed evidence until the
                # emergency deadline.
                if b"\0" in raw:
                    raise TransitionError(
                        "replacement old CSV bootstrap contains NUL")
                try:
                    raw.decode("ascii")
                except UnicodeDecodeError as exc:
                    raise TransitionError(
                        "replacement old CSV bootstrap is not ASCII") from exc
                csv_state = self.p32._bootstrap_thermal_csv_state(raw)
                if csv_state == "incomplete":
                    expected_header = (
                        ",".join(self.p32.THERMAL_FIELDS) + "\n").encode(
                            "ascii")
                    if not expected_header.startswith(raw):
                        raise TransitionError(
                            "replacement old CSV incomplete header prefix "
                            "is invalid")
                if csv_state in ("incomplete", "header"):
                    last_error = TransitionError(
                        "replacement old CSV has no complete sample: " +
                        csv_state)
                    time.sleep(0.05)
                    continue
                if csv_state != "row":
                    raise TransitionError(
                        "replacement old CSV bootstrap state is terminal: " +
                        str(csv_state))
                # The frozen P32 strict parser is deliberately outside the
                # readiness retry handler.  A complete malformed row, schema
                # mismatch, or implausible value is substantive and terminal.
                rows = self.p32._parse_thermal_csv(raw)
                try:
                    # This loop discovers only filesystem/reader readiness.
                    # It never reads or accepts the root child's /proc identity;
                    # the first such authority is the one-shot root helper below.
                    if self._i2c_readers() != (pid,) or \
                            self._fuser(self.plan.old_csv) != (pid,):
                        raise TransitionError(
                            "replacement old sampler readers are not ready")
                    csv_binding = file_binding(
                        self.plan.old_csv, with_hash=False)
                    pid_binding = file_binding(
                        self.plan.old_pid_file, with_hash=True)
                    if (csv_binding["uid"], csv_binding["mode"],
                            csv_binding["nlink"]) != (0, 0o444, 1) or \
                            (pid_binding["uid"], pid_binding["mode"],
                             pid_binding["nlink"], pid_binding["sha256"]) != (
                                0, 0o444, 1,
                                sha256_bytes(raw_pid)):
                        raise TransitionError(
                            "replacement old sampler evidence binding is not exact")
                    restored = {
                        "boot_id": boot_id,
                        "cmdline": list(self.plan.replacement_old_argv),
                        "cmdline_sha256":
                            self.plan.replacement_old_cmdline_sha256,
                        "csv_initial_size": csv_binding["size"],
                        "csv_live_identity": {
                            key: csv_binding[key] for key in (
                                "device", "gid", "inode", "mode", "nlink",
                                "uid")},
                        "csv_path": str(self.plan.old_csv),
                        "first_sample_monotonic_s": rows[0]["monotonic_s"],
                        "launcher_session": process.pid,
                        "launcher_start_tick": launcher_start_tick,
                        "pid": pid, "pid_binding": pid_binding,
                        "pid_path": str(self.plan.old_pid_file),
                        "source_sha256": self.plan.old_source_sha256,
                        "source_binding": source_binding,
                        "source_path": str(self.plan.legacy_sampler),
                        "start_tick": 0,
                    }
                    break
                except (OSError, TransitionError, ValueError) as exc:
                    last_error = exc
                    time.sleep(0.05)
            if restored is not None:
                pid = int(restored["pid"])
                launcher_returncode_before_identity = process.poll()
                if launcher_returncode_before_identity not in (None, 0):
                    raise TransitionError(
                        "replacement sudo launcher exited nonzero")
                replacement_receipt = self._inspect_replacement_identities(
                    child_pid=pid,
                    child_start_tick=0,
                    launcher_pid=process.pid,
                    launcher_start_tick=launcher_start_tick,
                    allowed_absence="launcher")
                child_target = replacement_receipt["targets"]["child"]
                launcher_target = replacement_receipt["targets"]["launcher"]
                if child_target.get("state") != "present":
                    raise TransitionError(
                        "replacement child disappeared at identity bootstrap")
                accepted_identity = child_target["identity"]
                if replacement_receipt["boot_id"] != boot_id:
                    raise TransitionError(
                        "replacement helper boot identity changed")
                self.replacement_identity_bootstrap = replacement_receipt
                restored.update({
                    "boot_id": replacement_receipt["boot_id"],
                    "cmdline": accepted_identity["argv"],
                    "cmdline_sha256": accepted_identity["cmdline_sha256"],
                    "identity_bootstrap": replacement_receipt,
                    "start_tick": accepted_identity["start_tick"],
                })
                self.legacy_fd_bootstrap = self._inspect_sampler_fds(
                    "legacy", pid=pid,
                    start_tick=accepted_identity["start_tick"],
                    boot_id=replacement_receipt["boot_id"],
                    argv=self.plan.replacement_old_argv,
                    csv_path=self.plan.old_csv)
                restored["fd_provenance_bootstrap"] = \
                    self.legacy_fd_bootstrap
                launcher_returncode = process.poll()
                if launcher_returncode not in (None, 0):
                    raise TransitionError(
                        "replacement sudo launcher exited nonzero")
                launcher_identity = launcher_target.get("identity") \
                    if launcher_target.get("state") == "present" else None
                if (launcher_returncode is None) != \
                        (launcher_identity is not None):
                    raise TransitionError(
                        "replacement launcher poll/identity state changed")
                if launcher_identity is not None:
                    if launcher_start_tick == 0:
                        raise TransitionError(
                            "replacement sudo launcher start identity was lost")
                    self._validate_replacement_launcher_identity(
                        launcher_identity, process.pid,
                        launcher_start_tick)
                restored["launcher_command"] = list(command)
                restored["launcher_command_sha256"] = sha256_bytes(
                    b"\0".join(value.encode("ascii")
                                for value in command) + b"\0")
                restored["launcher_identity"] = launcher_identity
                restored["launcher_returncode_before_identity"] = \
                    launcher_returncode_before_identity
                restored["launcher_returncode_at_accept"] = \
                    launcher_returncode
                restored["launcher_roster"] = \
                    self._replacement_launcher_roster(restored)
                return restored
        except BaseException:
            try:
                self.p32._kill_owned_process_session(
                    process, launcher_start_tick, boot_id,
                    self.plan.root, self.design)
            except BaseException as cleanup:
                raise TransitionError(
                    "old restart bootstrap and exact-session cleanup failed: %s" %
                    cleanup) from cleanup
            raise
        try:
            self.p32._kill_owned_process_session(
                process, launcher_start_tick, boot_id,
                self.plan.root, self.design)
        except BaseException as cleanup:
            raise TransitionError(
                "old restart timeout and exact-session cleanup failed: %s" %
                cleanup) from cleanup
        raise TransitionError(
            "replacement old sampler bootstrap failed: %s" % last_error)

    def restore_old(
            self, archive: Optional[Mapping[str, object]]) -> Mapping[str, object]:
        if not self.candidate_cleanup_complete:
            # Still independently enforce empty I2C below; this flag makes an
            # omitted cleanup attempt visible and fail closed.
            raise TransitionError("old recovery was entered before candidate cleanup")
        quiesce = self._quiesce_old_for_recovery()
        sealed_archive = self._ensure_recovery_archive(archive)
        restored = self._launch_old()
        restored["original_quiesce"] = quiesce
        restored["archive_record"] = sealed_archive
        restored["archived_pre_dry_sha256"] = \
            sealed_archive["binding"]["sha256"]
        restored["archived_pre_dry_path"] = sealed_archive["path"]
        return restored

    def _revalidate_restored_old(
            self, restored: Mapping[str, object]) -> Dict[str, object]:
        pid = restored.get("pid")
        if type(pid) is not int or pid <= 1:
            raise TransitionError("restored old PID is malformed")
        receipt = self._inspect_replacement_identities(
            child_pid=pid,
            child_start_tick=int(restored.get("start_tick", 0)),
            launcher_pid=int(restored.get("launcher_session", 0)),
            launcher_start_tick=int(restored.get("launcher_start_tick", 0)),
            allowed_absence="launcher")
        child = receipt["targets"]["child"]
        if child.get("state") != "present":
            raise TransitionError(
                "restored old sampler disappeared before audit binding")
        identity = child["identity"]
        if identity["argv"] != list(self.plan.replacement_old_argv) or \
                identity["cmdline_sha256"] != \
                    self.plan.replacement_old_cmdline_sha256 or \
                identity["cmdline_sha256"] != restored.get("cmdline_sha256") or \
                identity["start_tick"] != restored.get("start_tick") or \
                receipt["boot_id"] != restored.get("boot_id") or \
                identity["session"] != restored.get("launcher_session") or \
                identity["process_group"] != restored.get("launcher_session") or \
                self._i2c_readers() != (pid,) or \
                self._fuser(self.plan.old_csv) != (pid,):
            raise TransitionError("restored old sampler changed before audit binding")
        pid_binding = file_binding(self.plan.old_pid_file, with_hash=True)
        source_binding = file_binding(self.plan.legacy_sampler, with_hash=True)
        expected_source = verify_code_seal(
            self.plan)["files"]["legacy"]["binding"]
        if pid_binding != restored.get("pid_binding") or \
                not isinstance(expected_source, dict) or \
                source_binding != expected_source or \
                source_binding["sha256"] != self.plan.old_source_sha256 or \
                restored.get("source_binding") != source_binding or \
                restored.get("source_path") != str(self.plan.legacy_sampler):
            raise TransitionError("restored old PID/source binding changed")
        bootstrap = restored.get("fd_provenance_bootstrap")
        if not isinstance(bootstrap, dict) or \
                bootstrap != self.legacy_fd_bootstrap:
            raise TransitionError(
                "restored old lacks legacy bootstrap FD provenance")
        self._wait_for_legacy_csv_growth(bootstrap)
        csv_raw = stable_file_bytes(self.plan.old_csv, attempts=3)
        csv_rows = self.p32._parse_thermal_csv(csv_raw)
        csv_binding = file_binding(self.plan.old_csv, with_hash=False)
        csv_identity = {
            key: csv_binding[key] for key in (
                "device", "gid", "inode", "mode", "nlink", "uid")}
        if csv_identity != restored.get("csv_live_identity"):
            raise TransitionError("restored old CSV inode changed before audit binding")
        initial_size = restored.get("csv_initial_size")
        if type(initial_size) is not int or initial_size <= 0 or \
                csv_binding["size"] < initial_size:
            raise TransitionError("restored old CSV size regressed before audit binding")
        if not csv_rows or len(csv_raw) < initial_size or \
                csv_binding["size"] < len(csv_raw):
            raise TransitionError("restored old CSV content changed before audit binding")
        csv_times = [float(row["monotonic_s"]) for row in csv_rows]
        if csv_times != sorted(csv_times) or \
                any(right <= left or right - left > 2.25
                    for left, right in zip(csv_times, csv_times[1:])) or \
                not 0.0 <= time.monotonic() - csv_times[-1] <= 2.25:
            raise TransitionError("restored old CSV cadence is not live")
        self.legacy_fd_replay = self._inspect_sampler_fds(
            "legacy", pid=pid, start_tick=identity["start_tick"],
            boot_id=receipt["boot_id"], argv=self.plan.replacement_old_argv,
            csv_path=self.plan.old_csv)
        self._replay_fd_growth(
            bootstrap, self.legacy_fd_replay, description="legacy")
        return {"csv_current_size": csv_binding["size"],
                "csv_stable_sample_count": len(csv_rows),
                "csv_stable_sha256": sha256_bytes(csv_raw),
                "csv_stable_size": len(csv_raw),
                "csv_stable_last_monotonic_s":
                    csv_rows[-1]["monotonic_s"],
                "identity": identity,
                "identity_receipt": receipt,
                "fd_provenance": {
                    "bootstrap": bootstrap,
                    "replay": self.legacy_fd_replay,
                },
                "source": {"binding": source_binding,
                           "loaded_equality_basis":
                               "root-custody-normal-path-and-exact-argv; "
                               "no-retained-source-fd-claimed",
                           "path": str(self.plan.legacy_sampler)},
                "launcher_roster":
                    self._replacement_launcher_roster(restored)}

    def _replay_external_state(
            self, candidate_accept: Optional[Mapping[str, object]],
            archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object],
    ) -> Dict[str, object]:
        """Re-open every safety-critical artifact and live ownership roster."""
        frozen_plan = verify_transition_plan(self.plan)
        if frozen_plan.get("tools") != self.tools:
            raise TransitionError(
                "frozen transition plan changed during external replay")
        plan_binding = file_binding(self.plan.plan_receipt, with_hash=True)
        if plan_binding["uid"] != self.plan.controller_uid or \
                plan_binding["mode"] != 0o444 or plan_binding["nlink"] != 1:
            raise TransitionError(
                "transition plan receipt binding changed during replay")
        archive_record = archive if archive is not None else \
            restored.get("archive_record")
        if not isinstance(archive_record, dict) or \
                not isinstance(archive_record.get("path"), str) or \
                not isinstance(archive_record.get("binding"), dict):
            raise TransitionError("external replay lacks the pre-dry archive")
        archive_path = Path(str(archive_record["path"]))
        forced_stop = archive_record.get("forced_stop")
        expected_archive_path = self.plan.old_unclean_archive \
            if forced_stop is True else self.plan.old_archive
        if type(forced_stop) is not bool or \
                archive_path != expected_archive_path:
            raise TransitionError("external replay archive destination changed")
        alternate_archive = self.plan.old_archive \
            if forced_stop else self.plan.old_unclean_archive
        if alternate_archive.exists() or alternate_archive.is_symlink():
            raise TransitionError(
                "alternate pre-dry archive appeared during external replay")
        archive_binding = file_binding(archive_path, with_hash=True)
        if archive_binding != archive_record["binding"] or \
                archive_binding["uid"] != 0 or \
                archive_binding["mode"] != 0o444 or \
                archive_binding["nlink"] != 1 or \
                restored.get("archived_pre_dry_path") != str(archive_path) or \
                restored.get("archived_pre_dry_sha256") != \
                    archive_binding["sha256"]:
            raise TransitionError("pre-dry archive changed during external replay")
        stale_result = None
        stale_pid = archive_record.get("stale_pid")
        if stale_pid is not None:
            if not forced_stop:
                raise TransitionError(
                    "graceful archive gained stale PID evidence")
            if not isinstance(stale_pid, dict) or \
                    stale_pid.get("path") != \
                        str(self.plan.old_stale_pid_archive) or \
                    not isinstance(stale_pid.get("binding"), dict):
                raise TransitionError("external replay stale PID receipt changed")
            stale_binding = file_binding(
                self.plan.old_stale_pid_archive, with_hash=True)
            if stale_binding != stale_pid["binding"] or \
                    stale_binding["uid"] != 0 or \
                    stale_binding["mode"] != 0o444 or \
                    stale_binding["nlink"] != 1:
                raise TransitionError("stale PID archive changed during replay")
            stale_result = {"binding": stale_binding,
                            "path": str(self.plan.old_stale_pid_archive)}
        elif self.plan.old_stale_pid_archive.exists() or \
                self.plan.old_stale_pid_archive.is_symlink():
            raise TransitionError(
                "unbound stale PID archive appeared during replay")
        source = file_binding(self.plan.old_source, with_hash=True)
        expected_source = self.old_preflight.get("source") \
            if isinstance(self.old_preflight, dict) else None
        if not isinstance(expected_source, dict) or source != expected_source or \
                source["sha256"] != self.plan.old_source_sha256:
            raise TransitionError("old source changed during external replay")
        candidate_result = self._replay_candidate_accept(
            candidate_accept, restored.get("pid"))
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed during final replay")
        return {
            "archive": {"binding": archive_binding,
                        "path": str(archive_path)},
            "candidate": candidate_result,
            "controller_interpreter": interpreter,
            "execute_runtime": runtime,
            "old_source": {"binding": source,
                           "path": str(self.plan.old_source)},
            "restored_old": self._revalidate_restored_old(restored),
            "stale_pid": stale_result,
            "thermal_parent": self._validate_thermal_parent(
                self.old_preflight.get("thermal_parent")
                if isinstance(self.old_preflight, dict) else None),
            "tools": tools,
            "transition_plan": {
                "binding": plan_binding,
                "path": str(self.plan.plan_receipt),
                "self_sha256_excluding_field":
                    frozen_plan["self_sha256_excluding_field"],
            },
        }

    def _replay_candidate_accept(
            self, candidate_accept: Optional[Mapping[str, object]],
            restored_pid: object,
    ) -> Optional[Dict[str, object]]:
        if self.candidate_owner is None:
            if candidate_accept is not None:
                raise TransitionError(
                    "candidate evidence exists without a launched owner")
            return None
        if not self.candidate_cleanup_complete or \
                type(restored_pid) is not int or restored_pid <= 1:
            raise TransitionError(
                "launched candidate lacks a completed ownership cleanup")
        current: Optional[Dict[str, object]] = None
        if candidate_accept is not None:
            current = dict(self.accept_candidate())
            if current != dict(candidate_accept):
                raise TransitionError(
                    "candidate final evidence changed after acceptance")
        observation = self._candidate_cleanup_observation(self.candidate_owner)
        # This replay runs after the old sampler has been restored, so global
        # I2C ownership must now be exactly that restored PID.  Candidate
        # identity/session/CSV ownership must remain absent independently.
        if observation.get("launcher_returncode") is None or \
                observation.get("candidate_pid_present") is not False or \
                observation.get("exact_identity_live") is not False or \
                observation.get("session_members") != [] or \
                observation.get("csv_readers") != [] or \
                observation.get("i2c_readers") != [restored_pid]:
            raise TransitionError(
                "candidate ownership reappeared during final replay")
        if candidate_accept is None:
            return {"acceptance": None, "ownership": observation}
        assert current is not None
        current["ownership"] = observation
        return current

    def publish_audit_binding(
            self, archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object], dry_accepted: bool,
            candidate_accept: Optional[Mapping[str, object]],
            receipt_chain_prefix: Mapping[str, object],
    ) -> Mapping[str, object]:
        if type(dry_accepted) is not bool or \
                (dry_accepted and not isinstance(candidate_accept, dict)):
            raise TransitionError(
                "audit acceptance claim lacks retained candidate evidence")
        if not isinstance(receipt_chain_prefix, dict) or \
                set(receipt_chain_prefix) != {"count", "head_sha256", "roster"} or \
                type(receipt_chain_prefix.get("count")) is not int or \
                receipt_chain_prefix["count"] <= 0 or \
                not isinstance(receipt_chain_prefix.get("head_sha256"), str) or \
                SHA256_RE.fullmatch(receipt_chain_prefix["head_sha256"]) is None or \
                not isinstance(receipt_chain_prefix.get("roster"), list) or \
                len(receipt_chain_prefix["roster"]) != \
                    receipt_chain_prefix["count"]:
            raise TransitionError("audit receipt-chain prefix is malformed")
        prepublication_replay = self._replay_external_state(
            candidate_accept, archive, restored)
        archive_record = archive if archive is not None else \
            restored.get("archive_record")
        if not isinstance(archive_record, dict) or \
                not isinstance(archive_record.get("path"), str) or \
                not isinstance(archive_record.get("binding"), dict):
            raise TransitionError("audit binding lacks the pre-dry archive")
        archive_path = str(archive_record["path"])
        if archive_path not in {
                str(self.plan.old_archive), str(self.plan.old_unclean_archive)}:
            raise TransitionError("audit archive destination changed")
        archive_binding = file_binding(Path(archive_path), with_hash=True)
        if archive_binding != archive_record["binding"] or \
                archive_binding["uid"] != 0 or \
                archive_binding["mode"] != 0o444 or \
                archive_binding["nlink"] != 1 or \
                restored.get("archived_pre_dry_path") != archive_path or \
                restored.get("archived_pre_dry_sha256") != \
                    archive_binding["sha256"]:
            raise TransitionError("pre-dry archive changed before audit binding")
        stale_pid = archive_record.get("stale_pid")
        if stale_pid is not None:
            if not isinstance(stale_pid, dict) or \
                    stale_pid.get("path") != \
                        str(self.plan.old_stale_pid_archive) or \
                    not isinstance(stale_pid.get("binding"), dict):
                raise TransitionError("stale PID archive receipt changed")
            stale_binding = file_binding(
                self.plan.old_stale_pid_archive, with_hash=True)
            if stale_binding != stale_pid["binding"] or \
                    stale_binding["uid"] != 0 or \
                    stale_binding["mode"] != 0o444 or \
                    stale_binding["nlink"] != 1:
                raise TransitionError("stale PID archive changed before audit binding")
        # Keep the live identity check after potentially lengthy archive hashing,
        # immediately before publishing the future-auditor handoff.
        tools_after_restore = self._verify_tools()
        interpreter_after_restore = verify_running_interpreter(
            tools_after_restore)
        live_revalidation = self._revalidate_restored_old(restored)
        value = sealed_record(AUDIT_BINDING_SCHEMA, {
            "archived_pre_dry": {
                "binding": archive_binding, "path": archive_path},
            "candidate_dry_accepted": dry_accepted,
            "candidate_accept": dict(candidate_accept)
                if candidate_accept is not None else None,
            "created_utc": utc_now(),
            "controller_interpreter": interpreter_after_restore,
            "live_revalidation": live_revalidation,
            "live_old_sampler": dict(restored),
            "prepublication_replay": prepublication_replay,
            "receipt_chain_prefix": dict(receipt_chain_prefix),
            "tools": tools_after_restore,
            "transition_id": self.plan.transition_id,
        })
        binding = write_new(
            self.plan.audit_binding, canonical_json(value),
            owner_uid=self.plan.controller_uid)
        replay = load_canonical(
            self.plan.audit_binding, AUDIT_BINDING_SCHEMA,
            "future audit binding")
        if replay != value:
            raise TransitionError("future audit binding did not replay")
        tools_after_audit = verify_tool_records(replay.get("tools"))
        interpreter_after_audit = verify_running_interpreter(
            tools_after_audit)
        if tools_after_audit != tools_after_restore or \
                tools_after_audit != self.tools or \
                interpreter_after_audit != interpreter_after_restore:
            raise TransitionError(
                "tool identities changed after audit publication")
        return {"binding": binding, "path": str(self.plan.audit_binding),
                "controller_interpreter_after_audit":
                    interpreter_after_audit,
                "tools_after_audit": tools_after_audit, "value": value}

    def final_replay(
            self, candidate_accept: Optional[Mapping[str, object]],
            archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object], audit: Mapping[str, object],
    ) -> Mapping[str, object]:
        """Replay all artifacts and ownership after audit publication."""
        if audit.get("path") != str(self.plan.audit_binding) or \
                not isinstance(audit.get("binding"), dict) or \
                not isinstance(audit.get("value"), dict):
            raise TransitionError("postpublication audit receipt is malformed")
        replay = self._replay_external_state(
            candidate_accept, archive, restored)
        audit_binding = file_binding(self.plan.audit_binding, with_hash=True)
        audit_value = load_canonical(
            self.plan.audit_binding, AUDIT_BINDING_SCHEMA,
            "postpublication future audit binding")
        if audit_binding != audit["binding"] or audit_value != audit["value"]:
            raise TransitionError("future audit binding changed after publication")
        if audit_value.get("prepublication_replay") != \
                audit["value"].get("prepublication_replay") or \
                audit_value.get("candidate_accept") != \
                (dict(candidate_accept)
                 if candidate_accept is not None else None):
            raise TransitionError("audit replay contract changed")
        return {**replay,
                "audit": {"binding": audit_binding,
                          "path": str(self.plan.audit_binding),
                          "value_sha256": audit_binding["sha256"]}}


def execute_transition(plan: TransitionPlan) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)
    # Retained below as a source-level historical model.  This unreachable
    # body cannot become an execution authority for the rescued live sampler.
    execute_runtime = verify_execute_runtime(plan)
    deadline = Deadline(plan.deadline_s, plan.recovery_reserve_s)
    value = verify_transition_plan(plan)
    if Path(__file__).resolve() != plan.controller.resolve():
        raise TransitionError("execute mode requires the frozen controller path")
    tools = value.get("tools")
    if not isinstance(tools, dict):
        raise TransitionError("verified transition plan lacks its tool ledger")
    verify_running_interpreter(tools, require_exact_path=True)
    p32 = load_verified_p32(plan.p32, plan.p32_sha256)
    journal = ReceiptJournal(
        plan.receipts, plan.transition_id, plan.controller_uid, deadline)
    backend = LiveBackend(plan, p32, deadline, tools, execute_runtime)
    lock_path = plan.root / "transition.lock"
    lock_fd = os.open(
        str(lock_path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
        getattr(os, "O_NOFOLLOW", 0),
        0o400)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.fsync(lock_fd)
        fsync_directory(plan.root)
        with DeferredSignals() as guard:
            controller = TransitionController(
                plan, backend, journal, deadline, guard)
            return controller.run()
    finally:
        os.close(lock_fd)


def _decode_target_argv(raw_hex: str) -> Tuple[str, ...]:
    try:
        raw = bytes.fromhex(raw_hex)
        value = json.loads(raw.decode("ascii"))
    except (ValueError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TransitionError("FD inspection target argv encoding is malformed") \
            from exc
    if not isinstance(value, list) or not value or \
            not all(isinstance(item, str) and item and "\0" not in item
                    for item in value) or canonical_json(value) != raw:
        raise TransitionError("FD inspection target argv is noncanonical")
    return tuple(value)


def inspection_environment(plan: TransitionPlan) -> Dict[str, str]:
    return {**execute_environment(plan), "PYTHONDONTWRITEBYTECODE": "1"}


def execute_fd_inspection_mode(
        plan: TransitionPlan, *, kind: str, pid: int, start_tick: int,
        boot_id: str, csv_path: Path, target_argv_hex: str,
) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)
    # Historical helper model only; lower pure capture/verifier functions are
    # exercised directly by source-only tests.
    if Path(__file__).resolve() != plan.controller.resolve() or \
            sha256_file(Path(__file__).resolve()) != plan.controller_sha256:
        raise TransitionError(
            "FD inspection mode requires the reviewed root-sealed controller")
    if os.environ != inspection_environment(plan):
        raise TransitionError("FD inspection helper environment changed")
    flags = {name: getattr(sys.flags, name) for name in EXECUTE_FLAG_CONTRACT}
    if flags != EXECUTE_FLAG_CONTRACT:
        raise TransitionError("FD inspection helper Python flags changed")
    target_argv = _decode_target_argv(target_argv_hex)
    expected_orig_argv = (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--inspect-sealed-sampler-fds", kind,
        "--target-pid", str(pid),
        "--target-start-tick", str(start_tick),
        "--target-boot-id", boot_id,
        "--target-csv", str(csv_path),
        "--target-argv-json-hex", target_argv_hex,
        "--expected-controller-sha256", plan.controller_sha256,
    )
    if tuple(sys.orig_argv) != expected_orig_argv:
        raise TransitionError("FD inspection helper argv changed")
    tools = capture_tool_records()
    verify_running_interpreter(tools, require_exact_path=True)
    return inspect_sampler_fd_provenance(
        plan, kind=kind, pid=pid, start_tick=start_tick,
        boot_id=boot_id, argv=target_argv, csv_path=csv_path)


def execute_identity_inspection_mode(
        plan: TransitionPlan, request: Mapping[str, object],
) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)
    # Historical root-only, isolated identity-helper model; unreachable in V6.
    if os.geteuid() != 0 or Path(__file__).resolve() != plan.controller.resolve() \
            or sha256_file(Path(__file__).resolve()) != plan.controller_sha256:
        raise TransitionError(
            "identity inspection requires the root-sealed root helper")
    canonical = identity_inspection_request(
        plan, profile=str(request.get("profile")),
        child_pid=request.get("child_pid"),
        child_start_tick=request.get("child_start_tick"),
        launcher_pid=request.get("launcher_pid"),
        launcher_start_tick=request.get("launcher_start_tick"),
        controller_pid=request.get("controller_pid"),
        allowed_absence=str(request.get("allowed_absence")))
    tools = capture_tool_records()
    static_runtime = identity_inspection_runtime_contract(
        plan, canonical, tools)
    if os.environ != static_runtime["environment"]:
        raise TransitionError("identity inspection helper environment changed")
    flags = {name: getattr(sys.flags, name) for name in EXECUTE_FLAG_CONTRACT}
    if flags != EXECUTE_FLAG_CONTRACT:
        raise TransitionError("identity inspection helper Python flags changed")
    if list(sys.orig_argv) != static_runtime["sys_orig_argv"]:
        raise TransitionError("identity inspection helper argv changed")
    interpreter = verify_running_interpreter(tools, require_exact_path=True)
    if interpreter != static_runtime["controller_interpreter"]:
        raise TransitionError("identity inspection interpreter changed")
    # This replay is deliberately before the first target-proc read.  The
    # complete receipt is carried in both helper runtime and result evidence.
    code_seal, captured = capture_process_identity_targets_under_seal(
        plan, canonical, tools)
    runtime = _runtime_with_code_seal(static_runtime, code_seal)
    result = sealed_record(IDENTITY_INSPECTION_SCHEMA, {
        **captured,
        "code_seal": code_seal,
        "helper_runtime": runtime,
        "pidfd_policy": IDENTITY_PIDFD_POLICY,
        "request": canonical,
        "tools": tools,
        "transition_id": plan.transition_id,
    })
    return verify_identity_inspection_receipt(
        plan, result, canonical, tools)


def execute_process_signal_mode(
        plan: TransitionPlan, request: Mapping[str, object],
) -> Dict[str, object]:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)
    # Historical held-pidfd signal/receipt model; unreachable in V6.
    if os.geteuid() != 0 or Path(__file__).resolve() != plan.controller.resolve() \
            or sha256_file(Path(__file__).resolve()) != plan.controller_sha256:
        raise TransitionError("process signal requires the root-sealed helper")
    canonical = process_signal_request(
        plan, request["identity_request"], target=str(request.get("target")),
        signum=request.get("signal"))
    tools = capture_tool_records()
    static_runtime = process_signal_runtime_contract(
        plan, canonical, tools)
    if os.environ != static_runtime["environment"]:
        raise TransitionError("process signal helper environment changed")
    flags = {name: getattr(sys.flags, name) for name in EXECUTE_FLAG_CONTRACT}
    if flags != EXECUTE_FLAG_CONTRACT:
        raise TransitionError("process signal helper Python flags changed")
    if list(sys.orig_argv) != static_runtime["sys_orig_argv"]:
        raise TransitionError("process signal helper argv changed")
    interpreter = verify_running_interpreter(tools, require_exact_path=True)
    if interpreter != static_runtime["controller_interpreter"]:
        raise TransitionError("process signal helper interpreter changed")
    code_seal, captured = capture_process_identity_targets_under_seal(
        plan, canonical["identity_request"], tools,
        signal_target=canonical["target"], signum=canonical["signal"])
    runtime = _runtime_with_code_seal(static_runtime, code_seal)
    result = sealed_record(PROCESS_SIGNAL_SCHEMA, {
        **captured,
        "code_seal": code_seal,
        "helper_runtime": runtime,
        "pidfd_policy": PROCESS_SIGNAL_PIDFD_POLICY,
        "request": canonical,
        "tools": tools,
        "transition_id": plan.transition_id,
    })
    # Do not perform a new pathname/seal/proc operation after the successful
    # signal.  Validate entirely against the full seal receipt captured before
    # the target was touched; the unprivileged caller independently replays
    # the current seal when it consumes stdout.
    return _verify_process_signal_receipt_with_code_seal(
        plan, result, canonical, tools, code_seal)


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--print-root-code-seal-stage-plan",
                       action="store_true")
    modes.add_argument("--prepare-sealed-transition", action="store_true")
    modes.add_argument("--execute-sealed-transition")
    modes.add_argument("--inspect-sealed-sampler-fds",
                       choices=("candidate", "legacy"))
    modes.add_argument("--inspect-sealed-process-identities",
                       choices=("original", "replacement"))
    modes.add_argument("--signal-sealed-process-identities",
                       choices=("original",))
    parser.add_argument("--expected-controller-sha256", required=True)
    parser.add_argument("--confirmation")
    parser.add_argument("--target-pid", type=int)
    parser.add_argument("--target-start-tick", type=int)
    parser.add_argument("--target-boot-id")
    parser.add_argument("--target-csv")
    parser.add_argument("--target-argv-json-hex")
    parser.add_argument("--identity-child-pid", type=int)
    parser.add_argument("--identity-child-start-tick", type=int)
    parser.add_argument("--identity-launcher-pid", type=int)
    parser.add_argument("--identity-launcher-start-tick", type=int)
    parser.add_argument("--identity-controller-pid", type=int)
    parser.add_argument("--identity-allowed-absence",
                        choices=tuple(IDENTITY_ABSENCE_POLICIES))
    parser.add_argument("--signal-target", choices=("child", "launcher"))
    parser.add_argument("--signal-number", type=int,
                        choices=(int(signal.SIGTERM), int(signal.SIGKILL)))
    args = parser.parse_args(argv)
    if SHA256_RE.fullmatch(args.expected_controller_sha256) is None:
        parser.error("expected controller SHA256 must be canonical lowercase hex")
    parser.error(SOURCE_ONLY_RETIREMENT)


def main(argv: Optional[Sequence[str]] = None) -> int:
    raise TransitionError(SOURCE_ONLY_RETIREMENT)


if __name__ == "__main__":
    raise SystemExit(main())
