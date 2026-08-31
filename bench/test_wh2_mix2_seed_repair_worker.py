#!/usr/bin/env python3
"""Bounded v9 native-worker/controller compatibility test at K=2."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bench import wh2_mix2_seed_repair_contract as subject


K2_BASE_PACKET_SEED = "0x9799409a"
K2_BASE_PRECODE_SEED = "0xecbb2b8a0e18da1e"
K2_EFFECTIVE_PACKET_SEED = "0x278c881b"
K2_EFFECTIVE_PRECODE_SEED = "0x7cae730f87b736db"
K2_SOURCE_SHA256 = \
    "91c0640fc5acf5a88835f4c71c88e55ad6a9b164069a14453dbbbc45c33ee043"
K2_SELECTION_TRACE_GOLDENS = (
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (13, "e5ce1cd184f9ff951dcfcf196f1656baa15fadb745b4dd94f05fb5d1e10375ae"),
    (13, "ff58d29f143b0d1d52f6a6f25b0c2ef53ed4a123898639918cf72091ea6d3b5a"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (12, "84bb304d569928ea6c624690bdfe1f972f49861821b606b4b7706d72c30cf2be"),
    (12, "f3b97ab9a14079525f550e1d6b68f16c2b9a31023abf691f0d850bd56f5c98d1"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (12, "21a4523d90f94baaaf7b5e195575d7ed66a5c72d83c9a5e19e2db6eb305b92f5"),
    (12, "3a113628e89ab786ce8d16165834d4014d8299671733ce4390db6fb009fd8169"),
    (14, "95ce951131bbd0ddc276c486e2f0ade6991328d794f6b4f17f026988a3b5f4c2"),
    (11, "3f11d8e542fb2873d6ba63c2cbf754b25f7e4d812209b0ef0773d1dbcbab0bb5"),
    (11, "0097a00119bd08309f9005b37373dd7646a6086d280e383f611f8b66b8774396"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (9, "0205e0970f104fede45f7cc062f3ec104a3e12498a05047384e90c4c6d514f4f"),
    (9, "b1330ef09dff6b3b7bd9af6a8bb1b98b845c608980264e9079db11d7d5bd8c9d"),
    (22, "65a956290dc94a53f7009673eea4350bdb1d6956ba76fffaef1b0c6bb3c3402b"),
    (14, "c53c720feef0da7e34ba75d26b636bc1a467fd931aed7349d1f27f4c7fef7757"),
    (14, "32dc33b76b4d035ab13dc7d90da37fa2a7c14fc7b31dd64c30c3f263ee222a48"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (9, "887f0c1f24c4e93f15584b5f1b6d5cca650800f0fec53a7bf3c625b7ed5d7b88"),
    (9, "515d0225d4353d31d902606ef33b928994a28538c430ce7a4ff4f187d6054ae8"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (9, "88ac07201b4c3d1197d58477321e0bf39861a995a0a2c223b78de2c891fe3ff7"),
    (9, "c6ee91ea3bffdd3544eb6c736a4a14253d3bd1fd2aba1c8feec8ccbeedac53fa"),
    (6, "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
    (10, "6d9a15c58549a3d55bdeda0741cbaae5cb743756d05ce4f30a6de9e34e0327c6"),
    (10, "cffcf16f52370f926662839ae1f6e48f9f8db2202f9957cab4d9de131f270d1f"),
)
K2_LOWER_FAILURE_GOLDENS = (
    (0, 0, "construct_failed", None),
    (0, 0, "construct_failed", None),
    (1, 0, "need_more_at_oh0", 1),
    (1, 0, "need_more_at_oh0", 1),
    (0, 0, "construct_failed", None),
    (0, 0, "construct_failed", None),
    (0, 0, "construct_failed", None),
    (0, 0, "construct_failed", None),
    (16, 5, "need_more_at_oh0", 1),
)


def parse_line(raw, context):
    if not raw.endswith(b"\n"):
        raise AssertionError("{} is not newline terminated".format(context))
    value = json.loads(raw)
    if raw != (subject.canonical_json(value) + "\n").encode("ascii"):
        raise AssertionError("{} is not canonical JSON".format(context))
    return value


def check_K2_construction_identity(value, attempt, context):
    if value["base_packet_seed"] != K2_BASE_PACKET_SEED or \
            value["base_precode_seed"] != K2_BASE_PRECODE_SEED:
        raise AssertionError(
            "{} did not derive its base from SelectSeedProfile(2,2)".format(
                context))
    base_packet = int(K2_BASE_PACKET_SEED[2:], 16)
    base_precode = int(K2_BASE_PRECODE_SEED[2:], 16)
    if value["effective_packet_seed"] != \
            subject._effective_packet_seed(base_packet, attempt) or \
            value["effective_precode_seed"] != \
            subject._effective_precode_seed(base_precode, attempt):
        raise AssertionError(
            "{} reported the wrong attempt step".format(context))
    if value["source_sha256"] != K2_SOURCE_SHA256:
        raise AssertionError(
            "{} used the wrong deterministic source".format(context))


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_wh2_mix2_seed_repair_worker.py WORKER")
    path = Path(sys.argv[1])
    worker_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    described = subprocess.run(
        [str(path), "--describe"], check=True, stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if described.stderr:
        raise AssertionError("worker description emitted stderr")
    description = parse_line(described.stdout, "worker description")
    contract = subject.load_contract()
    if (description.get("validation_roster_schema") !=
            subject.VALIDATION_ROSTER_SCHEMA or
            description.get("validation_roster_sha256") !=
            subject.VALIDATION_ROSTER_SHA256 or
            subject.validation_roster_sha256() !=
            subject.VALIDATION_ROSTER_SHA256):
        raise AssertionError(
            "worker description did not preflight-bind the final V roster")
    permuted_roots = list(subject.VALIDATION_ROOTS)
    permuted_roots[0], permuted_roots[1] = (
        permuted_roots[1], permuted_roots[0])
    if subject.validation_roster_sha256(permuted_roots) == \
            description["validation_roster_sha256"]:
        raise AssertionError("validation-root permutation kept the roster hash")
    worker = subject.verify_worker(
        path, worker_sha256, description["source_git_commit"], contract)
    process = None
    try:
        process = subprocess.Popen(
            ["/proc/self/fd/{}".format(worker.descriptor), "--worker"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, pass_fds=(worker.descriptor,),
            start_new_session=True)
        if (process.stdin is None or process.stdout is None or
                process.stderr is None):
            raise AssertionError("persistent worker pipes are unavailable")
        process.stdin.write(b"D 0 2\n")
        process.stdin.flush()
        derivation = parse_line(
            process.stdout.readline(subject.MAX_RECORD_BYTES + 1),
            "K=2 derivation")
        selected = subject.validate_derivation_record(
            derivation, contract, worker.sha256, expected_K=2)
        if selected != 9:
            raise AssertionError("K=2 did not select the first joint pass, 9")
        check_K2_construction_identity(derivation, selected, "K=2 derivation")
        if (derivation["effective_packet_seed"] !=
                K2_EFFECTIVE_PACKET_SEED or
                derivation["effective_precode_seed"] !=
                K2_EFFECTIVE_PRECODE_SEED):
            raise AssertionError("K=2 attempt-9 effective seeds changed")
        observed_selection_goldens = tuple(
            (cell["attempted_candidates"], cell["trace_sha256"])
            for cell in derivation["selected_successes"])
        if observed_selection_goldens != K2_SELECTION_TRACE_GOLDENS:
            raise AssertionError("K=2 27-cell selection trace identity changed")
        observed_lower_goldens = tuple(
            (cell["cell_ordinal"], cell["root_index"], cell["outcome"],
             cell["decoded_extra"])
            for cell in derivation["lower_attempt_failure_witnesses"])
        if observed_lower_goldens != K2_LOWER_FAILURE_GOLDENS:
            raise AssertionError("K=2 lower-attempt witnesses changed")
        process.stdin.write(b"Q\n")
        process.stdin.flush()
        process.stdin.close()
        if process.wait(timeout=30.0) != 0 or process.stderr.read():
            raise AssertionError("persistent worker did not exit cleanly")
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait(timeout=5.0)
        os.close(worker.descriptor)

    # Keep the reserved validation roots unobserved in this ordinary test.
    # This only reaches command parsing: the missing attempt is rejected before
    # ValidationRecord can construct a trace or probe any holdout cell.
    malformed = subprocess.run(
        [str(path), "--worker"], input=b"V 0 2\n",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if (malformed.returncode == 0 or malformed.stdout or
            b"malformed command" not in malformed.stderr):
        raise AssertionError("malformed V command was not rejected pre-probe")

    # The 65th byte is rejected by the bounded reader, before command parsing
    # can dispatch any construction or validation work.
    oversized = subprocess.run(
        [str(path), "--worker"], input=b"D" + b"0" * 64 + b"\n",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if (oversized.returncode == 0 or oversized.stdout or
            b"command exceeds 64-byte limit" not in oversized.stderr):
        raise AssertionError("oversized command was not rejected pre-probe")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
