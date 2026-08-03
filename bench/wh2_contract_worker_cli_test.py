#!/usr/bin/env python3
"""Focused CLI goldens for native WH2 recovery and qualified timing work."""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import unittest


DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v1"
RAW_DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v3"
DESCRIPTOR_SCHEMA = "wirehair.wh2.native-arm-descriptor.v1"
RECOVERY_SCHEMA = "wirehair.wh2.native-recovery-record.v1"
RAW_RECOVERY_SCHEMA = "wirehair.wh2.native-recovery-record.v3"
TIMING_SCHEMA = "wirehair.wh2.native-timing-record.v4"
TIMING_QUALIFICATION_SCHEMA = (
    "wirehair.wh2.native-timing-qualification-record.v1"
)
NATIVE_WORK_SCHEMA = "wirehair.wh2.native-work.v1"
TIMING_QUALIFICATION_MESSAGE_DOMAIN = (
    b"wirehair.wh2.timing-qualification-message.v1\0"
)
TIMING_PROXY_SCHEMA = "wirehair.wh2.native-timing-proxy-witness.v2"
TIMING_PROXY_DOMAIN = b"wirehair.wh2.timing-proxy-domain.v2\0"
SOURCE_DOMAIN = b"wirehair.wh2.source.v1\0"
RAW_REALIZED_SCHEMA = "wirehair.wh2.raw-realized-construction.v2"
REALIZED_SCHEMA = "wirehair.wh2.realized-construction.v1"
ZERO_SHA256 = "0" * 64
SHA256 = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT = re.compile(r"^[0-9a-f]{40}$")
UINT64_MASK = (1 << 64) - 1
LOSS_RETRY_STRIDE = 0x9E3779B97F4A7C15
RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_CANONICAL = (
    '{"excluded_inputs":["K","block_bytes","production_peel_seed",'
    '"production_dense_seed","production_seed_fixups"],'
    '"inputs":["precode_attempt","packet_attempt"],'
    '"packet_attempt_stride":"0x9e3779b9",'
    '"packet_base_seed":"0x4ec72102",'
    '"precode_attempt_stride":"0x9e3779b97f4a7c15",'
    '"precode_base_seed":"0x487468302aad7105",'
    '"schema":"wirehair.wh2.uniform-raw-seed-schedule.v1"}'
)
RAW_SEED_SCHEDULE_SHA256 = (
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7"
)

# Exact K=8, 64-byte structure-only witness golden.  These values bind the
# serialization and native seed resolution, not merely SHA syntax.  Any
# intentional serialization change must bump the witness schema and domain.
TIMING_PROXY_K8_B64_GOLDEN = {
    "K": 8,
    "block_bytes": 64,
    "construction_attempt": 0,
    "normalized_structure_sha256":
        "cd7311fd7b8bea997428148e6331b8cc2d76cc52585d276f69fc65ad8327e20c",
    "production_timing_configuration_sha256":
        "899dec1ce827dd5f85deff01ab3c88af8f180948bb95b5b79134c160e8eedea4",
    "production_timing_equation_system_sha256":
        "994637c26be404e790c6232bcb6c97dbc4576c0afce9a96f92e54b4f99a86434",
    "production_timing_packet_seed": "0xf57e654d",
    "production_timing_precode_seed": "0x94c06c665dbaad3a",
    "raw_recovery_packet_seed": "0x4ec72102",
    "raw_recovery_precode_seed": "0x487468302aad7105",
    "seeds_differ": True,
}
TIMING_PROXY_CELLS_SHA256 = (
    "6701e65f9a81c3906d6c4f9367e15881526bd85ce3e46b79439edc23497611ca"
)

CONTROLS = [
    {
        "arm": "wirehair2_head",
        "arm_descriptor_sha256":
            "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a",
        "codec": "wirehair2_certified",
    },
    {
        "arm": "wirehair1",
        "arm_descriptor_sha256":
            "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d",
        "codec": "wirehair1",
    },
]

CANDIDATES = {
    "d12-h11-periodic": {
        "arm": "wirehair2_raw_d12_h11_periodic",
        "arm_descriptor_sha256":
            "91d7c1a558e1cf93b002fcf2062b7657d301faca03972215495bdf2429499e90",
        "codec": "wirehair2_experiment",
        "dense_rows": 12,
        "heavy_rows": 11,
        "heavy_family": "periodic-cauchy",
    },
    "d12-h13-periodic": {
        "arm": "wirehair2_raw_d12_h13_periodic",
        "arm_descriptor_sha256":
            "7c7889747a97ac160726b807fb03349344d49d4bec84c9e8220aa4689b00d2ca",
        "codec": "wirehair2_experiment",
        "dense_rows": 12,
        "heavy_rows": 13,
        "heavy_family": "periodic-cauchy",
    },
    "d13-h12-periodic": {
        "arm": "wirehair2_raw_d13_h12_periodic",
        "arm_descriptor_sha256":
            "c70e0f57bb8d7783fa29b0decbed5da5058a8eb532d57d540f72108e114f091a",
        "codec": "wirehair2_experiment",
        "dense_rows": 13,
        "heavy_rows": 12,
        "heavy_family": "periodic-cauchy",
    },
}

RAW_CONTROL = {
    "arm": "wirehair2_raw_d12_h12_periodic",
    "arm_descriptor_sha256":
        "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11",
    "codec": "wirehair2_experiment",
}

CANDIDATE_REALIZED_CELL_ZERO = {
    "d12-h11-periodic":
        "156e7063353cde2ad42b0ea517af0f5c941feed497e60d931e2e5ee08878717c",
    "d12-h13-periodic":
        "ea2fe8cc925a523a11def53fd10b1749d496987df23c573592c7144adaa7b225",
    "d13-h12-periodic":
        "945a94dfbe9b18d66d8da38c3555728b40766fc1c5869da8193ddda8b91794c3",
}

CANDIDATE_DISCRIMINATORS = {
    # These cells distinguish each raw structure from the common raw D12/H12
    # control under the same uniform seed schedule.
    "d12-h11-periodic": {
        "cell": 5,
        "K": 8,
        "control": ("success", 0),
        "candidate": ("construct_failed", None),
    },
    "d12-h13-periodic": {
        "cell": 92,
        "K": 4,
        "control": ("success", 0),
        "candidate": ("success", 1),
    },
    "d13-h12-periodic": {
        "cell": 1,
        "K": 3,
        "control": ("construct_failed", None),
        "candidate": ("success", 0),
    },
}

TIMING_CANDIDATE_RECOVERY_GOLDENS = (
    (1, 3, "construct_failed", None),
    (4, 6, "success", 0),
    (3, 5, "success", 0),
)

CONTROL_REALIZED_CELL_ZERO = {
    "wirehair2_head":
        "83808aa3698d36652c6f2abe16d1d02b8710c9e3eac1ec93ea0928cb209b6d68",
    "wirehair1":
        "bd17e847c7cf9108e7d135c09a41878bd96afa760b189141df046eb1ea1da183",
}

RAW_CONTROL_REALIZED_CELL_ZERO = (
    "e6cbc57f44a58f7f7a3b1a907ab0749cd27b5e8e4c5946cfb3f262fa67ce64d6"
)

TIMING_CANDIDATE = {
    "arm": "wirehair2_dense_two07_basis_v1",
    "arm_descriptor_sha256":
        "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388",
    "codec": "wirehair2_experiment",
}

NATIVE_RECORD_FIELDS = frozenset((
    "schema", "ordinal", "cpu", "worker_pid", "started_monotonic_ns",
    "finished_monotonic_ns", "worker_process_start_ticks",
    "worker_binary_sha256", "message_sha256", "work_sha256", "payload",
))
TIMING_QUALIFICATION_PAYLOAD_FIELDS = frozenset((
    "base_cell_sha256", "candidate_count", "cell_sha256",
    "loss_retry_offset", "loss_seed", "ordinal", "packet_count",
    "trace_sha256", "wirehair1_decoded_extra", "wirehair1_outcome",
    "wirehair2_head_decoded_extra", "wirehair2_head_outcome",
))
TIMING_PAYLOAD_FIELDS = frozenset((
    "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
    "replicate", "base_seed_attempt", "base_loss_seed", "base_cell_sha256",
    "loss_retry_offset", "loss_seed", "fixed_received_overhead",
    "receive_overhead_cap", "invocations_per_slot", "interleave_policy",
    "panel_kind", "scope", "left_arm", "right_arm", "order",
    "left_outcome", "right_outcome", "left_decoded_extra",
    "right_decoded_extra", "elapsed_ns", "cell_sha256", "trace_sha256",
    "left_binary_sha256", "right_binary_sha256",
    "left_arm_descriptor_sha256", "right_arm_descriptor_sha256",
    "left_construction_attempt", "right_construction_attempt",
    "left_realized_construction_sha256",
    "right_realized_construction_sha256", "left_repair_map_sha256",
    "right_repair_map_sha256",
))

TIMING_BASE_CELL_ZERO = {
    "K": 8,
    "band": "2-100",
    "base_loss_seed": "0x2d0f28c7e7e786b2",
    "base_seed_attempt": 0,
    "block_bytes": 64,
    "fixed_received_overhead": 4,
    "interleave_policy": "self-counterbalanced-repeat-major-v1",
    "invocations_per_slot": 8192,
    "loss_ppm": 100000,
    "phase": "development",
    "receive_overhead_cap": 256,
    "replicate": 0,
    "schedule": "iid",
}
TIMING_BASE_CELL_ZERO_JSON = (
    '{"K":8,"band":"2-100","base_loss_seed":"0x2d0f28c7e7e786b2",'
    '"base_seed_attempt":0,"block_bytes":64,"fixed_received_overhead":4,'
    '"interleave_policy":"self-counterbalanced-repeat-major-v1",'
    '"invocations_per_slot":8192,"loss_ppm":100000,'
    '"phase":"development","receive_overhead_cap":256,"replicate":0,'
    '"schedule":"iid"}'
)
TIMING_BASE_CELL_ZERO_SHA256 = (
    "158931044b779851eca21bd3112bce79877e6e59f5896aebf698466789c47814"
)
TIMING_SOURCE_MESSAGE_SHA256 = (
    "359d10cd28719881523fd74723b803ac067830be9233c2f75325a6bd0a8324f3"
)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def sha256_json(value):
    return hashlib.sha256(canonical_json(value).encode("ascii")).hexdigest()


def qualified_timing_cell(retry_offset):
    cell = dict(TIMING_BASE_CELL_ZERO)
    cell["base_cell_sha256"] = TIMING_BASE_CELL_ZERO_SHA256
    cell["loss_retry_offset"] = retry_offset
    base_seed = int(TIMING_BASE_CELL_ZERO["base_loss_seed"], 16)
    cell["loss_seed"] = "0x{:016x}".format(
        (base_seed + retry_offset * LOSS_RETRY_STRIDE) & UINT64_MASK)
    return cell


def native_work_sha256(evidence_kind, ordinal, cell_sha256):
    return sha256_json({
        "cell_sha256": cell_sha256,
        "evidence_kind": evidence_kind,
        "ordinal": ordinal,
        "phase": "development",
        "schema": NATIVE_WORK_SCHEMA,
    })


def generic_realized_sha256(arm, K, block_bytes, construction_attempt):
    return sha256_json({
        "K": K,
        "arm_descriptor_sha256": arm["arm_descriptor_sha256"],
        "block_bytes": block_bytes,
        "codec": arm["codec"],
        "construction_attempt": construction_attempt,
        "schema": REALIZED_SCHEMA,
    })


def splitmix64_next(state):
    state = (state + LOSS_RETRY_STRIDE) & UINT64_MASK
    value = state
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & UINT64_MASK
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & UINT64_MASK
    return state, value ^ (value >> 31)


def deterministic_timing_message_sha256(base_cell):
    """Independently reproduce SourceSeedFromCellJson plus source bytes."""
    base_json = canonical_json(base_cell).encode("ascii")
    source_digest = hashlib.sha256(SOURCE_DOMAIN + base_json).digest()
    source_seed = int.from_bytes(source_digest[:8], "big")
    K = base_cell["K"]
    block_bytes = base_cell["block_bytes"]
    state = source_seed ^ ((K * 0xD6E8FEB86659FD93) & UINT64_MASK) ^ \
        ((block_bytes * 0xA0761D6478BD642F) & UINT64_MASK)
    source = bytearray()
    while len(source) < K * block_bytes:
        state, word = splitmix64_next(state)
        source.extend(word.to_bytes(8, "little"))
    return hashlib.sha256(source[:K * block_bytes]).hexdigest()


def canonical_descriptor(arm, codec, transform):
    value = {
        "arm": arm,
        "codec": codec,
        "equation_transform": transform,
        "schema": DESCRIPTOR_SCHEMA,
    }
    return json.dumps(value, separators=(",", ":"))


def candidate_arm_record(candidate):
    """Description v3 binds raw structure and seed policy independently."""
    return {
        "arm": candidate["arm"],
        "arm_descriptor_sha256": candidate["arm_descriptor_sha256"],
        "codec": candidate["codec"],
        "construction_seed_basis": RAW_SEED_BASIS,
        "dense_anchor_layout": candidate.get(
            "dense_anchor_layout",
            "two07" if candidate["arm"] == TIMING_CANDIDATE["arm"]
            else "disabled"),
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
    }


def candidate_control_record(control):
    result = dict(control)
    if control["codec"] == "wirehair1":
        result["construction_seed_basis"] = "not-applicable"
    else:
        result["construction_seed_basis"] = "production-profile-v1"
    result["dense_anchor_layout"] = \
        "not-applicable" if control["codec"] == "wirehair1" else "disabled"
    result["seed_schedule_sha256"] = ZERO_SHA256
    return result


def raw_realized_object(payload, codec):
    return {
        "K": payload["K"],
        "arm": payload["arm"],
        "arm_descriptor_sha256": payload["arm_descriptor_sha256"],
        "binary_dense_rows": payload["binary_dense_rows"],
        "block_bytes": payload["block_bytes"],
        "codec": codec,
        "construction_seed_basis": payload["construction_seed_basis"],
        "dense_anchor_layout": payload["dense_anchor_layout"],
        "dense_identity_corner": payload["dense_identity_corner"],
        "effective_packet_seed": payload["effective_packet_seed"],
        "effective_precode_seed": payload["effective_precode_seed"],
        "gf256_heavy_rows": payload["gf256_heavy_rows"],
        "heavy_family": payload["heavy_family"],
        "mix_count": payload["mix_count"],
        "packet_attempt": payload["packet_attempt"],
        "precode_attempt": payload["precode_attempt"],
        "schema": RAW_REALIZED_SCHEMA,
        "seed_schedule_sha256": payload["seed_schedule_sha256"],
        "source_hits": payload["source_hits"],
        "staircase": payload["staircase"],
    }


def raw_realized_sha256(payload, codec):
    canonical = json.dumps(
        raw_realized_object(payload, codec),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


class ContractWorkerCliTest(unittest.TestCase):
    worker = None

    def run_worker(self, *arguments, stdin=None):
        assert self.worker is not None
        return subprocess.run(
            [str(self.worker), *arguments],
            input=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
            check=False,
        )

    def description(self, *arguments):
        result = self.run_worker(*arguments)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertTrue(result.stdout.endswith("\n"))
        self.assertEqual(result.stdout.count("\n"), 1)
        return json.loads(result.stdout), result.stdout[:-1]

    def assert_exact_description(self, value, raw, arms, binary_sha256,
                                 source_git_commit,
                                 schema=DESCRIPTION_SCHEMA):
        expected = {
            "arms": arms,
            "binary_sha256": binary_sha256,
            "schema": schema,
            "source_git_commit": source_git_commit,
        }
        self.assertEqual(value, expected)
        self.assertEqual(raw, json.dumps(expected, separators=(",", ":")))

    def test_descriptions_are_canonical_closed_and_unique(self):
        baseline, baseline_raw = self.description("--describe")
        executable_hash = hashlib.sha256(self.worker.read_bytes()).hexdigest()
        self.assertEqual(baseline["binary_sha256"], executable_hash)
        self.assertRegex(baseline["binary_sha256"], SHA256)
        self.assertRegex(baseline["source_git_commit"], GIT_COMMIT)
        self.assert_exact_description(
            baseline,
            baseline_raw,
            CONTROLS + [TIMING_CANDIDATE],
            executable_hash,
            baseline["source_git_commit"],
        )

        observed_hashes = {
            arm["arm_descriptor_sha256"]
            for arm in CONTROLS + [TIMING_CANDIDATE]
        }
        self.assertEqual(RAW_SEED_BASIS, "uniform-raw-v1")
        self.assertEqual(
            hashlib.sha256(
                RAW_SEED_SCHEDULE_CANONICAL.encode("ascii")).hexdigest(),
            RAW_SEED_SCHEDULE_SHA256,
        )
        raw_transform = "d12-h12-periodic"
        raw_descriptor = canonical_descriptor(
            RAW_CONTROL["arm"], RAW_CONTROL["codec"], raw_transform)
        raw_descriptor_hash = hashlib.sha256(
            raw_descriptor.encode("ascii")).hexdigest()
        self.assertEqual(
            RAW_CONTROL["arm_descriptor_sha256"], raw_descriptor_hash)
        self.assertNotIn(raw_descriptor_hash, observed_hashes)
        observed_hashes.add(raw_descriptor_hash)

        for candidate_id, candidate in CANDIDATES.items():
            with self.subTest(candidate=candidate_id):
                self.assertEqual(
                    candidate_id,
                    "d{}-h{}-periodic".format(
                        candidate["dense_rows"], candidate["heavy_rows"]),
                )
                self.assertEqual(
                    candidate["heavy_family"], "periodic-cauchy")
                description, raw = self.description(
                    "--describe-recovery-candidate", candidate_id)
                self.assert_exact_description(
                    description,
                    raw,
                    [candidate_control_record(control) for control in CONTROLS]
                    + [candidate_arm_record(RAW_CONTROL),
                       candidate_arm_record(candidate)],
                    executable_hash,
                    baseline["source_git_commit"],
                    RAW_DESCRIPTION_SCHEMA,
                )
                transform = candidate_id
                descriptor = canonical_descriptor(
                    candidate["arm"], candidate["codec"], transform)
                descriptor_hash = hashlib.sha256(
                    descriptor.encode("ascii")).hexdigest()
                self.assertEqual(
                    candidate["arm_descriptor_sha256"], descriptor_hash)
                self.assertNotIn(descriptor_hash, observed_hashes)
                observed_hashes.add(descriptor_hash)

        two07_description, two07_raw = self.description(
            "--describe-recovery-candidate", "two07")
        self.assert_exact_description(
            two07_description,
            two07_raw,
            [candidate_control_record(control) for control in CONTROLS] +
            [candidate_arm_record(RAW_CONTROL),
             candidate_arm_record(TIMING_CANDIDATE)],
            executable_hash,
            baseline["source_git_commit"],
            RAW_DESCRIPTION_SCHEMA,
        )
        self.assertEqual(
            two07_description["arms"][3]["arm_descriptor_sha256"],
            baseline["arms"][2]["arm_descriptor_sha256"])

    def test_structure_only_timing_proxy_witness_is_closed_and_canonical(self):
        result = self.run_worker("--emit-timing-proxy-witness")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        self.assertTrue(result.stdout.endswith("\n"))
        witness = json.loads(result.stdout)
        self.assertEqual(
            result.stdout[:-1],
            json.dumps(witness, sort_keys=True, separators=(",", ":")))
        self.assertEqual(witness["schema"], TIMING_PROXY_SCHEMA)
        self.assertEqual(
            witness["proof_scope"],
            "d12-disabled-structure-under-production-timing-seed-policy-v1")
        self.assertEqual(witness["evidence_phase"], "development")
        self.assertEqual(witness["construction_attempts"], [0])
        self.assertEqual(
            witness["raw_recovery_reference_arm_descriptor_sha256"],
            RAW_CONTROL["arm_descriptor_sha256"])
        self.assertEqual(witness["raw_recovery_seed_basis"], RAW_SEED_BASIS)
        self.assertEqual(
            witness["raw_recovery_seed_schedule_sha256"],
            RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(
            witness["timing_seed_policy_arms"],
            ["wirehair2_head", TIMING_CANDIDATE["arm"]])
        coordinates = [
            {"K": K, "block_bytes": block_bytes,
             "construction_attempt": 0}
            for block_bytes in (64, 1280)
            for K in (8, 32, 100, 128, 512, 1000, 2048, 5000,
                      8192, 20000, 32768, 64000)
        ]
        self.assertEqual(
            witness["witness_domain_sha256"],
            hashlib.sha256(
                TIMING_PROXY_DOMAIN +
                canonical_json(coordinates).encode("ascii")).hexdigest())
        self.assertEqual(len(witness["cells"]), 24)
        self.assertEqual(witness["cells"][0], TIMING_PROXY_K8_B64_GOLDEN)
        self.assertEqual(
            hashlib.sha256(canonical_json(witness["cells"]).encode(
                "ascii")).hexdigest(),
            TIMING_PROXY_CELLS_SHA256)
        for cell, coordinate in zip(witness["cells"], coordinates):
            self.assertEqual(
                {key: cell[key] for key in coordinate}, coordinate)
            self.assertTrue(cell["seeds_differ"])
            self.assertNotEqual(
                cell["production_timing_precode_seed"],
                cell["raw_recovery_precode_seed"])
            self.assertNotEqual(
                cell["production_timing_packet_seed"],
                cell["raw_recovery_packet_seed"])
            for field in (
                    "normalized_structure_sha256",
                    "production_timing_configuration_sha256",
                    "production_timing_equation_system_sha256"):
                self.assertRegex(cell[field], SHA256)

    def test_unknown_candidate_and_noncanonical_arity_reject(self):
        for arguments in (
            ("--describe-recovery-candidate", "unknown"),
            ("--describe-recovery-candidate",),
            ("--describe-recovery-candidate", "d12-h11-periodic", "extra"),
            ("--recovery-candidate-worker", "unknown", "0"),
            ("--recovery-candidate-worker", "d12-h11-periodic"),
        ):
            with self.subTest(arguments=arguments):
                result = self.run_worker(*arguments)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(result.stdout, "")
                self.assertIn("usage: wirehair_wh2_contract_worker", result.stderr)

    @unittest.skipUnless(
        hasattr(os, "sched_getaffinity"), "Linux CPU affinity is required")
    def test_candidate_worker_is_recovery_only_and_runs_each_transform(self):
        affinity = os.sched_getaffinity(0)
        self.assertTrue(affinity)
        cpu = str(min(affinity))

        clean = self.run_worker(
            "--recovery-candidate-worker",
            "d12-h11-periodic",
            cpu,
            stdin="Q\n",
        )
        self.assertEqual(clean.returncode, 0, clean.stderr)
        self.assertEqual(clean.stdout, "")
        self.assertEqual(clean.stderr, "")

        for command in ("L 0 0\n", "T 0 0\n"):
            with self.subTest(rejected_command=command.strip()):
                rejected = self.run_worker(
                    "--recovery-candidate-worker",
                    "d12-h11-periodic",
                    cpu,
                    stdin=command,
                )
                self.assertEqual(rejected.returncode, 1)
                self.assertEqual(rejected.stdout, "")
                self.assertIn(
                    "recovery candidate worker rejects qualification/timing "
                    "jobs",
                    rejected.stderr,
                )

        baseline = self.run_worker(
            "--worker",
            cpu,
            stdin="R 0 0\nR 0 1\nQ\n",
        )
        self.assertEqual(baseline.returncode, 0, baseline.stderr)
        self.assertEqual(baseline.stderr, "")
        baseline_records = [
            json.loads(line) for line in baseline.stdout.splitlines()]
        self.assertEqual(len(baseline_records), 2)
        baseline_realized = {
            record["payload"]["arm"]:
                record["payload"]["realized_construction_sha256"]
            for record in baseline_records
        }
        self.assertEqual(baseline_realized, CONTROL_REALIZED_CELL_ZERO)

        timing_candidate = self.run_worker(
            "--worker",
            cpu,
            stdin="".join(
                "R {} 2\n".format(cell)
                for cell, _, _, _ in TIMING_CANDIDATE_RECOVERY_GOLDENS) +
                "Q\n",
        )
        self.assertEqual(
            timing_candidate.returncode, 0, timing_candidate.stderr)
        self.assertEqual(timing_candidate.stderr, "")
        timing_candidate_records = [
            json.loads(line)
            for line in timing_candidate.stdout.splitlines()]
        self.assertEqual(
            len(timing_candidate_records),
            len(TIMING_CANDIDATE_RECOVERY_GOLDENS))
        self.assertEqual(
            [(record["payload"]["K"], record["payload"]["outcome"])
             for record in timing_candidate_records],
            [(K, outcome)
             for _, K, outcome, _ in TIMING_CANDIDATE_RECOVERY_GOLDENS],
        )
        for record, (_, K, outcome, decoded_extra) in zip(
                timing_candidate_records,
                TIMING_CANDIDATE_RECOVERY_GOLDENS):
            payload = record["payload"]
            self.assertEqual(payload["arm"], TIMING_CANDIDATE["arm"])
            self.assertEqual(
                payload["arm_descriptor_sha256"],
                TIMING_CANDIDATE["arm_descriptor_sha256"])
            self.assertEqual(payload["K"], K)
            self.assertEqual(payload["outcome"], outcome)
            self.assertEqual(payload["decoded_extra"], decoded_extra)
            self.assertEqual(
                payload["construction_attempt"],
                payload["base_seed_attempt"])
            self.assertEqual(
                payload["realized_construction_sha256"],
                generic_realized_sha256(
                    TIMING_CANDIDATE, payload["K"],
                    payload["block_bytes"],
                    payload["construction_attempt"]))
            self.assertEqual(payload["repair_map_sha256"], ZERO_SHA256)

        candidate_realized = set()
        for candidate_id, candidate in CANDIDATES.items():
            with self.subTest(candidate=candidate_id):
                discriminator = CANDIDATE_DISCRIMINATORS[candidate_id]
                discriminator_cell = discriminator["cell"]
                result = self.run_worker(
                    "--recovery-candidate-worker",
                    candidate_id,
                    cpu,
                    stdin=(
                        "R 0 0\nR 0 1\nR 0 2\nR 0 3\n"
                        "R {0} 2\nR {0} 3\nQ\n"
                    ).format(discriminator_cell),
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                lines = result.stdout.splitlines()
                self.assertEqual(len(lines), 6)
                records = [json.loads(line) for line in lines]
                self.assertEqual(
                    lines,
                    [json.dumps(record, sort_keys=True, separators=(",", ":"))
                     for record in records],
                )
                for control_record in records[:2]:
                    self.assertEqual(control_record["schema"], RECOVERY_SCHEMA)
                    self.assertNotIn(
                        "construction_seed_basis", control_record["payload"])
                control_realized = {
                    record["payload"]["arm"]:
                        record["payload"]["realized_construction_sha256"]
                    for record in records[:2]
                }
                self.assertEqual(control_realized, baseline_realized)

                raw_control_record = records[2]
                self.assertEqual(
                    raw_control_record["schema"], RAW_RECOVERY_SCHEMA)
                self.assertEqual(raw_control_record["ordinal"], 2)
                raw_control_payload = raw_control_record["payload"]
                self.assertEqual(
                    raw_control_payload["arm"], RAW_CONTROL["arm"])
                self.assertEqual(
                    raw_control_payload["arm_descriptor_sha256"],
                    RAW_CONTROL["arm_descriptor_sha256"],
                )
                self.assertEqual(
                    raw_control_payload["construction_attempt"],
                    raw_control_payload["base_seed_attempt"],
                )
                self.assertEqual(
                    raw_control_payload["precode_attempt"],
                    raw_control_payload["construction_attempt"],
                )
                self.assertEqual(
                    raw_control_payload["packet_attempt"],
                    raw_control_payload["construction_attempt"],
                )
                self.assertEqual(
                    raw_control_payload["construction_seed_basis"],
                    RAW_SEED_BASIS,
                )
                self.assertEqual(
                    raw_control_payload["seed_schedule_sha256"],
                    RAW_SEED_SCHEDULE_SHA256,
                )
                self.assertEqual(raw_control_payload["binary_dense_rows"], 12)
                self.assertEqual(raw_control_payload["gf256_heavy_rows"], 12)
                self.assertEqual(
                    raw_control_payload["dense_anchor_layout"], "disabled")
                self.assertEqual(
                    raw_control_payload["heavy_family"], "periodic-cauchy")
                self.assertRegex(
                    raw_control_payload["effective_precode_seed"],
                    r"^0x[0-9a-f]{16}$",
                )
                self.assertRegex(
                    raw_control_payload["effective_packet_seed"],
                    r"^0x[0-9a-f]{8}$",
                )
                self.assertEqual(
                    raw_control_payload["realized_construction_sha256"],
                    RAW_CONTROL_REALIZED_CELL_ZERO,
                )
                self.assertEqual(
                    raw_realized_sha256(
                        raw_control_payload, RAW_CONTROL["codec"]),
                    RAW_CONTROL_REALIZED_CELL_ZERO,
                )
                added_fields = {
                    "construction_seed_basis", "seed_schedule_sha256",
                    "precode_attempt", "packet_attempt",
                    "effective_precode_seed", "effective_packet_seed",
                    "staircase", "binary_dense_rows", "gf256_heavy_rows",
                    "source_hits", "dense_identity_corner", "heavy_family",
                    "mix_count", "dense_anchor_layout",
                }
                self.assertEqual(
                    set(raw_control_payload),
                    set(records[0]["payload"]) | added_fields,
                )

                record = records[3]
                self.assertEqual(record["schema"], RAW_RECOVERY_SCHEMA)
                self.assertEqual(record["ordinal"], 3)
                payload = record["payload"]
                self.assertEqual(payload["arm"], candidate["arm"])
                self.assertEqual(
                    payload["arm_descriptor_sha256"],
                    candidate["arm_descriptor_sha256"],
                )
                self.assertEqual(
                    payload["construction_attempt"],
                    payload["base_seed_attempt"],
                )
                self.assertEqual(payload["precode_attempt"], 0)
                self.assertEqual(payload["packet_attempt"], 0)
                self.assertEqual(
                    payload["effective_precode_seed"],
                    "0x487468302aad7105",
                )
                self.assertEqual(
                    payload["effective_packet_seed"], "0x4ec72102")
                self.assertEqual(
                    payload["binary_dense_rows"], candidate["dense_rows"])
                self.assertEqual(
                    payload["gf256_heavy_rows"], candidate["heavy_rows"])
                self.assertEqual(payload["dense_anchor_layout"], "disabled")
                self.assertEqual(payload["source_hits"], 2)
                self.assertEqual(payload["staircase"], 2)
                self.assertFalse(payload["dense_identity_corner"])
                self.assertEqual(payload["mix_count"], 3)
                self.assertEqual(
                    payload["heavy_family"], candidate["heavy_family"])
                realized = payload["realized_construction_sha256"]
                self.assertEqual(
                    realized, CANDIDATE_REALIZED_CELL_ZERO[candidate_id])
                self.assertEqual(
                    raw_realized_sha256(payload, candidate["codec"]), realized)
                self.assertNotIn(realized, candidate_realized)
                self.assertNotIn(realized, baseline_realized.values())
                self.assertNotEqual(realized, RAW_CONTROL_REALIZED_CELL_ZERO)
                candidate_realized.add(realized)

                realized_object = raw_realized_object(
                    payload, candidate["codec"])
                for field, original in realized_object.items():
                    with self.subTest(candidate=candidate_id, field=field):
                        changed = dict(realized_object)
                        if isinstance(original, bool):
                            changed[field] = not original
                        elif isinstance(original, int):
                            changed[field] = original + 1
                        else:
                            changed[field] = original + "-tampered"
                        changed_hash = hashlib.sha256(json.dumps(
                            changed,
                            sort_keys=True,
                            separators=(",", ":"),
                        ).encode("ascii")).hexdigest()
                        self.assertNotEqual(changed_hash, realized)

                raw_control_discriminator = records[4]["payload"]
                candidate_discriminator = records[5]["payload"]
                self.assertEqual(records[4]["schema"], RAW_RECOVERY_SCHEMA)
                self.assertEqual(records[5]["schema"], RAW_RECOVERY_SCHEMA)
                for discriminator_record, arm_codec in (
                        (records[4], RAW_CONTROL["codec"]),
                        (records[5], candidate["codec"])):
                    discriminator_payload = discriminator_record["payload"]
                    attempt = discriminator_payload["base_seed_attempt"]
                    self.assertEqual(
                        discriminator_payload["precode_attempt"], attempt)
                    self.assertEqual(
                        discriminator_payload["packet_attempt"], attempt)
                    self.assertEqual(
                        discriminator_payload["effective_precode_seed"],
                        "0x{:016x}".format(
                            (0x487468302AAD7105
                             + attempt * 0x9E3779B97F4A7C15)
                            & ((1 << 64) - 1)),
                    )
                    self.assertEqual(
                        discriminator_payload["effective_packet_seed"],
                        "0x{:08x}".format(
                            (0x4EC72102 + attempt * 0x9E3779B9)
                            & ((1 << 32) - 1)),
                    )
                    self.assertEqual(
                        raw_realized_sha256(
                            discriminator_payload, arm_codec),
                        discriminator_payload[
                            "realized_construction_sha256"],
                    )
                self.assertEqual(
                    raw_control_discriminator["arm"], RAW_CONTROL["arm"])
                self.assertEqual(
                    candidate_discriminator["arm"], candidate["arm"])
                self.assertEqual(
                    raw_control_discriminator["K"], discriminator["K"])
                self.assertEqual(
                    candidate_discriminator["K"], discriminator["K"])
                self.assertEqual(
                    (raw_control_discriminator["outcome"],
                     raw_control_discriminator["decoded_extra"]),
                    discriminator["control"],
                )
                self.assertEqual(
                    (candidate_discriminator["outcome"],
                     candidate_discriminator["decoded_extra"]),
                    discriminator["candidate"],
                )

    @unittest.skipUnless(
        hasattr(os, "sched_getaffinity"), "Linux CPU affinity is required")
    def test_timing_qualification_records_bind_exact_retry_evidence(self):
        affinity = os.sched_getaffinity(0)
        self.assertTrue(affinity)
        cpu = str(min(affinity))
        result = self.run_worker(
            "--worker", cpu, stdin="L 0 0\nL 0 1\nQ\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2)
        records = [json.loads(line) for line in lines]
        self.assertEqual(
            lines, [canonical_json(record) for record in records])

        expected_message = hashlib.sha256(
            TIMING_QUALIFICATION_MESSAGE_DOMAIN +
            TIMING_BASE_CELL_ZERO_JSON.encode("ascii")
        ).hexdigest()
        self.assertEqual(
            canonical_json(TIMING_BASE_CELL_ZERO), TIMING_BASE_CELL_ZERO_JSON)
        self.assertEqual(
            sha256_json(TIMING_BASE_CELL_ZERO),
            TIMING_BASE_CELL_ZERO_SHA256)
        self.assertEqual(
            expected_message,
            "ed140dcfb6a14da3adc461a7717d68318aacbdca29acdc1fe1464fa36a59509c")

        expected_payloads = [
            {
                "base_cell_sha256": TIMING_BASE_CELL_ZERO_SHA256,
                "candidate_count": 305,
                "cell_sha256":
                    "e85871b9c09ac946f6c8f6e0e077bd758d4862bfb5dddd73f6ee05e3ba35f2d7",
                "loss_retry_offset": 0,
                "loss_seed": "0x2d0f28c7e7e786b2",
                "ordinal": 0,
                "packet_count": 264,
                "trace_sha256":
                    "768c7e52d3bfc315b9d9eb44287766f9401904062a5595be439f8ba6b20c90aa",
                "wirehair1_decoded_extra": 0,
                "wirehair1_outcome": "success",
                "wirehair2_head_decoded_extra": 4,
                "wirehair2_head_outcome": "success",
            },
            {
                "base_cell_sha256": TIMING_BASE_CELL_ZERO_SHA256,
                "candidate_count": 302,
                "cell_sha256":
                    "4c216d6ca9416585f777980fc1ffd0a9f162644606d5f0468accec353470cf37",
                "loss_retry_offset": 1,
                "loss_seed": "0xcb46a281673202c7",
                "ordinal": 0,
                "packet_count": 264,
                "trace_sha256":
                    "3b1b36be4d0f88d0da4dd6ca02b1fecf43c7909ee63d1b80df40443422a2fdc9",
                "wirehair1_decoded_extra": 0,
                "wirehair1_outcome": "success",
                "wirehair2_head_decoded_extra": 4,
                "wirehair2_head_outcome": "success",
            },
        ]
        expected_work = (
            "27ac11c0f5dddb46078f3f467684cd69a1cec7802c758a8449a2c2793f330adc",
            "864e59732fc73fbb61e241aa04a31d24f971bdd26516760c4beaea4bdd289d16",
        )
        executable_hash = hashlib.sha256(self.worker.read_bytes()).hexdigest()
        for retry, (record, payload, work_sha256) in enumerate(zip(
                records, expected_payloads, expected_work)):
            with self.subTest(retry=retry):
                self.assertEqual(set(record), NATIVE_RECORD_FIELDS)
                self.assertEqual(record["schema"], TIMING_QUALIFICATION_SCHEMA)
                self.assertEqual(record["ordinal"], retry)
                self.assertEqual(record["cpu"], int(cpu))
                self.assertEqual(record["worker_binary_sha256"], executable_hash)
                self.assertEqual(record["message_sha256"], expected_message)
                self.assertEqual(record["payload"], payload)
                self.assertEqual(
                    set(record["payload"]),
                    TIMING_QUALIFICATION_PAYLOAD_FIELDS)
                self.assertEqual(
                    record["payload"]["cell_sha256"],
                    sha256_json(qualified_timing_cell(retry)))
                self.assertEqual(
                    native_work_sha256(
                        "timing_qualification", retry,
                        record["payload"]["cell_sha256"]),
                    work_sha256)
                self.assertEqual(record["work_sha256"], work_sha256)
                self.assertGreater(record["worker_pid"], 0)
                self.assertGreater(record["worker_process_start_ticks"], 0)
                self.assertLessEqual(
                    record["started_monotonic_ns"],
                    record["finished_monotonic_ns"])
        self.assertEqual(records[0]["worker_pid"], records[1]["worker_pid"])
        self.assertEqual(
            records[0]["message_sha256"], records[1]["message_sha256"])
        self.assertNotEqual(
            records[0]["payload"]["cell_sha256"],
            records[1]["payload"]["cell_sha256"])
        self.assertNotEqual(
            records[0]["payload"]["trace_sha256"],
            records[1]["payload"]["trace_sha256"])
        self.assertNotEqual(records[0]["work_sha256"], records[1]["work_sha256"])

    @unittest.skipUnless(
        hasattr(os, "sched_getaffinity"), "Linux CPU affinity is required")
    def test_timing_v4_retry_changes_trace_but_not_base_source(self):
        affinity = os.sched_getaffinity(0)
        self.assertTrue(affinity)
        cpu = str(min(affinity))
        # Packed T item is panel_index * 256 + retry_offset.  Both commands
        # therefore select panel zero while changing only retry zero to one.
        result = self.run_worker(
            "--worker", cpu, stdin="T 0 0\nT 0 1\nQ\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2)
        records = [json.loads(line) for line in lines]
        self.assertEqual(lines, [canonical_json(record) for record in records])

        expected_source_message = deterministic_timing_message_sha256(
            TIMING_BASE_CELL_ZERO)
        self.assertEqual(expected_source_message, TIMING_SOURCE_MESSAGE_SHA256)
        expected_cells = (
            "e85871b9c09ac946f6c8f6e0e077bd758d4862bfb5dddd73f6ee05e3ba35f2d7",
            "4c216d6ca9416585f777980fc1ffd0a9f162644606d5f0468accec353470cf37",
        )
        expected_traces = (
            "768c7e52d3bfc315b9d9eb44287766f9401904062a5595be439f8ba6b20c90aa",
            "3b1b36be4d0f88d0da4dd6ca02b1fecf43c7909ee63d1b80df40443422a2fdc9",
        )
        expected_loss_seeds = (
            "0x2d0f28c7e7e786b2", "0xcb46a281673202c7")
        expected_work = (
            "69efde6d79ea571b5b8e4dea365fb71609fcba4908b800aaf1496b5a1aef41f8",
            "f1a4cf9c23c158a28508e09b37d4a6dd28cb718470391bf5b904f61edf4acf46",
        )
        executable_hash = hashlib.sha256(self.worker.read_bytes()).hexdigest()
        for retry, record in enumerate(records):
            with self.subTest(retry=retry):
                payload = record["payload"]
                self.assertEqual(set(record), NATIVE_RECORD_FIELDS)
                self.assertEqual(set(payload), TIMING_PAYLOAD_FIELDS)
                self.assertEqual(record["schema"], TIMING_SCHEMA)
                # The result ordinal is cell * 11 + panel, independent of retry.
                self.assertEqual(record["ordinal"], 0)
                self.assertEqual(record["cpu"], int(cpu))
                self.assertEqual(record["worker_binary_sha256"], executable_hash)
                self.assertEqual(record["message_sha256"], expected_source_message)
                self.assertEqual(payload["K"], 8)
                self.assertEqual(payload["phase"], "development")
                self.assertEqual(payload["band"], "2-100")
                self.assertEqual(payload["block_bytes"], 64)
                self.assertEqual(payload["loss_ppm"], 100000)
                self.assertEqual(payload["schedule"], "iid")
                self.assertEqual(payload["replicate"], 0)
                self.assertEqual(payload["base_seed_attempt"], 0)
                self.assertEqual(
                    payload["base_cell_sha256"], TIMING_BASE_CELL_ZERO_SHA256)
                self.assertEqual(
                    payload["base_loss_seed"],
                    TIMING_BASE_CELL_ZERO["base_loss_seed"])
                self.assertEqual(payload["loss_retry_offset"], retry)
                self.assertEqual(payload["loss_seed"], expected_loss_seeds[retry])
                self.assertEqual(payload["cell_sha256"], expected_cells[retry])
                self.assertEqual(
                    payload["cell_sha256"],
                    sha256_json(qualified_timing_cell(retry)))
                self.assertEqual(payload["trace_sha256"], expected_traces[retry])
                self.assertEqual(payload["fixed_received_overhead"], 4)
                self.assertEqual(payload["receive_overhead_cap"], 256)
                self.assertEqual(
                    payload["interleave_policy"],
                    "self-counterbalanced-repeat-major-v1")
                self.assertEqual(payload["invocations_per_slot"], 8192)
                self.assertEqual(payload["order"], "BAAB")
                self.assertEqual(payload["panel_kind"], "AA")
                self.assertEqual(payload["scope"], "decoder_solve")
                self.assertEqual(payload["left_arm"], "wirehair2_head")
                self.assertEqual(payload["right_arm"], "wirehair2_head")
                self.assertEqual(
                    payload["left_arm_descriptor_sha256"],
                    CONTROLS[0]["arm_descriptor_sha256"])
                self.assertEqual(
                    payload["right_arm_descriptor_sha256"],
                    CONTROLS[0]["arm_descriptor_sha256"])
                self.assertEqual(payload["left_construction_attempt"], 0)
                self.assertEqual(payload["right_construction_attempt"], 0)
                self.assertEqual(
                    payload["left_realized_construction_sha256"],
                    "ca6c4a450d9d0abfa024d354fe436590829bce1d0235a0ec97fb2faf3e3faf32")
                self.assertEqual(
                    payload["right_realized_construction_sha256"],
                    payload["left_realized_construction_sha256"])
                self.assertEqual(payload["left_repair_map_sha256"], ZERO_SHA256)
                self.assertEqual(payload["right_repair_map_sha256"], ZERO_SHA256)
                self.assertEqual(payload["left_outcome"], "success")
                self.assertEqual(payload["right_outcome"], "success")
                self.assertEqual(payload["left_decoded_extra"], 4)
                self.assertEqual(payload["right_decoded_extra"], 4)
                self.assertEqual(payload["left_binary_sha256"], executable_hash)
                self.assertEqual(payload["right_binary_sha256"], executable_hash)
                self.assertEqual(len(payload["elapsed_ns"]), 8)
                self.assertTrue(all(
                    type(value) is int and value > 0
                    for value in payload["elapsed_ns"]))
                self.assertEqual(
                    native_work_sha256("timing", 0, expected_cells[retry]),
                    expected_work[retry])
                self.assertEqual(record["work_sha256"], expected_work[retry])
        self.assertEqual(records[0]["message_sha256"],
                         records[1]["message_sha256"])
        self.assertNotEqual(records[0]["payload"]["loss_seed"],
                            records[1]["payload"]["loss_seed"])
        self.assertNotEqual(records[0]["payload"]["trace_sha256"],
                            records[1]["payload"]["trace_sha256"])
        self.assertNotEqual(records[0]["payload"]["cell_sha256"],
                            records[1]["payload"]["cell_sha256"])
        self.assertNotEqual(records[0]["work_sha256"],
                            records[1]["work_sha256"])

    @unittest.skipUnless(
        hasattr(os, "sched_getaffinity"), "Linux CPU affinity is required")
    def test_timing_candidate_executes_enabled_two_anchor_codec(self):
        affinity = os.sched_getaffinity(0)
        self.assertTrue(affinity)
        cpu = str(min(affinity))
        # Cell 48 is replicate-zero K=2048/B=64.  Packed items 4*256 and
        # 7*256 select the candidate decoder-solve A/A and candidate/head A/B
        # panels at retry zero.
        result = self.run_worker(
            "--worker", cpu, stdin="T 48 1024\nT 48 1792\nQ\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2)
        records = [json.loads(line) for line in lines]
        self.assertEqual(lines, [canonical_json(record) for record in records])
        expected_realized = generic_realized_sha256(
            TIMING_CANDIDATE, 2048, 64, 0)
        expected_rights = (
            ("AA", TIMING_CANDIDATE, expected_realized),
            ("AB", CONTROLS[0],
             "323d9bd9919c764e5f8a0a40148f31fb"
             "0607964f1e8505b36bb0a8121e35636f"),
        )
        for record, (panel_kind, right_arm, right_realized) in zip(
                records, expected_rights):
            payload = record["payload"]
            self.assertEqual(record["schema"], TIMING_SCHEMA)
            self.assertEqual(payload["K"], 2048)
            self.assertEqual(payload["block_bytes"], 64)
            self.assertEqual(payload["panel_kind"], panel_kind)
            self.assertEqual(payload["scope"], "decoder_solve")
            self.assertEqual(payload["left_arm"], TIMING_CANDIDATE["arm"])
            self.assertEqual(payload["right_arm"], right_arm["arm"])
            self.assertEqual(
                payload["left_arm_descriptor_sha256"],
                TIMING_CANDIDATE["arm_descriptor_sha256"])
            self.assertEqual(
                payload["right_arm_descriptor_sha256"],
                right_arm["arm_descriptor_sha256"])
            self.assertEqual(payload["left_construction_attempt"], 0)
            self.assertEqual(payload["right_construction_attempt"], 0)
            self.assertEqual(
                payload["left_realized_construction_sha256"],
                expected_realized)
            self.assertEqual(
                payload["right_realized_construction_sha256"],
                right_realized)
            self.assertEqual(payload["left_outcome"], "success")
            self.assertEqual(payload["right_outcome"], "success")
            self.assertEqual(payload["left_decoded_extra"], 4)
            self.assertEqual(payload["right_decoded_extra"], 4)
            self.assertEqual(len(payload["elapsed_ns"]), 8)
            self.assertTrue(all(
                type(value) is int and value > 0
                for value in payload["elapsed_ns"]))

    @unittest.skipUnless(
        hasattr(os, "sched_getaffinity"), "Linux CPU affinity is required")
    def test_timing_commands_reject_malformed_and_out_of_range_indexes(self):
        affinity = os.sched_getaffinity(0)
        self.assertTrue(affinity)
        cpu = str(min(affinity))
        cases = (
            ("L 0 00\n", "malformed command"),
            ("L 0 -1\n", "malformed command"),
            ("L 0 0 extra\n", "malformed command"),
            ("X 0 0\n", "malformed command"),
            ("L 2304 0\n", "qualification index is outside"),
            ("L 0 256\n", "qualification index is outside"),
            ("T 2304 0\n", "timing job index is outside"),
            # panel 11 is out of range: packed item = 11 * 256 + retry 0.
            ("T 0 2816\n", "timing job index is outside"),
        )
        for command, diagnostic in cases:
            with self.subTest(command=command.strip()):
                result = self.run_worker(
                    "--worker", cpu, stdin=command)
                self.assertEqual(result.returncode, 1)
                self.assertEqual(result.stdout, "")
                self.assertIn(diagnostic, result.stderr)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "usage: wh2_contract_worker_cli_test.py WORKER")
    worker = Path(sys.argv[1]).resolve()
    if not worker.is_file():
        raise SystemExit("worker executable does not exist: {}".format(worker))
    ContractWorkerCliTest.worker = worker
    unittest.main(argv=[sys.argv[0]], verbosity=2)


if __name__ == "__main__":
    main()
