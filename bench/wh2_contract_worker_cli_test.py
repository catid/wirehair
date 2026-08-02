#!/usr/bin/env python3
"""Focused CLI goldens for the native WH2 recovery-candidate worker."""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import unittest


DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v1"
RAW_DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v2"
DESCRIPTOR_SCHEMA = "wirehair.wh2.native-arm-descriptor.v1"
RECOVERY_SCHEMA = "wirehair.wh2.native-recovery-record.v1"
RAW_RECOVERY_SCHEMA = "wirehair.wh2.native-recovery-record.v2"
RAW_REALIZED_SCHEMA = "wirehair.wh2.raw-realized-construction.v1"
ZERO_SHA256 = "0" * 64
SHA256 = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT = re.compile(r"^[0-9a-f]{40}$")
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
        "d9b77ab20e2d06aa4a9664c85f8f020dbec35249ad1f911a916c7ccae8937a49",
    "d12-h13-periodic":
        "a89c6b491b5f36bfe18b6ff21d6f327cb45b9e0c48b0f5b63571a002571ce8bf",
    "d13-h12-periodic":
        "943fffd7354809c211d93da6303616216287eaa554400d726199bffaeccf323a",
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

GENERIC_IDENTITY_DISCRIMINATORS = ((1, 3), (4, 6), (3, 5))

CONTROL_REALIZED_CELL_ZERO = {
    "wirehair2_head":
        "83808aa3698d36652c6f2abe16d1d02b8710c9e3eac1ec93ea0928cb209b6d68",
    "wirehair1":
        "bd17e847c7cf9108e7d135c09a41878bd96afa760b189141df046eb1ea1da183",
}

RAW_CONTROL_REALIZED_CELL_ZERO = (
    "7562bed772f60e6cc3cfc7f1d9448e1c6fa88f976af003dfdfd4cc8425da43d6"
)

IDENTITY = {
    "arm": "wirehair2_identity",
    "arm_descriptor_sha256":
        "87271a16c290b8e3e14ff76fae52ca07e94ec9c31d4188b9ceee3740c299c4ac",
    "codec": "wirehair2_experiment",
}


def canonical_descriptor(arm, codec, transform):
    value = {
        "arm": arm,
        "codec": codec,
        "equation_transform": transform,
        "schema": DESCRIPTOR_SCHEMA,
    }
    return json.dumps(value, separators=(",", ":"))


def candidate_arm_record(candidate):
    """Description v2 binds raw structure and seed policy independently."""
    return {
        "arm": candidate["arm"],
        "arm_descriptor_sha256": candidate["arm_descriptor_sha256"],
        "codec": candidate["codec"],
        "construction_seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
    }


def candidate_control_record(control):
    result = dict(control)
    if control["codec"] == "wirehair1":
        result["construction_seed_basis"] = "not-applicable"
    else:
        result["construction_seed_basis"] = "production-profile-v1"
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
            CONTROLS + [IDENTITY],
            executable_hash,
            baseline["source_git_commit"],
        )

        observed_hashes = {
            arm["arm_descriptor_sha256"] for arm in CONTROLS + [IDENTITY]
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

        timing = self.run_worker(
            "--recovery-candidate-worker",
            "d12-h11-periodic",
            cpu,
            stdin="T 0 0\n",
        )
        self.assertEqual(timing.returncode, 1)
        self.assertEqual(timing.stdout, "")
        self.assertIn("recovery candidate worker rejects timing jobs",
                      timing.stderr)

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

        identity = self.run_worker(
            "--worker",
            cpu,
            stdin="".join(
                "R {} 2\n".format(cell)
                for cell, _ in GENERIC_IDENTITY_DISCRIMINATORS) + "Q\n",
        )
        self.assertEqual(identity.returncode, 0, identity.stderr)
        self.assertEqual(identity.stderr, "")
        identity_records = [
            json.loads(line) for line in identity.stdout.splitlines()]
        self.assertEqual(
            len(identity_records), len(GENERIC_IDENTITY_DISCRIMINATORS))
        self.assertEqual(
            [(record["payload"]["K"], record["payload"]["outcome"])
             for record in identity_records],
            [(K, "success")
             for _, K in GENERIC_IDENTITY_DISCRIMINATORS],
        )

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
                    "mix_count",
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
