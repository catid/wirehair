#!/usr/bin/env python3
"""Adversarial tests for the path-only WH2 architecture selector."""

from __future__ import annotations

import copy
import errno
import os
from pathlib import Path
import stat
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_select_development_architecture as subject
from test_wh2_benchmark_contract import (
    architecture_selection_receipt, write_timing_qualification,
)


class PathOnlySelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = contract_api.load_contract()
        cls.qualification_temp = tempfile.TemporaryDirectory()
        cls.qualification = write_timing_qualification(
            Path(cls.qualification_temp.name), cls.contract)[0]

    @classmethod
    def tearDownClass(cls) -> None:
        cls.qualification_temp.cleanup()

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.receipt = architecture_selection_receipt(
            self.contract, self.qualification)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_receipt(self, path: Path, value=None) -> None:
        if value is None:
            value = self.receipt
        path.write_bytes(
            (contract_api.canonical_json(value) + "\n").encode("utf-8"))

    def test_atomic_publication_reopens_and_seals_exact_receipt(self) -> None:
        output = self.root / "development-selection.json"
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt):
            handle = subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertIs(type(handle), contract_api.ArchitectureSelection)
        self.assertEqual(handle.selection_sha256,
                         self.receipt["selection_sha256"])
        self.assertEqual(
            subject._published_receipt(self.contract, output)[0], self.receipt)
        with self.assertRaises(subject.SelectionError):
            subject._publish_receipt(output, self.receipt)

    def test_no_winner_publishes_nothing(self) -> None:
        output = self.root / "development-selection.json"
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=None), self.assertRaises(subject.SelectionError):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertFalse(output.exists())

    def test_outer_publication_transaction_rolls_back_terminal_failures(
            self) -> None:
        forged = copy.deepcopy(self.receipt)
        forged["recovery_run_summary_sha256"] = "9" * 64
        cases = (
            ("reopen", mock.patch.object(
                subject, "_published_receipt",
                side_effect=subject.SelectionError("injected reopen failure"))),
            ("compare", mock.patch.object(
                subject, "_published_receipt", side_effect=lambda *_args,
                **_kwargs: (
                    forged, subject._receipt_fingerprint(output.stat()),
                    output.read_bytes()))),
            ("seal", mock.patch.object(
                subject.contract_api,
                "_seal_validated_architecture_selection",
                side_effect=contract_api.ContractError(
                    "injected seal failure"))),
        )
        for name, terminal_patch in cases:
            output = self.root / (name + "-selection.json")
            with self.subTest(failure=name), mock.patch.object(
                    subject, "recompute_development_architecture_selection",
                    return_value=self.receipt), terminal_patch, \
                    self.assertRaises((subject.SelectionError,
                                       contract_api.ContractError)):
                subject.select_development_architecture(
                    self.contract, self.root / "recovery",
                    self.root / "work", self.root / "timing", output)
            self.assertFalse(output.exists())

    def test_terminal_identity_recheck_rejects_same_byte_replacement(
            self) -> None:
        output = self.root / "development-selection.json"
        replacement = self.root / "same-byte-replacement.json"
        self.write_receipt(replacement)
        original_seal = \
            contract_api._seal_validated_architecture_selection

        def replace_during_seal(contract, value):
            os.replace(str(replacement), str(output))
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                    subject.contract_api,
                    "_seal_validated_architecture_selection",
                    side_effect=replace_during_seal), \
                self.assertRaises(subject.SelectionError):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        # The rollback must not unlink an attacker-installed different inode,
        # but the transaction must never report success for it.
        self.assertTrue(output.exists())

    def test_terminal_fingerprint_recheck_rejects_in_place_mutation(
            self) -> None:
        output = self.root / "development-selection.json"
        original_seal = \
            contract_api._seal_validated_architecture_selection

        def mutate_during_seal(contract, value):
            with output.open("r+b", buffering=0) as stream:
                stream.seek(0)
                stream.write(b"[")
                os.fsync(stream.fileno())
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                    subject.contract_api,
                    "_seal_validated_architecture_selection",
                    side_effect=mutate_during_seal), \
                self.assertRaises(subject.SelectionError):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        # In-place mutation touched the exact transaction inode, so rollback
        # removes it instead of leaving corrupt published output behind.
        self.assertFalse(output.exists())

    def test_publication_rolls_back_destination_after_parent_fsync_error(
            self) -> None:
        output = self.root / "development-selection.json"
        original_fsync = os.fsync

        def failing_directory_fsync(descriptor: int) -> None:
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise OSError(errno.EIO, "injected directory fsync failure")
            original_fsync(descriptor)

        with mock.patch.object(
                subject.os, "fsync", side_effect=failing_directory_fsync), \
                self.assertRaises(subject.SelectionError):
            subject._publish_receipt(output, self.receipt)
        self.assertFalse(output.exists())

    def test_publication_rolls_back_destination_after_identity_error(
            self) -> None:
        output = self.root / "development-selection.json"
        original_stat = os.stat
        substituted = False

        def mismatching_stat(path, *args, **kwargs):
            nonlocal substituted
            result = original_stat(path, *args, **kwargs)
            if (path == output.name and kwargs.get("dir_fd") is not None and
                    not substituted):
                substituted = True
                return SimpleNamespace(
                    st_dev=result.st_dev, st_ino=result.st_ino + 1)
            return result

        with mock.patch.object(
                subject.os, "stat", side_effect=mismatching_stat), \
                self.assertRaises(subject.SelectionError):
            subject._publish_receipt(output, self.receipt)
        self.assertTrue(substituted)
        self.assertFalse(output.exists())

    def test_publication_rolls_back_if_parent_identity_changes(self) -> None:
        parent = self.root / "publish-parent"
        parent.mkdir()
        moved = self.root / "moved-publish-parent"
        output = parent / "development-selection.json"
        original_verify = subject._verify_real_parent
        checks = 0

        def replace_parent(path, identity):
            nonlocal checks
            checks += 1
            if checks == 2:
                parent.rename(moved)
                parent.mkdir()
            original_verify(path, identity)

        with mock.patch.object(
                subject, "_verify_real_parent", side_effect=replace_parent), \
                self.assertRaises(subject.SelectionError):
            subject._publish_receipt(output, self.receipt)
        self.assertFalse(output.exists())
        self.assertFalse((moved / output.name).exists())

    def test_authoritative_load_recomputes_exact_receipt(self) -> None:
        output = self.root / "development-selection.json"
        self.write_receipt(output)
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt):
            handle = subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertEqual(handle.selection_sha256,
                         self.receipt["selection_sha256"])

        forged = copy.deepcopy(self.receipt)
        forged["recovery_run_summary_sha256"] = "a" * 64
        unsigned = {key: value for key, value in forged.items()
                    if key != "selection_sha256"}
        forged["selection_sha256"] = contract_api.sha256_json(unsigned)
        self.write_receipt(output, forged)
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), \
                self.assertRaises(subject.SelectionError):
            subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)

    def test_authoritative_load_rejects_post_read_same_byte_replacement(
            self) -> None:
        output = self.root / "development-selection.json"
        replacement = self.root / "same-byte-replacement.json"
        self.write_receipt(output)
        self.write_receipt(replacement)
        original_seal = \
            contract_api._seal_validated_architecture_selection

        def replace_during_seal(contract, value):
            os.replace(str(replacement), str(output))
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                    subject.contract_api,
                    "_seal_validated_architecture_selection",
                    side_effect=replace_during_seal), \
                self.assertRaises(subject.SelectionError):
            subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertTrue(output.exists())

    def test_authoritative_load_rejects_post_read_in_place_mutation(
            self) -> None:
        output = self.root / "development-selection.json"
        self.write_receipt(output)
        original_seal = \
            contract_api._seal_validated_architecture_selection

        def mutate_during_seal(contract, value):
            with output.open("r+b", buffering=0) as stream:
                stream.seek(0)
                stream.write(b"[")
                os.fsync(stream.fileno())
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                    subject.contract_api,
                    "_seal_validated_architecture_selection",
                    side_effect=mutate_during_seal), \
                self.assertRaises(subject.SelectionError):
            subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertTrue(output.exists())

    def test_recovery_loader_rejects_independent_work_join_drift(self) -> None:
        work_summary = "a" * 64
        work_result = "b" * 64
        work_domain = "c" * 64
        witness = {"strict": "witness"}
        witness_hash = contract_api.sha256_json(witness)
        campaign = {
            "freeze": {
                "source_git_commit": "1" * 40,
                "timing_proxy_witness_sha256": witness_hash,
                "work_rank_summary_sha256": work_summary,
                "work_rank_result_stream_sha256": work_result,
                "work_rank_domain_sha256": work_domain,
            },
            "rows": [],
            "summary": {
                "summary_sha256": "d" * 64,
                "timing_proxy_witness_sha256": witness_hash,
                "work_rank_summary_sha256": work_summary,
                "work_rank_result_stream_sha256": work_result,
                "work_rank_domain_sha256": work_domain,
            },
            "receipt": {
                "result_stream_sha256": "e" * 64,
                "receipt_sha256": "f" * 64,
            },
            "timing_proxy_witness": witness,
        }
        binding = {
            "work_rank_summary_sha256": work_summary,
            "work_rank_result_stream_sha256": work_result,
            "work_rank_domain_sha256": work_domain,
        }
        patches = (
            mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=campaign),
            mock.patch.object(
                subject.work_api, "load_completed_work_screen",
                return_value={}),
            mock.patch.object(
                subject, "_revalidate_recovery_summary", return_value={}),
        )
        with patches[0], patches[1], patches[2], mock.patch.object(
                subject.recovery_api, "_bind_work_rank_identities",
                return_value={**binding,
                              "work_rank_domain_sha256": "9" * 64}), \
                self.assertRaises(subject.SelectionError):
            subject._load_recovery_and_work(
                self.contract, self.root / "recovery", self.root / "work")
        with patches[0], patches[1], patches[2], mock.patch.object(
                subject.recovery_api, "_bind_work_rank_identities",
                side_effect=subject.recovery_api.RecoveryRunnerError(
                    "injected identity join mismatch")), \
                self.assertRaises(subject.SelectionError):
            subject._load_recovery_and_work(
                self.contract, self.root / "recovery", self.root / "work")

    def timing_loader_bundle(self):
        arms = [
            {"arm": "wirehair2_head", "codec": "wirehair2_certified",
             "binary_sha256": "a" * 64,
             "arm_descriptor_sha256": "b" * 64},
            {"arm": "wirehair1", "codec": "wirehair1",
             "binary_sha256": "a" * 64,
             "arm_descriptor_sha256": "c" * 64},
            {"arm": "wirehair2_dense_two07_basis_v1",
             "codec": "wirehair2_experiment", "binary_sha256": "a" * 64,
             "arm_descriptor_sha256":
                "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388"},
        ]
        freeze = {
            "schema": contract_api.FREEZE_SCHEMA,
            "source_git_commit": "1" * 40,
            "arm_roster": [arm["arm"] for arm in arms],
            "arms": arms,
        }
        summary = {"validated": True}
        receipt = {
            "result_stream_sha256": "d" * 64,
            "receipt_sha256": "e" * 64,
            "timing_qualification_execution_receipt_sha256": "f" * 64,
        }
        run_summary = {
            "summary_sha256": "1" * 64,
            "timing_freeze_sha256":
                contract_api.freeze_manifest_sha256(freeze),
            "timing_architecture_artifact_sha256":
                contract_api.architecture_artifact_sha256(freeze),
            "timing_result_sha256": receipt["result_stream_sha256"],
            "timing_execution_receipt_sha256": receipt["receipt_sha256"],
            "timing_validator_summary_sha256":
                contract_api.sha256_json(summary),
            "timing_qualification_execution_receipt_sha256":
                receipt["timing_qualification_execution_receipt_sha256"],
        }
        return {
            "directory": str(self.root), "directory_identity": (1, 2),
            "run_summary": run_summary, "freeze": freeze,
            "summary": summary, "execution_receipt": receipt,
            "timing_qualification": self.qualification,
        }

    def test_timing_loader_rejects_interface_and_terminal_hash_drift(
            self) -> None:
        valid = self.timing_loader_bundle()
        loader_name = "load_completed_timing_screen"
        with mock.patch.object(
                subject.timing_api, loader_name, return_value=valid,
                create=True):
            loaded = subject._load_timing(
                self.contract, self.root / "timing")
        self.assertEqual(
            loaded["provenance"]["timing_result_stream_sha256"], "d" * 64)
        mutations = []
        changed = dict(valid)
        changed["extra"] = True
        mutations.append(("interface", changed))
        for field in (
                "timing_result_sha256", "timing_execution_receipt_sha256",
                "timing_qualification_execution_receipt_sha256",
                "timing_validator_summary_sha256", "timing_freeze_sha256",
                "timing_architecture_artifact_sha256"):
            changed = dict(valid)
            changed["run_summary"] = copy.deepcopy(valid["run_summary"])
            changed["run_summary"][field] = "9" * 64
            mutations.append((field, changed))
        for name, changed in mutations:
            with self.subTest(mutation=name), mock.patch.object(
                    subject.timing_api, loader_name, return_value=changed,
                    create=True), self.assertRaises(subject.SelectionError):
                subject._load_timing(
                    self.contract, self.root / "timing")

    def test_receipt_reader_rejects_symlink_fifo_oversize_and_directory(
            self) -> None:
        real = self.root / "real.json"
        self.write_receipt(real)
        symlink = self.root / "symlink.json"
        symlink.symlink_to(real)
        fifo = self.root / "receipt.fifo"
        os.mkfifo(str(fifo))
        oversize = self.root / "oversize.json"
        oversize.write_bytes(
            b"x" * (subject.MAX_SELECTION_RECEIPT_BYTES + 1))
        for name, path in (
                ("symlink", symlink), ("fifo", fifo),
                ("oversize", oversize), ("directory", self.root)):
            with self.subTest(kind=name), \
                    self.assertRaises(subject.SelectionError):
                subject._published_receipt(self.contract, path)
        with self.assertRaises(subject.SelectionError):
            subject._publish_receipt(Path("/"), self.receipt)

    def test_receipt_reader_rejects_late_directory_entry_replacement(
            self) -> None:
        path = self.root / "development-selection.json"
        replacement = self.root / "replacement.json"
        self.write_receipt(path)
        self.write_receipt(replacement)
        original_read = os.read
        replaced = False

        def racing_read(descriptor: int, count: int) -> bytes:
            nonlocal replaced
            data = original_read(descriptor, count)
            if data and not replaced:
                replaced = True
                os.replace(str(replacement), str(path))
            return data

        with mock.patch.object(subject.os, "read", side_effect=racing_read), \
                self.assertRaises(subject.SelectionError):
            subject._published_receipt(self.contract, path)
        self.assertTrue(replaced)

    def test_receipt_reader_rejects_parent_identity_replacement(self) -> None:
        parent = self.root / "receipt-parent"
        parent.mkdir()
        moved = self.root / "moved-receipt-parent"
        path = parent / "development-selection.json"
        self.write_receipt(path)
        original_read = os.read
        replaced = False

        def racing_read(descriptor: int, count: int) -> bytes:
            nonlocal replaced
            data = original_read(descriptor, count)
            if data and not replaced:
                replaced = True
                parent.rename(moved)
                parent.mkdir()
            return data

        with mock.patch.object(subject.os, "read", side_effect=racing_read), \
                self.assertRaises(subject.SelectionError):
            subject._published_receipt(self.contract, path)
        self.assertTrue(replaced)

    def test_summary_only_cli_is_not_exposed(self) -> None:
        with self.assertRaises(SystemExit):
            subject._parser().parse_args([
                "--recovery-summary", "recovery.json",
                "--timing-summary", "timing.json",
                "--output", str(self.root / "selection.json"),
            ])


if __name__ == "__main__":
    unittest.main()
