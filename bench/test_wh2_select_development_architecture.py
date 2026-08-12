#!/usr/bin/env python3
"""Adversarial tests for the path-only WH2 architecture selector."""

from __future__ import annotations

import copy
import errno
import os
from pathlib import Path
import stat
import tempfile
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

    def publish_receipt(self, path: Path, value=None):
        if value is None:
            value = self.receipt
        transaction = subject._ReceiptPublicationTransaction(path)
        try:
            transaction.open_parent()
            selection = subject._prepare_receipt(
                transaction, self.contract, value)
            subject._publish_receipt(transaction)
            return selection
        finally:
            transaction.close()

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
        original_bytes = output.read_bytes()
        with self.assertRaises(subject.SelectionError):
            self.publish_receipt(output)
        self.assertEqual(output.read_bytes(), original_bytes)

    def test_no_winner_publishes_nothing(self) -> None:
        output = self.root / "development-selection.json"
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=None), self.assertRaises(subject.SelectionError):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertFalse(output.exists())

    def test_validator_and_seal_failures_leave_no_public_name(self) -> None:
        cases = (
            ("validator", "validate_selection_receipt"),
            ("seal", "_seal_validated_architecture_selection"),
        )
        for name, target in cases:
            output = self.root / (name + "-selection.json")
            with self.subTest(failure=name), mock.patch.object(
                    subject, "recompute_development_architecture_selection",
                    return_value=self.receipt), mock.patch.object(
                    subject.contract_api, target,
                    side_effect=contract_api.ContractError(
                        "injected {} failure".format(name))), \
                    self.assertRaises(subject.SelectionError):
                subject.select_development_architecture(
                    self.contract, self.root / "recovery",
                    self.root / "work", self.root / "timing", output)
            self.assertFalse(output.exists())

    def test_malformed_seal_handle_is_rejected_before_publication(self) -> None:
        output = self.root / "malformed-seal-selection.json"
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.contract_api,
                "_seal_validated_architecture_selection",
                return_value=None), self.assertRaisesRegex(
                    subject.SelectionError, "malformed handle"):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertFalse(output.exists())

    def test_fd_ownership_cleanup_closes_each_descriptor_once(self) -> None:
        class InjectedAsyncFailure(BaseException):
            pass

        parent_output = self.root / "parent-ownership.json"
        transaction = subject._ReceiptPublicationTransaction(parent_output)
        original_open = os.open
        original_close = os.close
        parent_fds = []
        closed_fds = []

        def record_parent_open(*args, **kwargs):
            descriptor = original_open(*args, **kwargs)
            parent_fds.append(descriptor)
            return descriptor

        def interrupt_first_fstat(_descriptor):
            raise InjectedAsyncFailure("after parent descriptor transfer")

        def record_close(descriptor):
            closed_fds.append(descriptor)
            return original_close(descriptor)

        with mock.patch.object(
                subject.os, "open", side_effect=record_parent_open), \
                mock.patch.object(
                    subject.os, "fstat", side_effect=interrupt_first_fstat), \
                mock.patch.object(subject.os, "close", side_effect=record_close), \
                self.assertRaises(InjectedAsyncFailure):
            transaction.open_parent()
        transaction.close()
        self.assertEqual(len(parent_fds), 1)
        self.assertEqual(closed_fds.count(parent_fds[0]), 1)

        output = self.root / "prepared-ownership.json"
        original_unnamed = subject._open_unnamed_receipt
        prepared_fds = []
        closed_fds = []

        def record_unnamed(current_transaction):
            descriptor = original_unnamed(current_transaction)
            prepared_fds.append(descriptor)
            return descriptor

        with mock.patch.object(
                subject, "_open_unnamed_receipt",
                side_effect=record_unnamed), mock.patch.object(
                subject.contract_api, "validate_selection_receipt",
                side_effect=contract_api.ContractError(
                    "after prepared descriptor transfer")), \
                mock.patch.object(subject.os, "close", side_effect=record_close), \
                self.assertRaises(subject.SelectionError):
            self.publish_receipt(output)
        self.assertEqual(len(prepared_fds), 1)
        self.assertEqual(closed_fds.count(prepared_fds[0]), 1)
        self.assertFalse(output.exists())

    def test_unnamed_stage_failure_leaves_no_directory_residue(self) -> None:
        output = self.root / "development-selection.json"
        before = tuple(self.root.iterdir())
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.contract_api, "validate_selection_receipt",
                side_effect=contract_api.ContractError(
                    "injected validation failure")), \
                self.assertRaises(subject.SelectionError):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertEqual(tuple(self.root.iterdir()), before)

    def test_receipt_is_not_visible_during_validation_or_sealing(self) -> None:
        output = self.root / "development-selection.json"
        original_validate = contract_api.validate_selection_receipt
        original_seal = \
            contract_api._seal_validated_architecture_selection
        observed = []

        def validate_while_unnamed(contract, value):
            observed.append(("validate", output.exists(), list(self.root.iterdir())))
            return original_validate(contract, value)

        def seal_while_unnamed(contract, value):
            observed.append(("seal", output.exists(), list(self.root.iterdir())))
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.contract_api, "validate_selection_receipt",
                side_effect=validate_while_unnamed), mock.patch.object(
                subject.contract_api,
                "_seal_validated_architecture_selection",
                side_effect=seal_while_unnamed):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertEqual(
            observed,
            [("validate", False, []), ("seal", False, []),
             ("validate", False, [])])
        self.assertTrue(output.exists())

    def test_link_return_interrupt_keeps_fully_valid_committed_receipt(
            self) -> None:
        output = self.root / "development-selection.json"
        original_link = os.link

        class InjectedAsyncFailure(BaseException):
            pass

        def link_then_interrupt(*args, **kwargs):
            original_link(*args, **kwargs)
            self.assertTrue(output.exists())
            raise InjectedAsyncFailure("lost publication acknowledgement")

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.os, "link", side_effect=link_then_interrupt), \
                self.assertRaises(InjectedAsyncFailure):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        published = subject._published_receipt(self.contract, output)[0]
        self.assertEqual(published, self.receipt)
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt):
            loaded = subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertIs(type(loaded), contract_api.ArchitectureSelection)
        with self.assertRaises(subject.SelectionError):
            self.publish_receipt(output)

    def test_link_return_oserror_still_keeps_valid_committed_receipt(
            self) -> None:
        output = self.root / "development-selection-oserror.json"
        original_link = os.link

        def link_then_error(*args, **kwargs):
            original_link(*args, **kwargs)
            raise OSError(errno.EIO, "lost link acknowledgement")

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.os, "link", side_effect=link_then_error), \
                self.assertRaisesRegex(
                    subject.SelectionError, "lost link acknowledgement"):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        published = subject._published_receipt(self.contract, output)[0]
        self.assertEqual(published, self.receipt)
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt):
            loaded = subject.load_authoritative_selection(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertIs(type(loaded), contract_api.ArchitectureSelection)

    def test_retained_unnamed_fd_has_exact_identity_and_content(self) -> None:
        output = self.root / "development-selection.json"
        expected_bytes = (
            contract_api.canonical_json(self.receipt) + "\n").encode("utf-8")
        original_link = os.link
        observed = []
        transaction = subject._ReceiptPublicationTransaction(output)
        try:
            transaction.open_parent()
            selection = subject._prepare_receipt(
                transaction, self.contract, self.receipt)
            prepared = transaction.require_prepared()

            def inspect_retained_fd(source, destination, **kwargs):
                descriptor = int(Path(source).name)
                opened = os.fstat(descriptor)
                retained_bytes, fingerprint = subject._read_unnamed_receipt(
                    descriptor, expected_identity=prepared.identity)
                self.assertTrue(stat.S_ISREG(opened.st_mode))
                self.assertEqual(opened.st_nlink, 0)
                self.assertEqual(descriptor, prepared.descriptor)
                self.assertEqual(destination, output.name)
                self.assertEqual(kwargs["dst_dir_fd"], transaction.parent_fd)
                self.assertTrue(kwargs["follow_symlinks"])
                self.assertEqual(retained_bytes, expected_bytes)
                self.assertEqual(fingerprint[0:2], prepared.identity)
                observed.append(prepared.identity)
                return original_link(source, destination, **kwargs)

            with mock.patch.object(
                    subject.os, "link", side_effect=inspect_retained_fd):
                subject._publish_receipt(transaction)
            self.assertIs(
                type(selection), contract_api.ArchitectureSelection)
        finally:
            transaction.close()
        self.assertEqual(observed, [(output.stat().st_dev, output.stat().st_ino)])
        self.assertEqual(output.read_bytes(), expected_bytes)

    def test_prelink_content_mutation_is_rejected_without_public_name(
            self) -> None:
        output = self.root / "mutated-selection.json"
        transaction = subject._ReceiptPublicationTransaction(output)
        try:
            transaction.open_parent()
            subject._prepare_receipt(
                transaction, self.contract, self.receipt)
            prepared = transaction.require_prepared()
            os.pwrite(prepared.descriptor, b"[", 0)
            os.fsync(prepared.descriptor)
            with self.assertRaisesRegex(
                    subject.SelectionError, "changed before final link"):
                subject._publish_receipt(transaction)
        finally:
            transaction.close()
        self.assertFalse(output.exists())

    def test_parent_path_replacement_immediately_prelink_is_rejected(
            self) -> None:
        parent = self.root / "publish-parent"
        parent.mkdir()
        moved = self.root / "moved-publish-parent"
        output = parent / "development-selection.json"
        original_seal = \
            contract_api._seal_validated_architecture_selection

        def replace_parent_before_link(contract, value):
            parent.rename(moved)
            parent.mkdir()
            return original_seal(contract, value)

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.contract_api,
                "_seal_validated_architecture_selection",
                side_effect=replace_parent_before_link), \
                self.assertRaisesRegex(
                    subject.SelectionError, "parent identity changed"):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertFalse(output.exists())
        self.assertFalse((moved / output.name).exists())

    def test_unsupported_tmpfile_and_proc_link_fail_without_a_name(self) -> None:
        output = self.root / "development-selection.json"
        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject, "_open_unnamed_receipt",
                side_effect=subject.SelectionError(
                    "injected O_TMPFILE failure")), mock.patch.object(
                subject.os, "link") as link, self.assertRaisesRegex(
                    subject.SelectionError, "O_TMPFILE failure"):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        link.assert_not_called()
        self.assertFalse(output.exists())

        with mock.patch.object(
                subject, "recompute_development_architecture_selection",
                return_value=self.receipt), mock.patch.object(
                subject.os, "link",
                side_effect=OSError(errno.EOPNOTSUPP,
                                    "injected proc-link failure")), \
                self.assertRaisesRegex(
                    subject.SelectionError, "proc-link failure"):
            subject.select_development_architecture(
                self.contract, self.root / "recovery", self.root / "work",
                self.root / "timing", output)
        self.assertFalse(output.exists())

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
            self.publish_receipt(Path("/"))

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
