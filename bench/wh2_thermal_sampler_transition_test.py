#!/usr/bin/env python3
"""Deterministic tests for the bounded thermal transition controller."""

from __future__ import annotations

from dataclasses import replace
from contextlib import redirect_stderr
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import signal
import stat
import sys
import tempfile
import unittest
from unittest import mock

import wh2_thermal_sampler_transition as transition


class FakeClock:
    def __init__(self) -> None:
        self.value = 100.0

    def __call__(self) -> float:
        return self.value


class FakeJournal:
    def __init__(self, fail_at=None, *, clock=None, advance_at=None) -> None:
        self.fail_at = fail_at
        self.clock = clock
        self.advance_at = advance_at
        self.records = []
        self.replay_calls = 0

    def record(self, phase, status, payload):
        self.records.append((phase, status, dict(payload)))
        if self.advance_at == (phase, status) and self.clock is not None:
            self.clock.value += 541.0
        if self.fail_at == (phase, status):
            raise transition.TransitionError("injected journal failure")
        return {"phase": phase, "status": status}

    def replay(self):
        self.replay_calls += 1
        if self.advance_at == (
                "receipt_chain", "replay-%d" % self.replay_calls) and \
                self.clock is not None:
            self.clock.value += 541.0
        if self.fail_at == ("receipt_chain", "replay"):
            raise transition.TransitionError("injected journal failure")
        return {
            "count": len(self.records),
            "head_sha256": "c" * 64,
            "roster": [
                {"sequence": index, "phase": phase, "status": status}
                for index, (phase, status, _payload) in enumerate(self.records)
            ],
        }


class FakeSignalGuard:
    def __init__(self) -> None:
        self.requested = False

    def raise_if_requested(self):
        if self.requested:
            raise transition.TransitionError("deferred controller signal: SIGTERM")


class FakeBackend:
    NORMAL = (
        "hard_preflight", "arm_recovery", "stop_old", "archive_old",
        "start_candidate", "exercise_candidate", "stop_candidate",
        "accept_candidate",
    )

    def __init__(self, fail_at=None, clock=None, advance_at=None,
                 *, forced_stop=False, forced_archive=False) -> None:
        self.fail_at = fail_at
        self.clock = clock
        self.advance_at = advance_at
        self.calls = []
        self.recovery_budget = None
        self.forced_stop = forced_stop
        self.forced_archive = forced_archive

    def _call(self, name, value=None):
        self.calls.append(name)
        if self.advance_at == name and self.clock is not None:
            self.clock.value += 500.0
        if self.fail_at == name:
            raise transition.TransitionError("injected " + name)
        return value if value is not None else {"phase": name}

    def hard_preflight(self):
        return self._call("hard_preflight")

    def arm_recovery(self):
        return self._call("arm_recovery")

    def stop_old(self):
        return self._call("stop_old", {"forced": self.forced_stop})

    def archive_old(self):
        return self._call("archive_old", {
            "binding": {"sha256": "a" * 64},
            "path": "/archive", "forced_stop": self.forced_archive,
        })

    def start_candidate(self):
        return self._call("start_candidate")

    def exercise_candidate(self):
        return self._call("exercise_candidate")

    def stop_candidate(self):
        return self._call("stop_candidate")

    def accept_candidate(self):
        return self._call("accept_candidate", {
            "sample_count": 6, "raw_sha256": "d" * 64})

    def cleanup_candidate(self):
        return self._call("cleanup_candidate")

    def begin_emergency_recovery(self, budget):
        self.recovery_budget = budget
        return {"budget_installed": True}

    def restore_old(self, archive):
        if archive is not None:
            assert archive["binding"]["sha256"] == "a" * 64
        return self._call("restore_old", {
            "pid": 7001,
            "csv_live_identity": {"device": 1, "inode": 2},
            "archived_pre_dry_sha256": "a" * 64,
        })

    def publish_audit_binding(self, archive, restored, dry_accepted,
                              candidate_accept, receipt_chain_prefix):
        assert restored["pid"] == 7001
        assert receipt_chain_prefix["count"] > 0
        if dry_accepted:
            assert candidate_accept["sample_count"] == 6
        return self._call("publish_audit_binding", {
            "dry_accepted": dry_accepted, "sha256": "b" * 64,
        })

    def final_replay(self, candidate_accept, archive, restored, audit):
        assert restored["pid"] == 7001
        assert audit["sha256"] == "b" * 64
        return self._call("final_replay", {
            "candidate_accept": candidate_accept,
            "archive": archive,
            "receipt": "replayed",
        })


class TransitionStateMachineTests(unittest.TestCase):
    def make_controller(self, backend, journal=None, clock=None,
                        signal_guard=None):
        clock = clock or FakeClock()
        deadline = transition.Deadline(540.0, 60.0, clock=clock)
        return transition.TransitionController(
            transition.TransitionPlan(controller_uid=os.geteuid()), backend,
            journal or FakeJournal(), deadline, signal_guard), clock

    def test_success_is_nonpromotional_and_restores_old(self):
        backend = FakeBackend()
        journal = FakeJournal()
        controller, _clock = self.make_controller(backend, journal)
        terminal = controller.run()
        self.assertTrue(terminal["success"])
        self.assertTrue(terminal["dry_accepted"])
        self.assertEqual(terminal["future_audit_binding"]["sha256"], "b" * 64)
        self.assertEqual(backend.calls, [
            *FakeBackend.NORMAL, "cleanup_candidate", "restore_old",
            "publish_audit_binding", "final_replay",
        ])
        self.assertIn(("recovery_arm", "completed", {"phase": "arm_recovery"}),
                      journal.records)
        self.assertEqual(terminal["candidate_accept"]["sample_count"], 6)
        self.assertEqual(terminal["final_replay"]["receipt"], "replayed")
        self.assertGreater(terminal["receipt_chain_prefix"]["count"], 0)
        checkpoints = [payload for phase, status, payload in journal.records
                       if (phase, status) == ("terminal", "started")]
        self.assertEqual(len(checkpoints), 1)
        self.assertIsNone(checkpoints[0]["transition_success"])
        self.assertNotIn("success", checkpoints[0])

    def test_forced_old_stop_vetoes_candidate_and_enters_recovery_immediately(self):
        backend = FakeBackend(forced_stop=True)
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "forced old stop vetoes"):
            controller.run()
        self.assertNotIn("archive_old", backend.calls)
        self.assertNotIn("start_candidate", backend.calls)
        self.assertEqual(backend.calls[:3], [
            "hard_preflight", "arm_recovery", "stop_old"])
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_unclean_archive_vetoes_candidate_before_launch(self):
        backend = FakeBackend(forced_archive=True)
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "forced old stop vetoes"):
            controller.run()
        self.assertIn("archive_old", backend.calls)
        self.assertNotIn("start_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_preflight_failure_never_enters_recovery(self):
        backend = FakeBackend(fail_at="hard_preflight")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "hard_preflight"):
            controller.run()
        self.assertEqual(backend.calls, ["hard_preflight"])

    def test_recovery_arm_failure_never_stops_old(self):
        backend = FakeBackend(fail_at="arm_recovery")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "arm_recovery"):
            controller.run()
        self.assertEqual(backend.calls, ["hard_preflight", "arm_recovery"])

    def test_every_post_arm_boundary_failure_enters_recovery(self):
        for failure in FakeBackend.NORMAL[2:]:
            with self.subTest(failure=failure):
                backend = FakeBackend(fail_at=failure)
                controller, _clock = self.make_controller(backend)
                with self.assertRaisesRegex(transition.TransitionError, failure):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)
                self.assertIn("publish_audit_binding", backend.calls)
                self.assertLess(
                    backend.calls.index("cleanup_candidate"),
                    backend.calls.index("restore_old"))

    def test_complete_receipt_failure_after_stop_still_recovers(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_stop", "completed"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_late_accept_completion_retains_result_but_is_not_accepted(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("candidate_accept", "completed"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertEqual(controller.candidate_accept["sample_count"], 6)
        self.assertFalse(controller.dry_accepted)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_stop_start_receipt_failure_leaves_old_untouched(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_stop", "started"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertNotIn("stop_old", backend.calls)
        self.assertNotIn("cleanup_candidate", backend.calls)
        self.assertNotIn("restore_old", backend.calls)

    def test_stop_start_receipt_overrun_leaves_old_untouched(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("old_stop", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertNotIn("stop_old", backend.calls)
        self.assertNotIn("cleanup_candidate", backend.calls)
        self.assertNotIn("restore_old", backend.calls)

    def test_post_stop_phase_start_overrun_recovers_without_running_action(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("candidate_start", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertNotIn("start_candidate", backend.calls)
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_started_receipt_failure_after_arm_still_recovers(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_archive", "started"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_every_post_stop_phase_receipt_failure_recovers(self):
        boundaries = [("old_stop", "completed")]
        for phase in (
                "old_archive", "candidate_start", "candidate_exercise",
                "candidate_stop", "candidate_accept"):
            boundaries.extend(((phase, "started"), (phase, "completed")))
        for boundary in boundaries:
            with self.subTest(boundary=boundary):
                backend = FakeBackend()
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=boundary))
                with self.assertRaisesRegex(transition.TransitionError,
                                            "journal failure"):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_every_post_stop_failed_receipt_failure_still_recovers(self):
        failures = {
            "stop_old": "old_stop",
            "archive_old": "old_archive",
            "start_candidate": "candidate_start",
            "exercise_candidate": "candidate_exercise",
            "stop_candidate": "candidate_stop",
            "accept_candidate": "candidate_accept",
        }
        for action, phase in failures.items():
            with self.subTest(action=action):
                backend = FakeBackend(fail_at=action)
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")))
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_every_recovery_receipt_failure_preserves_action_order(self):
        for boundary in (
                ("emergency_recovery", "started"),
                ("candidate_cleanup", "started"),
                ("candidate_cleanup", "completed"),
                ("old_restore", "started"),
                ("old_restore", "completed"),
                ("audit_binding", "started"),
                ("audit_binding", "completed"),
                ("final_replay", "started"),
                ("final_replay", "completed"),
                ("emergency_recovery", "completed"),
                ("terminal", "started")):
            with self.subTest(boundary=boundary):
                backend = FakeBackend()
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=boundary))
                with self.assertRaisesRegex(transition.TransitionError,
                                            "journal failure"):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)
                self.assertLess(
                    backend.calls.index("cleanup_candidate"),
                    backend.calls.index("restore_old"))

    def test_every_recovery_failed_receipt_failure_preserves_safe_progress(self):
        cases = (
            ("cleanup_candidate", "candidate_cleanup", True, True),
            ("restore_old", "old_restore", True, False),
            ("publish_audit_binding", "audit_binding", True, True),
            ("final_replay", "final_replay", True, True),
        )
        for action, phase, restore_expected, audit_expected in cases:
            with self.subTest(action=action):
                backend = FakeBackend(fail_at=action)
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")))
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertEqual("restore_old" in backend.calls,
                                 restore_expected)
                self.assertEqual("publish_audit_binding" in backend.calls,
                                 audit_expected)

    def test_cleanup_failure_does_not_skip_restore_attempt(self):
        backend = FakeBackend(fail_at="cleanup_candidate")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "cleanup_candidate"):
            controller.run()
        self.assertIn("restore_old", backend.calls)

    def test_restore_failure_is_terminal_and_skips_audit_publication(self):
        backend = FakeBackend(fail_at="restore_old")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "restore_old"):
            controller.run()
        self.assertNotIn("publish_audit_binding", backend.calls)

    def test_failed_terminal_receipt_failure_follows_completed_recovery(self):
        backend = FakeBackend(fail_at="accept_candidate")
        controller, _clock = self.make_controller(
            backend, FakeJournal(fail_at=("terminal", "started")))
        with self.assertRaises(transition.TransitionError):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_audit_binding_failure_is_terminal_after_restore(self):
        backend = FakeBackend(fail_at="publish_audit_binding")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "publish_audit_binding"):
            controller.run()
        self.assertIn("restore_old", backend.calls)

    def test_signal_received_during_recovery_is_not_lost(self):
        guard = FakeSignalGuard()
        backend = FakeBackend()
        original = backend.publish_audit_binding

        def signal_during_audit(archive, restored, dry_accepted,
                                candidate_accept, receipt_chain_prefix):
            result = original(
                archive, restored, dry_accepted, candidate_accept,
                receipt_chain_prefix)
            guard.requested = True
            return result

        backend.publish_audit_binding = signal_during_audit
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal, signal_guard=guard)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "deferred controller signal"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        terminal = [record for record in journal.records
                    if record[0] == "terminal"]
        self.assertEqual(terminal[-1][1], "started")
        self.assertIsNone(terminal[-1][2]["transition_success"])

    def test_normal_deadline_failure_enters_recovery_reserve(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock, advance_at="stop_old")
        controller, _clock = self.make_controller(backend, clock=clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "normal deadline"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_expired_scoring_deadline_fails_but_never_vetoes_recovery(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock, advance_at="stop_old")
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        clock.value = 200.0
        backend.advance_at = "stop_old"
        # Advance beyond the full 540 s at the first post-arm operation.
        backend._call = self._absolute_deadline_call(backend, clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "scoring absolute deadline exhausted"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)
        self.assertIsNotNone(backend.recovery_budget)
        self.assertEqual(backend.recovery_budget.started, clock.value)
        self.assertEqual(
            backend.recovery_budget.absolute,
            clock.value + controller.plan.emergency_recovery_s)
        deadline_receipts = [record for record in journal.records
                             if record[0] == "scoring_deadline"]
        self.assertEqual(len(deadline_receipts), 1)
        self.assertEqual(deadline_receipts[0][1], "failed")
        self.assertTrue(deadline_receipts[0][2]["recovery_actions_continue"])

    def test_terminal_checkpoint_write_crossing_absolute_forces_failure(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("terminal", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "scoring absolute deadline exhausted"):
            controller.run()
        self.assertGreaterEqual(clock.value, controller.deadline.absolute)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("final_replay", backend.calls)
        checkpoints = [payload for phase, _status, payload in journal.records
                       if phase == "terminal"]
        self.assertTrue(checkpoints)
        self.assertTrue(all("success" not in payload for payload in checkpoints))
        self.assertTrue(all(payload["transition_success"] is None
                            for payload in checkpoints))
        self.assertFalse(any(phase == "terminal_verdict"
                             for phase, _status, _payload in journal.records))

    def test_terminal_receipt_replay_crossing_absolute_forces_failure(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("receipt_chain", "replay-2"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "scoring absolute deadline exhausted"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("final_replay", backend.calls)
        self.assertGreaterEqual(
            controller.scoring_evidence_completed_monotonic_s,
            controller.deadline.absolute)

    def test_terminal_uses_one_post_replay_scoring_observation(self):
        class ContradictionClock:
            def __init__(self):
                self.armed = False
                self.calls_after_arm = 0

            def __call__(self):
                if not self.armed:
                    return 100.0
                self.calls_after_arm += 1
                return 100.0 if self.calls_after_arm == 1 else 641.0

        class ArmAfterFinalReplayJournal(FakeJournal):
            def replay(inner_self):
                result = super().replay()
                if inner_self.replay_calls == 2:
                    clock.armed = True
                return result

        clock = ContradictionClock()
        journal = ArmAfterFinalReplayJournal()
        controller, _clock = self.make_controller(
            FakeBackend(), journal=journal, clock=clock)
        terminal = controller.run()
        self.assertTrue(terminal["success"])
        self.assertFalse(terminal["scoring_deadline"]["exhausted"])
        self.assertEqual(clock.calls_after_arm, 1)

    def test_emergency_budget_exhaustion_is_terminal_after_safe_attempts(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock)
        original = backend.restore_old

        def exhaust_during_restore(archive):
            result = original(archive)
            clock.value += 61.0
            return result

        backend.restore_old = exhaust_during_restore
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "emergency recovery budget exhausted"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)
        emergency_receipts = [record for record in journal.records
                              if record[0] == "emergency_recovery"]
        self.assertEqual([record[1] for record in emergency_receipts],
                         ["started", "failed"])
        self.assertTrue(emergency_receipts[-1][2]["exhausted"])

    def test_emergency_setup_failure_does_not_veto_cleanup_or_restore(self):
        backend = FakeBackend()

        def fail_after_install(budget):
            backend.recovery_budget = budget
            raise transition.TransitionError("injected recovery tool replay")

        backend.begin_emergency_recovery = fail_after_install
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "recovery tool replay"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_emergency_and_scoring_failed_receipt_errors_do_not_veto_restore(self):
        cases = ("emergency_recovery", "scoring_deadline")
        for phase in cases:
            with self.subTest(phase=phase):
                clock = FakeClock()
                backend = FakeBackend()
                if phase == "emergency_recovery":
                    def fail_after_install(budget):
                        backend.recovery_budget = budget
                        raise transition.TransitionError(
                            "injected emergency setup")
                    backend.begin_emergency_recovery = fail_after_install
                else:
                    def expire_after_stop(name, value=None):
                        result = FakeBackend._call(backend, name, value)
                        if name == "stop_old":
                            clock.value += 600.0
                        return result
                    backend._call = expire_after_stop
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")), clock)
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_live_backend_installs_emergency_budget_before_tool_replay(self):
        backend = object.__new__(transition.LiveBackend)
        backend.recovery_budget = None
        backend._verify_tools = mock.Mock(side_effect=
            transition.TransitionError("injected tool mismatch"))
        budget = transition.EmergencyRecoveryBudget(60.0)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "tool mismatch"):
            transition.LiveBackend.begin_emergency_recovery(backend, budget)
        self.assertIs(backend.recovery_budget, budget)

    def test_recovery_waits_never_use_the_expired_scoring_deadline(self):
        for method in (
                transition.LiveBackend._quiesce_old_for_recovery,
                transition.LiveBackend._launch_old):
            source = inspect.getsource(method)
            self.assertNotIn("deadline.absolute", source)
            self.assertIn("_recovery_wait_deadline", source)
        clock = FakeClock()
        scoring = transition.Deadline(10.0, 2.0, clock=clock)
        clock.value = scoring.absolute + 25.0
        emergency = scoring.start_emergency_recovery(60.0)
        backend = object.__new__(transition.LiveBackend)
        backend.recovery_budget = emergency
        self.assertGreater(
            transition.LiveBackend._recovery_wait_deadline(backend, 5.0),
            scoring.absolute)
        clock.value = emergency.absolute + 1.0
        extended = transition.LiveBackend._recovery_wait_deadline(
            backend, 20.0, minimum_safety_wait_s=5.0)
        self.assertEqual(extended, clock.value + 5.0)
        self.assertEqual(emergency.receipt()["safety_extension_count"], 1)
        self.assertTrue(emergency.receipt()["exhausted"])

    def test_emergency_receipt_uses_one_consistent_clock_observation(self):
        values = iter((100.0, 159.0, 161.0))
        calls = []

        def clock():
            value = next(values)
            calls.append(value)
            return value

        budget = transition.EmergencyRecoveryBudget(60.0, clock=clock)
        receipt = budget.receipt()
        self.assertFalse(receipt["exhausted"])
        self.assertEqual(receipt["observed_monotonic_s"], 159.0)
        self.assertEqual(receipt["remaining_s"], 1.0)
        self.assertEqual(calls, [100.0, 159.0])

    @staticmethod
    def _absolute_deadline_call(backend, clock):
        original = FakeBackend._call.__get__(backend, FakeBackend)

        def call(name, value=None):
            result = original(name, value)
            if name == "stop_old":
                clock.value += 600.0
            return result

        return call


class PrimitiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        os.chmod(self.root, 0o700)
        self.uid = os.geteuid()

    def tearDown(self):
        self.temp.cleanup()

    def test_write_new_is_exclusive_sealed_and_single_link(self):
        path = self.root / "receipt.json"
        binding = transition.write_new(
            path, b"evidence\n", owner_uid=self.uid)
        self.assertEqual(binding["sha256"], hashlib.sha256(b"evidence\n").hexdigest())
        self.assertEqual(binding["mode"], 0o444)
        self.assertEqual(binding["nlink"], 1)
        with self.assertRaises(FileExistsError):
            transition.write_new(path, b"replacement\n", owner_uid=self.uid)

    def test_receipt_journal_is_exclusive_canonical_and_hash_chained(self):
        clock = FakeClock()
        deadline = transition.Deadline(540.0, 60.0, clock=clock)
        journal = transition.ReceiptJournal(
            self.root, "test-transition", self.uid, deadline)
        first = journal.record("alpha", "started", {"value": 1})
        second = journal.record("alpha", "completed", {"value": 2})
        first_path = self.root / "0000-alpha-started.json"
        second_path = self.root / "0001-alpha-completed.json"
        self.assertEqual(
            transition.load_canonical(
                first_path, transition.PHASE_SCHEMA, "first phase"),
            first)
        self.assertEqual(
            transition.load_canonical(
                second_path, transition.PHASE_SCHEMA, "second phase"),
            second)
        self.assertEqual(second["previous_receipt_sha256"],
                         transition.file_binding(
                             first_path, with_hash=True)["sha256"])
        chain = journal.replay()
        self.assertEqual(chain["count"], 2)
        self.assertEqual(chain["head_sha256"], transition.file_binding(
            second_path, with_hash=True)["sha256"])
        os.chmod(first_path, 0o644)
        first_path.write_bytes(first_path.read_bytes().replace(b'"value":1',
                                                               b'"value":9'))
        os.chmod(first_path, 0o444)
        with self.assertRaises(transition.TransitionError):
            journal.replay()
        with self.assertRaises(FileExistsError):
            transition.ReceiptJournal(
                self.root, "test-transition", self.uid, deadline).record(
                    "alpha", "started", {})

    def test_receipt_replay_rejects_fifo_before_path_read(self):
        clock = FakeClock()
        journal = transition.ReceiptJournal(
            self.root, "test-transition", self.uid,
            transition.Deadline(540.0, 60.0, clock=clock))
        fifo = self.root / "0000-alpha-started.json"
        os.mkfifo(fifo)
        journal.sequence = 1
        with self.assertRaisesRegex(
                transition.TransitionError, "not a regular file"):
            journal.replay()

    def test_mpstat_receipt_requires_three_finite_exact_pair_intervals(self):
        def raw(intervals):
            return json.dumps({
                "sysstat": {"hosts": [{"statistics": [
                    {"cpu-load": loads} for loads in intervals
                ]}]}
            }).encode("ascii")

        good = [
            [{"cpu": "60", "idle": 99.0},
             {"cpu": "124", "idle": 100.0}]
            for _ in range(3)
        ]
        self.assertEqual(
            transition.parse_mpstat_idle_receipt(raw(good), (60, 124)),
            {60: 99.0, 124: 100.0})
        malformed = {
            "short": good[:2],
            "duplicate": [good[0] + [{"cpu": "60", "idle": 100.0}],
                          good[1], good[2]],
            "missing": [[{"cpu": "60", "idle": 100.0}], good[1], good[2]],
            "nonfinite-string": [
                [{"cpu": "60", "idle": "NaN"},
                 {"cpu": "124", "idle": 100.0}], good[1], good[2]],
            "too-busy": [
                [{"cpu": "60", "idle": 96.0},
                 {"cpu": "124", "idle": 100.0}], good[1], good[2]],
        }
        for name, intervals in malformed.items():
            with self.subTest(name=name):
                with self.assertRaises(transition.TransitionError):
                    transition.parse_mpstat_idle_receipt(
                        raw(intervals), (60, 124))
        nan_token = raw(good).replace(b"99.0", b"NaN", 1)
        with self.assertRaises(transition.TransitionError):
            transition.parse_mpstat_idle_receipt(nan_token, (60, 124))

    def test_rename_noreplace_preserves_exact_binding(self):
        source = self.root / "live.csv"
        source.write_bytes(b"raw\n")
        os.chmod(source, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        destination = self.root / "archive.csv"
        observed = transition.rename_noreplace(
            source, destination, binding, parent_uid=self.uid)
        self.assertEqual(observed, binding)
        self.assertFalse(source.exists())

    def test_rename_noreplace_refuses_collision_without_mutation(self):
        source = self.root / "live.csv"
        destination = self.root / "archive.csv"
        source.write_bytes(b"raw\n")
        destination.write_bytes(b"old\n")
        os.chmod(source, 0o444)
        os.chmod(destination, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        with self.assertRaises(FileExistsError):
            transition.rename_noreplace(
                source, destination, binding, parent_uid=self.uid)
        self.assertEqual(source.read_bytes(), b"raw\n")
        self.assertEqual(destination.read_bytes(), b"old\n")

    def test_rename_noreplace_refuses_changed_source(self):
        source = self.root / "live.csv"
        source.write_bytes(b"raw\n")
        os.chmod(source, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        os.chmod(source, 0o644)
        source.write_bytes(b"bad\n")
        os.chmod(source, 0o444)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "source binding changed"):
            transition.rename_noreplace(
                source, self.root / "archive.csv", binding,
                parent_uid=self.uid)

    def test_file_binding_fifo_is_nonblocking_and_rejected(self):
        fifo = self.root / "fifo"
        os.mkfifo(fifo)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "not a regular file"):
            transition.file_binding(fifo, with_hash=True)

    def test_pidfd_program_opens_descriptor_before_proc_identity_reads(self):
        program = transition.PIDFD_SIGNAL_PROGRAM
        self.assertLess(program.index("pidfd = libc.syscall(434"),
                        program.index("stat_raw = open"))
        self.assertIn("expected_tick", program)
        self.assertIn("expected_cmdline_sha", program)
        self.assertIn("expected_exe_device", program)
        self.assertIn("expected_exe_inode", program)
        self.assertIn("pidfd_send_signal", program)
        self.assertNotIn("os.kill", program)
        compile(program, "<pidfd-old-sampler>", "exec")
        launcher = transition.LAUNCHER_PIDFD_SIGNAL_PROGRAM
        self.assertLess(launcher.index("pidfd = libc.syscall(434"),
                        launcher.index("stat_raw = open"))
        self.assertIn("expected_uids", launcher)
        self.assertIn("expected_exe_device", launcher)
        self.assertIn("expected_exe_inode", launcher)
        self.assertNotIn("os.kill", launcher)
        compile(launcher, "<pidfd-old-launcher>", "exec")

    def test_background_launcher_allows_only_absent_or_exact_session_leader(self):
        proc_root = self.root / "proc"
        proc_root.mkdir()
        self.assertEqual(
            transition.capture_owned_session_leader(
                7001, proc_root=proc_root),
            0)
        leader = proc_root / "7001"
        leader.mkdir()
        (leader / "stat").write_bytes(b"synthetic")
        with mock.patch.object(transition, "_parse_proc_stat", return_value={
                "process_group": 7001, "session": 7001,
                "start_tick": 123}):
            self.assertEqual(
                transition.capture_owned_session_leader(
                    7001, proc_root=proc_root),
                123)
        with mock.patch.object(transition, "_parse_proc_stat", return_value={
                "process_group": 7002, "session": 7001,
                "start_tick": 123}):
            with self.assertRaisesRegex(
                    transition.TransitionError, "exact session"):
                transition.capture_owned_session_leader(
                    7001, proc_root=proc_root)

    def test_prepare_freezes_exact_bytes_and_plan_replays(self):
        candidate = self.root / "candidate.py"
        p32 = self.root / "p32.py"
        controller = self.root / "controller.py"
        candidate.write_bytes(b"candidate\n")
        p32.write_bytes(b"p32\n")
        controller.write_bytes(b"controller\n")
        for path in (candidate, p32, controller):
            os.chmod(path, 0o444)
        plan = replace(
            transition.TransitionPlan(),
            root=self.root / "dry", controller_uid=self.uid,
            candidate_sha256=hashlib.sha256(b"candidate\n").hexdigest(),
            p32_sha256=hashlib.sha256(b"p32\n").hexdigest(),
            controller_sha256=hashlib.sha256(b"controller\n").hexdigest())
        prepared = transition.prepare_transition(
            plan, controller_source=controller,
            candidate_source=candidate, p32_source=p32)
        self.assertEqual(prepared["value"]["transition_id"], plan.transition_id)
        self.assertEqual(prepared["value"]["emergency_recovery_s"], 60.0)
        self.assertEqual(
            prepared["value"]["python_isolation"]["candidate_argv_prefix"],
            ["/usr/bin/python3.12", "-I", "-S", "-B"])
        self.assertEqual(
            set(prepared["value"]["tools"]),
            {"env", "fuser", "mpstat", "python3", "sudo", "taskset",
             "timeout"})
        replay = transition.verify_transition_plan(plan)
        self.assertEqual(replay, prepared["value"])
        with self.assertRaisesRegex(transition.TransitionError,
                                    "plan contract changed"):
            transition.verify_transition_plan(
                replace(plan, candidate_cpu=plan.candidate_cpu + 1))
        with self.assertRaisesRegex(transition.TransitionError,
                                    "old-sampler plan changed"):
            transition.verify_transition_plan(
                replace(plan, old_csv_inode=plan.old_csv_inode + 1))
        changed_tools = json.loads(json.dumps(prepared["value"]["tools"]))
        changed_tools["env"]["binding"]["inode"] += 1
        with mock.patch.object(
                transition, "capture_tool_records",
                return_value=changed_tools):
            with self.assertRaisesRegex(
                    transition.TransitionError, "tool identity changed"):
                transition.verify_transition_plan(plan)

    def test_prepare_rejects_unreviewed_controller_digest_before_creating_root(self):
        candidate = self.root / "candidate-unreviewed.py"
        p32 = self.root / "p32-unreviewed.py"
        controller = self.root / "controller-unreviewed.py"
        for path, raw in ((candidate, b"candidate\n"), (p32, b"p32\n"),
                          (controller, b"controller\n")):
            path.write_bytes(raw)
            os.chmod(path, 0o444)
        dry = self.root / "unreviewed-dry"
        plan = replace(
            transition.TransitionPlan(), root=dry, controller_uid=self.uid,
            candidate_sha256=hashlib.sha256(b"candidate\n").hexdigest(),
            p32_sha256=hashlib.sha256(b"p32\n").hexdigest(),
            controller_sha256="f" * 64)
        with self.assertRaisesRegex(
                transition.TransitionError, "externally reviewed SHA256"):
            transition.prepare_transition(
                plan, controller_source=controller,
                candidate_source=candidate, p32_source=p32)
        self.assertFalse(dry.exists())

    def test_execute_runtime_requires_exact_env_flags_and_orig_argv(self):
        plan = replace(
            transition.TransitionPlan(), root=self.root / "runtime",
            controller_sha256="a" * 64)
        environment = transition.execute_environment(plan)
        orig_argv = transition.expected_execute_orig_argv(plan)
        flags = dict(transition.EXECUTE_FLAG_CONTRACT)
        receipt = transition.verify_execute_runtime(
            plan, observed_environment=environment,
            observed_orig_argv=orig_argv, observed_flags=flags)
        self.assertEqual(receipt["command"][:2], ["/usr/bin/env", "-i"])
        self.assertIn("/usr/bin/python3.12", receipt["command"])
        self.assertEqual(receipt["sys_orig_argv"][1:4], ["-I", "-S", "-B"])
        bad_environment = dict(environment, PYTHONPATH="/attacker")
        with self.assertRaisesRegex(transition.TransitionError, "environment"):
            transition.verify_execute_runtime(
                plan, observed_environment=bad_environment,
                observed_orig_argv=orig_argv, observed_flags=flags)
        bad_flags = dict(flags, no_site=0)
        with self.assertRaisesRegex(transition.TransitionError, "sys.flags"):
            transition.verify_execute_runtime(
                plan, observed_environment=environment,
                observed_orig_argv=orig_argv, observed_flags=bad_flags)
        with self.assertRaisesRegex(transition.TransitionError, "orig_argv"):
            transition.verify_execute_runtime(
                plan, observed_environment=environment,
                observed_orig_argv=orig_argv[:2] + orig_argv[3:],
                observed_flags=flags)

    def test_hardcoded_tool_contract_matches_current_exact_binaries(self):
        records = transition.capture_tool_records()
        self.assertEqual(set(records), {
            "env", "fuser", "mpstat", "python3", "sudo", "taskset",
            "timeout"})
        expected = {name: (path, digest)
                    for name, path, digest in transition.TOOL_CONTRACT}
        for name, record in records.items():
            self.assertEqual(record["path"], expected[name][0])
            self.assertEqual(record["sha256"], expected[name][1])
            self.assertEqual(record["binding"]["sha256"], expected[name][1])
            self.assertEqual(record["binding"]["uid"], 0)
            self.assertFalse(record["binding"]["mode"] & 0o022)

    def test_tool_replay_rejects_same_path_binding_change(self):
        expected = transition.capture_tool_records()
        changed = json.loads(json.dumps(expected))
        changed["taskset"]["binding"]["inode"] += 1
        with mock.patch.object(
                transition, "capture_tool_records", return_value=changed):
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "sealed tool identity changed: taskset"):
                transition.verify_tool_records(expected)
        wrong_python = json.loads(json.dumps(expected))
        wrong_python["python3"]["binding"]["inode"] += 1
        with self.assertRaisesRegex(
                transition.TransitionError,
                "controller interpreter differs"):
            transition.verify_running_interpreter(wrong_python)
        with mock.patch.object(sys, "executable", "/usr/bin/python3"):
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "controller interpreter differs"):
                transition.verify_running_interpreter(
                    expected, require_exact_path=True)

    def test_parse_arguments_requires_both_live_tokens(self):
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                transition.parse_arguments(["--execute-sealed-transition",
                                            transition.TRANSITION_ID])
        args = transition.parse_arguments([
            "--execute-sealed-transition", transition.TRANSITION_ID,
            "--expected-controller-sha256", "a" * 64,
            "--confirmation", transition.EXECUTE_CONFIRMATION,
        ])
        self.assertEqual(args.execute_sealed_transition,
                         transition.TRANSITION_ID)

    def test_thermal_parent_requires_the_expected_controller_owner(self):
        paths = {
            "old_csv": self.root / "thermal.csv",
            "old_pid_file": self.root / "sampler.pid",
            "old_archive": self.root / "archive.csv",
            "old_unclean_archive": self.root / "unclean.csv",
            "old_stale_pid_archive": self.root / "stale.pid",
            "audit_binding": self.root / "audit.json",
        }
        backend = object.__new__(transition.LiveBackend)
        backend.plan = replace(
            transition.TransitionPlan(), controller_uid=self.uid + 1,
            **paths)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "trust boundary"):
            transition.LiveBackend._validate_thermal_parent(backend)
        backend.plan = replace(backend.plan, controller_uid=self.uid)
        os.chmod(self.root, 0o750)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "trust boundary"):
            transition.LiveBackend._validate_thermal_parent(backend)


class ExclusiveReaderAndAuditTests(unittest.TestCase):
    def test_candidate_cleanup_before_start_allows_only_exact_old_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        result = transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(result["i2c_readers_after"], [backend.plan.old_pid])

    def test_candidate_cleanup_allows_exact_old_reader_to_exit_between_probes(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        readers = iter(((backend.plan.old_pid,), ()))
        backend._i2c_readers = lambda: next(readers)
        backend._cleanup_proof_pause = lambda: None
        result = transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(result["i2c_readers_after"], [])

    def test_candidate_cleanup_before_start_rejects_unknown_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        backend._i2c_readers = lambda: (99123,)
        with self.assertRaisesRegex(transition.TransitionError, "unowned"):
            transition.LiveBackend.cleanup_candidate(backend)

    @staticmethod
    def _owned_candidate_backend(*, kill_error=None, observations_empty=True):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.design = {"test": True}
        backend.recovery_budget = None
        backend.candidate_cleanup_complete = False
        backend.p32 = mock.Mock()
        owner = mock.Mock()
        owner.process.pid = 7701
        owner.process.poll.side_effect = [0, 0]
        owner.pid = 7702
        owner.identity = {"pid": owner.pid}
        owner.csv_part = Path("/synthetic/candidate.csv")
        backend.candidate_owner = owner
        if kill_error is not None:
            backend.p32._kill_owned_sampler_session.side_effect = kill_error
        backend.p32._sampler_evidence_paths.return_value = {
            "csv": Path("/synthetic/candidate-final.csv")}
        backend.p32.process_identity_matches.return_value = \
            not observations_empty
        backend._owned_session_members = lambda _session: \
            () if observations_empty else (owner.pid,)
        backend._candidate_pid_present = lambda _pid: \
            not observations_empty
        backend._i2c_readers = lambda: \
            () if observations_empty else (owner.pid,)
        backend._fuser = lambda _path: \
            () if observations_empty else (owner.pid,)
        backend._cleanup_proof_pause = lambda: None
        return backend, owner

    def test_kill_helper_error_with_independent_empty_proof_allows_restore(self):
        helper_error = transition.TransitionError("injected kill transport error")
        backend, _owner = self._owned_candidate_backend(
            kill_error=helper_error, observations_empty=True)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "independently proven empty ownership"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)

    def test_kill_helper_error_without_empty_proof_forbids_second_reader(self):
        helper_error = transition.TransitionError("injected kill transport error")
        backend, _owner = self._owned_candidate_backend(
            kill_error=helper_error, observations_empty=False)
        with self.assertRaisesRegex(
                transition.TransitionError, "independent empty-ownership proof failed"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertFalse(backend.candidate_cleanup_complete)
        backend._quiesce_old_for_recovery = mock.Mock()
        with self.assertRaisesRegex(
                transition.TransitionError, "before candidate cleanup"):
            transition.LiveBackend.restore_old(backend, None)
        backend._quiesce_old_for_recovery.assert_not_called()

    def test_exited_candidate_launcher_still_invokes_exact_session_kill(self):
        backend, owner = self._owned_candidate_backend(
            observations_empty=True)
        result = transition.LiveBackend.cleanup_candidate(backend)
        backend.p32._kill_owned_sampler_session.assert_called_once_with(
            owner, backend.plan.root, backend.design)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(len(result["empty_ownership_observations"]), 2)

    def test_cleanup_probe_separation_survives_expired_emergency_budget(self):
        backend, owner = self._owned_candidate_backend(
            observations_empty=True)
        clock = FakeClock()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(
            1.0, clock=clock)
        clock.value += 2.0
        backend._cleanup_proof_pause = \
            transition.LiveBackend._cleanup_proof_pause.__get__(
                backend, transition.LiveBackend)
        with mock.patch.object(transition.time, "sleep") as sleep:
            result = transition.LiveBackend.cleanup_candidate(backend)
        sleep.assert_called_once_with(0.05)
        backend.p32._kill_owned_sampler_session.assert_called_once_with(
            owner, backend.plan.root, backend.design)
        self.assertEqual(len(result["empty_ownership_observations"]), 2)
        self.assertTrue(backend.candidate_cleanup_complete)

    def test_candidate_accept_replay_rejects_changed_final_artifacts(self):
        backend = object.__new__(transition.LiveBackend)
        backend.candidate_owner = object()
        backend.candidate_cleanup_complete = True
        backend.accept_candidate = lambda: {
            "raw_sha256": "b" * 64, "sample_count": 6}
        with self.assertRaisesRegex(
                transition.TransitionError, "changed after acceptance"):
            transition.LiveBackend._replay_candidate_accept(
                backend, {"raw_sha256": "a" * 64, "sample_count": 6},
                7001)

    def test_candidate_replay_allows_only_restored_old_i2c_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.candidate_owner = object()
        backend.candidate_cleanup_complete = True
        accepted = {"raw_sha256": "a" * 64, "sample_count": 6}
        backend.accept_candidate = lambda: dict(accepted)
        base = {
            "candidate_pid_present": False,
            "csv_readers": [],
            "exact_identity_live": False,
            "i2c_readers": [7001],
            "launcher_returncode": 0,
            "session_members": [],
        }
        backend._candidate_cleanup_observation = lambda _owner: dict(base)
        replay = transition.LiveBackend._replay_candidate_accept(
            backend, accepted, 7001)
        self.assertEqual(replay["ownership"]["i2c_readers"], [7001])
        rejected_replay = transition.LiveBackend._replay_candidate_accept(
            backend, None, 7001)
        self.assertIsNone(rejected_replay["acceptance"])
        self.assertEqual(
            rejected_replay["ownership"]["i2c_readers"], [7001])
        for readers in ([], [7999], [7001, 7999]):
            with self.subTest(readers=readers):
                backend._candidate_cleanup_observation = lambda _owner, r=readers: {
                    **base, "i2c_readers": r}
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "ownership reappeared"):
                    transition.LiveBackend._replay_candidate_accept(
                        backend, accepted, 7001)

    def test_candidate_final_roster_rejects_unbound_extra_artifact(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            segments = root / "segments"
            segments.mkdir()
            paths = {
                "csv": segments / "segment0000.thermal.csv",
                "receipt": segments /
                    "segment0000.thermal.sampler-receipt.json",
                "validation": segments /
                    "segment0000.thermal.validation.jsonl",
            }
            for path in paths.values():
                path.write_bytes(b"sealed\n")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(transition.TransitionPlan(), root=root)
            self.assertEqual(
                set(transition.LiveBackend._candidate_artifact_roster(
                    backend, paths)), set(paths.values()))
            (segments / "unbound-extra").write_bytes(b"extra\n")
            with self.assertRaisesRegex(
                    transition.TransitionError, "roster changed"):
                transition.LiveBackend._candidate_artifact_roster(
                    backend, paths)

    def test_candidate_acceptance_replays_retained_thermal_summary(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            segments = root / "segments"
            segments.mkdir()
            paths = {
                "csv": segments / "segment0000.thermal.csv",
                "receipt": segments /
                    "segment0000.thermal.sampler-receipt.json",
                "validation": segments /
                    "segment0000.thermal.validation.jsonl",
            }
            for path in paths.values():
                path.write_bytes(b"sealed\n")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(transition.TransitionPlan(), root=root)
            backend.p32 = mock.Mock()
            backend.design = {}
            backend.candidate_timing_start = 100.0
            backend.candidate_benchmark_end = 101.0
            backend.candidate_terminal = {
                "thermal_summary": {"cpu_max_c": 1.0}}
            backend.p32._sampler_evidence_paths.return_value = paths
            backend.p32.validate_sampler_terminal_evidence.return_value = (
                {}, b"validation\n")
            backend.p32._parse_thermal_csv.return_value = ({},)
            backend.p32._parse_thermal_validation.return_value = ({},)
            backend.p32.validate_thermal_interval.return_value = {
                "cpu_max_c": 2.0}
            with self.assertRaisesRegex(
                    transition.TransitionError, "summary did not replay"):
                transition.LiveBackend.accept_candidate(backend)

    def test_cleanup_inspection_error_never_authorizes_old_restart(self):
        backend, _owner = self._owned_candidate_backend(
            kill_error=transition.TransitionError("kill helper error"),
            observations_empty=True)
        backend._i2c_readers = mock.Mock(
            side_effect=transition.TransitionError("fuser inspection error"))
        with self.assertRaisesRegex(
                transition.TransitionError, "empty-ownership proof failed"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertFalse(backend.candidate_cleanup_complete)

    def test_old_restart_cleans_exact_session_after_unexpected_capture_error(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            source = root / "old.py"
            source.write_bytes(b"old\n")
            os.chmod(source, 0o444)
            plan = replace(
                transition.TransitionPlan(), root=root,
                old_source=source, old_csv=root / "thermal.csv",
                old_pid_file=root / "sampler.pid",
                old_source_sha256=hashlib.sha256(b"old\n").hexdigest())
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend.design = {"tools": backend.tools}
            backend.p32 = mock.Mock()
            backend.deadline = transition.Deadline(540.0, 60.0)
            backend.recovery_budget = \
                transition.EmergencyRecoveryBudget(60.0)
            backend.old_preflight = {
                "source": transition.file_binding(source, with_hash=True)}
            backend._i2c_readers = lambda: ()
            backend._environment = lambda: {"PATH": "/usr/bin:/bin"}
            backend._require_old_launcher_absent = lambda: None
            process = mock.Mock(pid=7701)
            with mock.patch.object(
                    transition.subprocess, "Popen", return_value=process), \
                    mock.patch.object(
                        transition, "capture_owned_session_leader",
                        side_effect=RuntimeError("injected capture error")):
                with self.assertRaisesRegex(RuntimeError,
                                            "injected capture error"):
                    transition.LiveBackend._launch_old(backend)
            backend.p32._kill_owned_process_session.assert_called_once_with(
                process, 0, mock.ANY, plan.root, backend.design)

    def test_recovery_quiesces_exact_original_launcher_before_restart(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
        backend._old_identity_lives = lambda: False
        backend._require_old_pid_absent = mock.Mock()
        launcher_liveness = iter((True, False, False))
        backend._old_launcher_identity_lives = lambda: next(
            launcher_liveness, False)
        backend._old_launcher_pidfd_signal = mock.Mock()
        backend._require_old_launcher_absent = mock.Mock()
        backend._raw_proc_pid_present = lambda _pid: False
        backend._i2c_readers = lambda: ()
        transition.LiveBackend._quiesce_old_for_recovery(backend)
        backend._old_launcher_pidfd_signal.assert_called_once_with(signal.SIGTERM)
        backend._require_old_launcher_absent.assert_called_once_with()

    def test_recovery_waits_for_raw_zombie_launcher_absence(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
        backend._old_identity_lives = lambda: False
        backend._old_launcher_identity_lives = lambda: False
        backend._require_old_pid_absent = mock.Mock()
        backend._require_old_launcher_absent = mock.Mock()
        launcher_presence = iter((True, False))
        backend._raw_proc_pid_present = lambda pid: (
            False if pid == backend.plan.old_pid else
            next(launcher_presence))
        backend._i2c_readers = lambda: ()
        with mock.patch.object(transition.time, "sleep") as sleep:
            transition.LiveBackend._quiesce_old_for_recovery(backend)
        sleep.assert_called_once_with(0.05)
        backend._require_old_launcher_absent.assert_called_once_with()

    def test_old_sampler_identity_rejects_wrong_executable_inode(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.tools = transition.capture_tool_records()
        identity = {
            "affinity": str(backend.plan.old_cpu),
            "argv": list(backend.plan.old_argv),
            "cmdline_sha256": backend.plan.old_cmdline_sha256,
            "executable": {"device": 1, "inode": 2},
            "pid": backend.plan.old_pid,
            "ppid": backend.plan.old_launcher_pid,
            "process_group": backend.plan.old_process_group,
            "session": backend.plan.old_session,
            "start_tick": backend.plan.old_start_tick,
            "uids": [0, 0, 0, 0],
        }
        with mock.patch.object(
                transition, "_proc_identity", return_value=identity), \
                mock.patch.object(
                    transition.Path, "read_text",
                    return_value=backend.plan.old_boot_id):
            with self.assertRaisesRegex(
                    transition.TransitionError, "changed: executable"):
                transition.LiveBackend._validate_old_child_identity(backend)

    def test_replacement_launcher_roster_rejects_unknown_session_member(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.tools = transition.capture_tool_records()
        backend._owned_session_members = lambda _session: (7001, 7002, 7003)
        command = transition.LiveBackend._replacement_launch_command(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "unknown member"):
            transition.LiveBackend._replacement_launcher_roster(backend, {
                "launcher_command": list(command),
                "launcher_command_sha256": hashlib.sha256(
                    b"\0".join(value.encode("ascii")
                                for value in command) + b"\0").hexdigest(),
                "launcher_session": 7001, "launcher_start_tick": 123,
                "pid": 7002, "start_tick": 456})

    def test_replacement_launcher_roster_rejects_unsealed_command(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.tools = transition.capture_tool_records()
        backend._owned_session_members = mock.Mock(return_value=(7002,))
        command = transition.LiveBackend._replacement_launch_command(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "identity is malformed"):
            transition.LiveBackend._replacement_launcher_roster(backend, {
                "launcher_command": [*command[:-1], "--changed"],
                "launcher_command_sha256": hashlib.sha256(
                    b"\0".join(value.encode("ascii")
                                for value in command) + b"\0").hexdigest(),
                "launcher_session": 7001, "launcher_start_tick": 0,
                "pid": 7002, "start_tick": 456})
        backend._owned_session_members.assert_not_called()

    def test_live_candidate_start_defensively_rejects_forced_archive(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.stop_record = {"forced": False}
        backend.archive_record = {"forced_stop": True}
        backend._i2c_readers = lambda: ()
        with self.assertRaisesRegex(
                transition.TransitionError, "lacks an empty-I2C archive seal"):
            transition.LiveBackend.start_candidate(backend)

    def test_csv_reader_preflight_rejects_any_extra_pid(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        backend._fuser = lambda _path: (backend.plan.old_pid, 919688)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "auditor or unknown"):
            transition.LiveBackend._require_exclusive_old_readers(backend)

    def test_csv_reader_preflight_accepts_only_old_writer(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        backend._fuser = lambda _path: (backend.plan.old_pid,)
        self.assertEqual(
            transition.LiveBackend._require_exclusive_old_readers(backend),
            {"csv_readers": [backend.plan.old_pid],
             "i2c_readers": [backend.plan.old_pid]})

    def test_final_old_identity_replays_owned_session_and_process_group(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            source = root / "old.py"
            csv_path = root / "thermal.csv"
            pid_path = root / "sampler.pid"
            source.write_bytes(b"old source\n")
            csv_path.write_bytes(b"sample\n")
            pid_path.write_bytes(b"7001\n")
            for path in (source, csv_path, pid_path):
                os.chmod(path, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_source=source, old_csv=csv_path, old_pid_file=pid_path,
                old_source_sha256=hashlib.sha256(b"old source\n").hexdigest())
            identity = {
                "boot_id": "boot",
                "cmdline": list(plan.replacement_old_argv),
                "cmdline_sha256": plan.replacement_old_cmdline_sha256,
                "pid": 7001,
                "process_group": 9000, "session_id": 9000,
                "start_tick": 123,
            }
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.p32 = mock.Mock()
            backend.old_preflight = {
                "source": transition.file_binding(source, with_hash=True)}
            backend.p32.capture_process_identity.return_value = identity
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "%.6f" % transition.time.monotonic()},)
            backend._i2c_readers = lambda: (7001,)
            backend._fuser = lambda _path: (7001,)
            csv_binding = transition.file_binding(csv_path, with_hash=False)
            restored = {
                "boot_id": "boot",
                "cmdline_sha256": plan.replacement_old_cmdline_sha256,
                "csv_initial_size": csv_binding["size"],
                "csv_live_identity": {
                    key: csv_binding[key] for key in (
                        "device", "gid", "inode", "mode", "nlink", "uid")},
                "launcher_session": 9000,
                "pid": 7001,
                "pid_binding": transition.file_binding(pid_path, with_hash=True),
                "start_tick": 123,
            }
            backend._replacement_launcher_roster = lambda _restored: [
                {"pid": 7001}]
            result = transition.LiveBackend._revalidate_restored_old(
                backend, restored)
            self.assertEqual(result["identity"], identity)
            backend.p32._parse_thermal_csv.return_value = ()
            with self.assertRaisesRegex(
                    transition.TransitionError, "CSV content changed"):
                transition.LiveBackend._revalidate_restored_old(
                    backend, restored)
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "%.6f" % transition.time.monotonic()},)
            bad_identity = dict(identity)
            bad_identity["process_group"] = 9001
            backend.p32.capture_process_identity.return_value = bad_identity
            with self.assertRaisesRegex(
                    transition.TransitionError, "changed before audit"):
                transition.LiveBackend._revalidate_restored_old(
                    backend, restored)

    def test_future_audit_binding_preserves_archive_sha_and_new_inode(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            archive_binding = transition.file_binding(archive, with_hash=True)
            archive_binding["uid"] = 0
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "csv_live_identity": {"device": 7, "inode": 99},
                "archived_pre_dry_sha256": archive_binding["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            real_file_binding = transition.file_binding

            def root_archive_binding(path, *, with_hash):
                binding = real_file_binding(path, with_hash=with_hash)
                if path == archive:
                    binding["uid"] = 0
                return binding

            with mock.patch.object(
                    transition, "file_binding", side_effect=root_archive_binding), \
                    mock.patch.object(
                        transition, "verify_running_interpreter",
                        return_value={"device": 1, "inode": 2,
                                      "argv_path": "/sealed/python",
                                      "resolved_path": "/sealed/python"}):
                result = transition.LiveBackend.publish_audit_binding(
                    backend,
                    {"path": str(archive), "binding": archive_binding},
                    restored, False, None,
                    {"count": 1, "head_sha256": "c" * 64,
                     "roster": [{"sequence": 0}]})
            value = transition.load_canonical(
                plan.audit_binding, transition.AUDIT_BINDING_SCHEMA,
                "test audit binding")
            self.assertEqual(value, result["value"])
            self.assertEqual(value["tools"], backend.tools)
            self.assertEqual(value["receipt_chain_prefix"]["head_sha256"],
                             "c" * 64)
            self.assertEqual(result["tools_after_audit"], backend.tools)
            self.assertEqual(
                value["archived_pre_dry"]["binding"]["sha256"],
                archive_binding["sha256"])
            self.assertEqual(value["live_old_sampler"]
                             ["csv_live_identity"]["inode"],
                             99)

    def test_audit_rejects_success_without_retained_candidate_acceptance(self):
        backend = object.__new__(transition.LiveBackend)
        backend._replay_external_state = mock.Mock()
        with self.assertRaisesRegex(
                transition.TransitionError, "lacks retained candidate"):
            transition.LiveBackend.publish_audit_binding(
                backend, None, {}, True, None,
                {"count": 1, "head_sha256": "c" * 64,
                 "roster": [{"sequence": 0}]})
        backend._replay_external_state.assert_not_called()

    def test_future_audit_binding_rehashes_the_archive(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            stale = transition.file_binding(archive, with_hash=True)
            os.chmod(archive, 0o644)
            archive.write_bytes(b"changed!\n")
            os.chmod(archive, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "archived_pre_dry_sha256": stale["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            with self.assertRaisesRegex(
                    transition.TransitionError, "archive changed"):
                transition.LiveBackend.publish_audit_binding(
                    backend, {"path": str(archive), "binding": stale},
                    restored, False, None,
                    {"count": 1, "head_sha256": "c" * 64,
                     "roster": [{"sequence": 0}]})
            self.assertFalse(plan.audit_binding.exists())

    def test_future_audit_binding_fails_if_post_write_tool_replay_changes(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            archive_binding = transition.file_binding(archive, with_hash=True)
            archive_binding["uid"] = 0
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            tools = transition.capture_tool_records()
            changed_tools = json.loads(json.dumps(tools))
            changed_tools["sudo"]["binding"]["inode"] += 1
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = tools
            backend._verify_tools = lambda: tools
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "archived_pre_dry_sha256": archive_binding["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            real_file_binding = transition.file_binding

            def root_archive_binding(path, *, with_hash):
                binding = real_file_binding(path, with_hash=with_hash)
                if path == archive:
                    binding["uid"] = 0
                return binding

            with mock.patch.object(
                    transition, "file_binding",
                    side_effect=root_archive_binding), \
                    mock.patch.object(
                        transition, "verify_running_interpreter",
                        return_value={"device": 1, "inode": 2,
                                      "argv_path": "/sealed/python",
                                      "resolved_path": "/sealed/python"}), \
                    mock.patch.object(
                        transition, "verify_tool_records",
                        return_value=changed_tools):
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "tool identities changed after audit"):
                    transition.LiveBackend.publish_audit_binding(
                        backend,
                        {"path": str(archive), "binding": archive_binding},
                        restored, False, None,
                        {"count": 1, "head_sha256": "c" * 64,
                         "roster": [{"sequence": 0}]})
            self.assertTrue(plan.audit_binding.exists())

    def test_final_replay_rejects_postpublication_audit_replacement(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            audit_path = root / "future.json"
            value = transition.sealed_record(
                transition.AUDIT_BINDING_SCHEMA, {
                    "candidate_accept": None,
                    "prepublication_replay": {"state": "clean"},
                })
            binding = transition.write_new(
                audit_path, transition.canonical_json(value),
                owner_uid=os.geteuid())
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                audit_binding=audit_path)
            backend._replay_external_state = lambda *_args: {"state": "clean"}
            audit = {"binding": binding, "path": str(audit_path),
                     "value": value}
            replay = transition.LiveBackend.final_replay(
                backend, None, None, {}, audit)
            self.assertEqual(replay["audit"]["binding"], binding)
            replacement = transition.sealed_record(
                transition.AUDIT_BINDING_SCHEMA, {
                    "candidate_accept": None,
                    "prepublication_replay": {"state": "changed"},
                })
            os.chmod(audit_path, 0o644)
            audit_path.write_bytes(transition.canonical_json(replacement))
            os.chmod(audit_path, 0o444)
            with self.assertRaisesRegex(
                    transition.TransitionError, "changed after publication"):
                transition.LiveBackend.final_replay(
                    backend, None, None, {}, audit)

    def test_recovery_discovers_archive_after_post_rename_exception(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv", old_pid_file=root / "sampler.pid",
                old_archive=archive,
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            backend._fuser = lambda _path: ()
            record = transition.LiveBackend._ensure_recovery_archive(
                backend, None)
            self.assertEqual(record["path"], str(archive))
            self.assertEqual(record["binding"]["sha256"],
                             hashlib.sha256(b"pre-dry\n").hexdigest())

    def test_recovery_rejects_archive_with_wrong_preflight_inode(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            binding = transition.file_binding(archive, with_hash=True)
            expected = dict(binding)
            expected["inode"] += 1
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv", old_pid_file=root / "sampler.pid",
                old_archive=archive,
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.old_preflight = {"csv": expected}
            backend.archive_record = None
            with self.assertRaisesRegex(
                    transition.TransitionError, "preflight CSV inode"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)

    def test_recovery_rejects_alternate_or_unbound_stale_archive_names(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            clean = root / "clean.csv"
            unclean = root / "unclean.csv"
            stale = root / "stale.pid"
            for path, raw in ((clean, b"pre-dry\n"),
                              (unclean, b"other\n")):
                path.write_bytes(raw)
                os.chmod(path, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv",
                old_pid_file=root / "sampler.pid", old_archive=clean,
                old_unclean_archive=unclean,
                old_stale_pid_archive=stale)
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            with self.assertRaisesRegex(
                    transition.TransitionError, "multiple pre-dry archives"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)
            unclean.unlink()
            stale.write_bytes(b"123\n")
            os.chmod(stale, 0o444)
            with self.assertRaisesRegex(
                    transition.TransitionError, "stale PID evidence"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)

    def test_recovery_archives_canonical_csv_and_stale_pid(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            csv_path = root / "thermal.csv"
            pid_path = root / "sampler.pid"
            csv_path.write_bytes(b"pre-dry\n")
            pid_path.write_bytes(b"123\n")
            os.chmod(csv_path, 0o444)
            os.chmod(pid_path, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=csv_path, old_pid_file=pid_path,
                old_archive=root / "archive.csv",
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            backend._fuser = lambda _path: ()
            record = transition.LiveBackend._ensure_recovery_archive(
                backend, None)
            self.assertEqual(record["path"], str(plan.old_unclean_archive))
            self.assertEqual(record["stale_pid"]["path"],
                             str(plan.old_stale_pid_archive))
            self.assertFalse(csv_path.exists())
            self.assertFalse(pid_path.exists())


if __name__ == "__main__":
    unittest.main()
