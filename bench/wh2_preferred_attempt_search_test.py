#!/usr/bin/env python3
"""Offline regression tests for the preferred-attempt holdout controller."""

from __future__ import annotations

import argparse
from contextlib import ExitStack, redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import tempfile
import unittest
from unittest import mock

import wh2_preferred_attempt_campaign as campaign_module
import wh2_preferred_attempt_search as subject
import wh2_rank_floor_two_anchor_screen as screen_module


class PreferredAttemptSearchTest(unittest.TestCase):
    def fake_repo(self, parent: Path) -> tuple[Path, str]:
        repo = parent / "source"
        remote = parent / "origin.git"
        subprocess.run(("git", "init", "-q", "--bare", str(remote)), check=True)
        subprocess.run(("git", "init", "-q", str(repo)), check=True)
        subprocess.run(
            ("git", "-C", str(repo), "config", "user.name", "WH2 Test"),
            check=True,
        )
        subprocess.run(
            (
                "git", "-C", str(repo), "config", "user.email",
                "wh2-test@example.invalid",
            ),
            check=True,
        )
        subprocess.run(
            ("git", "-C", str(repo), "commit", "--allow-empty", "-qm", "fixture"),
            check=True,
        )
        subprocess.run(
            ("git", "-C", str(repo), "remote", "add", "origin", str(remote)),
            check=True,
        )
        commit = subprocess.check_output(
            ("git", "-C", str(repo), "rev-parse", "HEAD"), text=True,
        ).strip()
        return repo, commit

    def fake_contract(self, repo: Path, commit: str) -> dict[str, object]:
        tools: dict[str, dict[str, str]] = {}
        for name in ("git", "ssh", "ssh-add"):
            executable = Path(shutil.which(name) or f"/usr/bin/{name}").resolve(
                strict=True)
            tools[name] = {
                "path": str(executable),
                "sha256": subject.common.sha256_file(executable),
            }
        return {
            "source_repo": str(repo),
            "source_commit": commit,
            "protocol_sha256": "1" * 64,
            "binary_sha256": "2" * 64,
            "drand_wrapper_sha256": "3" * 64,
            "drand_client_bundle_sha256": "4" * 64,
            "drand_package_sha256": "5" * 64,
            "drand_package_lock_sha256": "6" * 64,
            "node": {
                "path": "/usr/bin/node", "version": subject.DRAND_NODE_VERSION,
                "sha256": "7" * 64,
            },
            "npm": {
                "path": "/usr/bin/npm", "version": subject.DRAND_NPM_VERSION,
                "sha256": "8" * 64,
            },
            "system_tools": tools,
        }

    def required_h1_tree(self, result: Path) -> None:
        for relative in subject.HOLDOUT_PHASES["h1"].required_files:
            path = result / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((relative + "\n").encode("ascii"))

    def latest_bound(self, latest_round: int = 1) -> dict[str, object]:
        origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
        return {
            "schema": "wirehair.wh2.verified_latest_bound.v1",
            "latest_round": latest_round,
            "canonical_beacon_sha256": "9" * 64,
            "consensus_origins": origins[:3],
            "raw_response_sha256": {
                origin: f"{index + 1:x}" * 64
                for index, origin in enumerate(origins[:3])
            },
            "wrapper_result_sha256": "a" * 64,
            "verified_ms": subject.now_utc_ms(),
        }

    def test_quicknet_formula_and_root_goldens(self) -> None:
        subject.validate_beacon_contract_kats()
        canonical = subject.canonical_beacon_bytes(subject.DRAND_KNOWN_ROUND)
        self.assertEqual(
            hashlib.sha256(canonical).hexdigest(),
            subject.DRAND_CANONICAL_KNOWN_SHA256,
        )
        self.assertEqual(
            subject.holdout_master_root("0" * 64, subject.DRAND_KNOWN_ROUND),
            subject.DRAND_MASTER_ROOT_KNOWN_SHA256,
        )
        mutated = dict(subject.DRAND_KNOWN_ROUND)
        mutated["randomness"] = "0" * 64
        with self.assertRaises(subject.CampaignError):
            subject.canonical_beacon_bytes(mutated)
        extended = dict(subject.DRAND_KNOWN_ROUND)
        extended["unexpected"] = "ignored entropy"
        with self.assertRaises(subject.CampaignError):
            subject.canonical_beacon_bytes(extended)
        seal_ms = subject.DRAND_GENESIS_TIME * 1000 + 1_000_000
        local_round = subject.target_beacon_round(seal_ms)
        latest_round = local_round + 100
        record, _encoded = subject.build_seal_record(
            "h1", "b" * 64,
            self.fake_contract(Path("/tmp"), "c" * 40),
            "git@github.com:fake/wirehair.git", seal_ms,
            self.latest_bound(latest_round),
        )
        self.assertEqual(record["round"], latest_round + 40)
        reseal, _encoded = subject.build_seal_record(
            "h1", "b" * 64,
            self.fake_contract(Path("/tmp"), "c" * 40),
            "git@github.com:fake/wirehair.git", seal_ms,
            self.latest_bound(latest_round), record["round"],
        )
        self.assertEqual(reseal["round"], record["round"] + 1)

    def test_verified_latest_is_only_a_hashed_round_bound(self) -> None:
        beacon = dict(subject.DRAND_KNOWN_ROUND)
        beacon["round"] = 123456
        origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
        canonical_sha256 = hashlib.sha256(
            subject.canonical_beacon_bytes(beacon)).hexdigest()
        observations = [{
            "origin": origin,
            "ok": True,
            "beacon": beacon,
            "canonical_sha256": canonical_sha256,
            "raw_response_sha256": f"{index + 1:x}" * 64,
        } for index, origin in enumerate(origins)]
        wave = {
            "schema": "wirehair.wh2.drand_quicknet_latest_wave.v1",
            "role": "preseal-clock-liveness-lower-bound-only",
            "chain_hash": subject.DRAND_CHAIN_HASH,
            "status": "QUORUM",
            "latest_round": beacon["round"],
            "beacon": beacon,
            "consensus_origins": origins[:3],
            "observations": observations,
        }
        completed = subprocess.CompletedProcess(
            ("node", "wrapper", "bundle", "latest-bound"), 0,
            json.dumps(wave).encode("ascii"), b"",
        )
        contract = {"node": {"path": "/usr/bin/node"}}
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(
                subject.subprocess, "run", return_value=completed))
            stack.enter_context(mock.patch.object(
                subject, "offline_verify_beacon", return_value={}))
            bound = subject.obtain_verified_latest_bound(
                contract, Path("/does/not/matter"))
        self.assertEqual(bound["latest_round"], beacon["round"])
        self.assertNotIn("beacon", bound)
        self.assertNotIn("randomness", json.dumps(bound))
        no_quorum = dict(wave)
        no_quorum["status"] = "TEMPORARY_NO_QUORUM"
        failed = subprocess.CompletedProcess(
            completed.args, 2, json.dumps(no_quorum).encode("ascii"), b"")
        with mock.patch.object(subject.subprocess, "run", return_value=failed):
            with self.assertRaises(subject.CampaignError):
                subject.obtain_verified_latest_bound(
                    contract, Path("/does/not/matter"))
        malformed = subprocess.CompletedProcess(
            completed.args, 0, b"\xff", b"")
        with mock.patch.object(subject.subprocess, "run", return_value=malformed):
            with self.assertRaises(subject.CampaignError):
                subject.obtain_verified_latest_bound(
                    contract, Path("/does/not/matter"))

    def test_agent_refresh_uses_frozen_absolute_ssh_tools(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="ssh-wh2-agent-") as temporary:
            root = Path(temporary)
            repo, commit = self.fake_repo(root)
            tools = self.fake_contract(repo, commit)["system_tools"]
            socket_path = root / "agent.999999"
            listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            listener.bind(str(socket_path))
            try:
                listed = subprocess.CompletedProcess(
                    (), 0, "2048 SHA256:test forwarded-key\n", "")
                authenticated = subprocess.CompletedProcess(
                    (), 1, "", "Hi! You've successfully authenticated.\n")
                with mock.patch.object(
                        subject.Path, "glob", return_value=[socket_path]), \
                     mock.patch.object(
                         subject.subprocess, "run",
                         side_effect=(listed, authenticated)) as run:
                    environment = subject.newest_github_agent_environment(tools)
                ssh = str(tools["ssh"]["path"])
                ssh_add = str(tools["ssh-add"]["path"])
                self.assertEqual(run.call_count, 2)
                self.assertEqual(run.call_args_list[0].args[0][0], ssh_add)
                self.assertEqual(run.call_args_list[1].args[0][0], ssh)
                self.assertEqual(environment["SSH_AUTH_SOCK"], str(socket_path))
                self.assertEqual(
                    environment["GIT_SSH_COMMAND"], subject.shlex.join((
                        ssh, "-o", "BatchMode=yes", "-o",
                        "ConnectTimeout=10")))
                self.assertEqual(environment["GIT_SSH_VARIANT"], "ssh")
                self.assertNotIn("GIT_SSH", environment)
            finally:
                listener.close()

    def test_fake_remote_tag_and_late_ordering(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-tag-test-") as temporary:
            repo, commit = self.fake_repo(Path(temporary))
            contract = self.fake_contract(repo, commit)
            tools = contract["system_tools"]
            environment = os.environ.copy()
            poison = Path(temporary) / "poison-path"
            poison.mkdir()
            fake_git = poison / "git"
            fake_git.write_text("#!/bin/sh\nexit 97\n", encoding="ascii")
            fake_git.chmod(0o755)
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=environment))
                stack.enter_context(mock.patch.dict(os.environ, {
                    "PATH": str(poison) + os.pathsep +
                    environment.get("PATH", ""),
                }))
                seal = "1" * 64
                record = {
                    "round": 123,
                    "round_time_ms": subject.now_utc_ms() + 60_000,
                    "source_commit": commit,
                }
                publication = subject.publish_seal_tag(
                    repo, "h1", seal, "2" * 64, record, tools)
                repeated = subject.publish_seal_tag(
                    repo, "h1", seal, "2" * 64, record, tools)
                self.assertEqual(publication["status"], "CONFIRMED")
                self.assertEqual(publication["tag_object"], repeated["tag_object"])
                late_seal = "3" * 64
                late_record = {
                    "round": 124,
                    "round_time_ms": subject.now_utc_ms() - 1,
                    "source_commit": commit,
                }
                with self.assertRaises(subject.CampaignError):
                    subject.publish_seal_tag(
                        repo, "h2", late_seal, "4" * 64, late_record, tools)
                late_name = (
                    "wh2-h12-kpreferred-h2-seal-" + late_seal[:16]
                )
                _raw, rows = subject.tag_remote_rows(
                    repo, late_name, environment, tools)
                self.assertEqual(rows, {})
                self.assertEqual(subprocess.run(
                    (
                        str(tools["git"]["path"]), "-C", str(repo),
                        "show-ref", "--verify", "--quiet",
                        f"refs/tags/{late_name}",
                    ),
                    check=False,
                ).returncode, 0)

                freeze = Path(temporary) / "freeze"
                (freeze / "frozen").mkdir(parents=True)
                prepare = {
                    "source_commit": commit,
                    "binary_sha256": "1" * 64,
                    "script_sha256": "2" * 64,
                    "contract_sha256": "3" * 64,
                    "protocol_sha256": "4" * 64,
                    "route_context_sha256": "5" * 64,
                    "staged_seal_sha256": "6" * 64,
                    "development_seed_ledger_sha256": "7" * 64,
                    "future_beacon_contract_sha256": "8" * 64,
                    "github_remote": str(Path(temporary) / "origin.git"),
                }
                subject.common.atomic_json(freeze / "prepare.json", prepare)
                prepare_sha = subject.common.sha256_file(freeze / "prepare.json")
                freeze_publication = subject.publish_r1_freeze_tag(
                    repo, prepare, prepare_sha, tools)
                publication_path = freeze / "frozen/r1_freeze_publication.json"
                subject.atomic_state_json(publication_path, freeze_publication)
                freeze_contract = {
                    "source_repo": str(repo), "source_commit": commit,
                    "system_tools": tools,
                }
                subject.verify_r1_freeze_publication(
                    freeze, prepare, freeze_contract)
                publication_bytes = publication_path.read_bytes()
                publication_path.unlink()
                with self.assertRaises(subject.CampaignError):
                    subject.verify_r1_freeze_publication(
                        freeze, prepare, freeze_contract)
                subject.common.atomic_write(publication_path, publication_bytes)
                mutated_prepare = dict(prepare)
                mutated_prepare["staged_seal_sha256"] = "9" * 64
                subject.common.atomic_json(
                    freeze / "prepare.json", mutated_prepare)
                with self.assertRaises(subject.CampaignError):
                    subject.verify_r1_freeze_publication(
                        freeze, mutated_prepare, freeze_contract)
                mutated_sha = subject.common.sha256_file(freeze / "prepare.json")
                replacement_publication = subject.publish_r1_freeze_tag(
                    repo, mutated_prepare, mutated_sha, tools)
                subject.atomic_state_json(
                    publication_path, replacement_publication)
                subject.verify_r1_freeze_publication(
                    freeze, mutated_prepare, freeze_contract)

    def test_staged_seal_is_immutable_root_of_trust(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-staged-anchor-") as temporary:
            result = Path(temporary)
            frozen = result / "frozen"
            frozen.mkdir()
            subject.common.atomic_json(frozen / "contract.json", {"version": 1})
            subject.common.atomic_write(frozen / "binary", b"original\n")
            staged_files = [frozen / "contract.json", frozen / "binary"]
            subject.common.atomic_write(
                frozen / "staged.sha256",
                subject.common.sha_manifest(result, staged_files),
            )
            prepare = {
                "staged_seal_sha256":
                    subject.common.sha256_file(frozen / "staged.sha256"),
                "contract_sha256":
                    subject.common.sha256_file(frozen / "contract.json"),
            }
            subject.verify_frozen_staged_anchor(result, prepare)
            alias = frozen / "binary-alias"
            (frozen / "binary").rename(alias)
            (frozen / "binary").symlink_to(alias.name)
            with self.assertRaises(subject.CampaignError):
                subject.verify_frozen_staged_anchor(result, prepare)
            (frozen / "binary").unlink()
            alias.rename(frozen / "binary")
            subject.verify_frozen_staged_anchor(result, prepare)
            subject.common.atomic_json(frozen / "contract.json", {"version": 2})
            subject.common.atomic_write(frozen / "binary", b"replacement\n")
            subject.common.atomic_write(
                frozen / "staged.sha256",
                subject.common.sha_manifest(result, staged_files),
            )
            with self.assertRaises(subject.CampaignError):
                subject.verify_frozen_staged_anchor(result, prepare)
            q0 = {"passed": True, "comparisons": 1}
            subject.common.atomic_json(frozen / "q0_identity.json", q0)
            q0_sha = subject.common.sha256_file(frozen / "q0_identity.json")
            contract = {
                "source_commit": "1" * 40,
                "binary_sha256": "2" * 64,
                "script_sha256": "3" * 64,
                "allk_sha256": "4" * 64,
                "helper_sha256": "5" * 64,
                "protocol_sha256": "6" * 64,
                "route_context_sha256": "7" * 64,
                "q0_identity_sha256": q0_sha,
                "drand_wrapper_sha256": "8" * 64,
                "drand_package_sha256": "9" * 64,
                "drand_package_lock_sha256": "a" * 64,
                "drand_client_bundle_sha256": "b" * 64,
                "node": {"version": subject.DRAND_NODE_VERSION},
                "npm": {"version": subject.DRAND_NPM_VERSION},
                "drand_offline_selftest": {"passed": True},
            }
            bound_prepare = {**contract, "q0_identity": q0}
            subject.verify_prepare_contract_bindings(
                frozen, bound_prepare, contract)
            for field in contract:
                mutated = dict(bound_prepare)
                mutated[field] = None
                with self.assertRaises(subject.CampaignError, msg=field):
                    subject.verify_prepare_contract_bindings(
                        frozen, mutated, contract)
            mutated_q0 = dict(bound_prepare)
            mutated_q0["q0_identity"] = {"passed": False}
            with self.assertRaises(subject.CampaignError):
                subject.verify_prepare_contract_bindings(
                    frozen, mutated_q0, contract)

    def test_phase_manifests_exclude_only_operational_locks(self) -> None:
        for phase, spec in subject.HOLDOUT_PHASES.items():
            with self.subTest(phase=phase), tempfile.TemporaryDirectory(
                    prefix=f"wh2-{phase}-manifest-locks-") as temporary:
                result = Path(temporary)
                for relative in spec.required_files:
                    path = result / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes((relative + "\n").encode("ascii"))
                before = subject.phase_manifest_bytes(result, spec)
                operational = (
                    result / "locks/development-control.lock",
                    result / "locks/h1-paired.lock",
                    result / "beacon/h1.lock",
                    result / "beacon/h2.lock",
                )
                for path in operational:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b"")
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                operational[0].write_bytes(b"operational metadata\n")
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                future_outputs = tuple(
                    relative for relative in spec.forbidden_before_seal
                    if relative == "timing" or relative.startswith("thermal/"))
                self.assertTrue(future_outputs)
                for relative in future_outputs:
                    path = result / relative
                    if relative == "timing":
                        path = path / "sealed-report.json"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b"future phase output\n")
                    self.assertTrue(any(
                        subject.is_relative_prefix(
                            path.relative_to(result).as_posix(), prefix)
                        for prefix in spec.excluded_prefixes))
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                (result / "unexpected.txt").write_bytes(b"not a lock\n")
                with self.assertRaises(subject.CampaignError):
                    subject.phase_manifest_bytes(result, spec)

    def test_development_preseal_requires_exact_terminal_manifests(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-development-readiness-") as temporary:
            result = Path(temporary)
            phase_files: set[Path] = set()
            for top in ("jobs", "ledgers", "completed", "thermal"):
                path = result / top / "artifact.bin"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes((top + "\n").encode("ascii"))
                phase_files.add(path)
            development = result / "development"
            development.mkdir()
            table_path = development / "frozen_table.json"
            table_path.write_bytes(b"{}\n")
            phase_path = development / "phase_complete.sha256"
            phase_path.write_bytes(subject.common.sha_manifest(
                result, phase_files))
            analysis_path = development / "analysis_complete.sha256"
            analysis_path.write_bytes(subject.common.sha_manifest(
                result, {table_path}))

            route_context = "7" * 64
            contract = {
                "source_commit": "1" * 40,
                "protocol_sha256": "2" * 64,
                "route_context_sha256": route_context,
            }
            prepare = {"staged_seal_sha256": "3" * 64}
            table_record = {
                "source_commit": contract["source_commit"],
                "staged_seal_sha256": prepare["staged_seal_sha256"],
                "protocol_sha256": contract["protocol_sha256"],
                "route_context_sha256": route_context,
                "phase_seal_sha256": subject.common.sha256_file(phase_path),
            }
            fake_campaign = mock.Mock()
            fake_campaign.load_canonical_object.return_value = table_record

            def verify() -> None:
                with ExitStack() as stack:
                    stack.enter_context(mock.patch.object(
                        subject, "load_frozen_python_module",
                        side_effect=(fake_campaign, object())))
                    stack.enter_context(mock.patch.object(
                        subject, "load_frozen_cohort_bins",
                        return_value=((4096,), {0: (4096,)})))
                    stack.enter_context(mock.patch.object(
                        subject, "load_development_table_and_caches",
                        return_value=({4096: 1}, {0: object()}, "4" * 64)))
                    subject.verify_development_seal_readiness(
                        result, prepare, contract)

            verify()
            (development / "late-artifact.bin").write_bytes(b"late\n")
            with self.assertRaises(subject.CampaignError):
                verify()

    def test_h1_terminal_record_and_h2_required_files_are_exact(self) -> None:
        record = subject.h1_controller_complete_record(
            campaign_module,
            phase_completion_sha256="1" * 64,
            analysis_complete_sha256="2" * 64,
            root_file_sha256="3" * 64,
            parent_table_sha256="4" * 64,
            frozen_table_sha256="5" * 64,
            route_cache_index_sha256="6" * 64,
            route_context_sha256="7" * 64,
            input_routed_K=123,
            retained_routed_K=45,
        )
        campaign_module.verify_sealed_record(
            record,
            "wirehair.wh2.h12_preferred_attempt.h1_controller_complete.v1")
        self.assertEqual(record["cache_count"], subject.BIN_COUNT)
        required = set(subject.HOLDOUT_PHASES["h2"].required_files)
        self.assertTrue({
            "h1/controller_complete.json",
            "h1/controller_complete.json.sha256",
            "h1/route_cache/index.json",
            "h1/route_cache/index.json.sha256",
        }.issubset(required))

    def test_distinct_cpu_and_dimm_thermal_limits(self) -> None:
        def row(monotonic: int, cpu_c: int, dimm_c: int) -> bytes:
            values = [
                "2026-07-18T00:00:00Z", str(monotonic), "100", "4000",
                str(cpu_c), *([str(dimm_c)] * 8), "0", "1", "1", "1",
                "0", "0",
            ]
            self.assertEqual(len(values), len(subject.THERMAL_FIELDS))
            return (",".join(values) + "\n").encode("ascii")

        with tempfile.TemporaryDirectory(
                prefix="wh2-separate-thermal-limits-") as temporary:
            root = Path(temporary)
            thermal = root / "thermal.csv"
            header = (",".join(subject.THERMAL_FIELDS) + "\n").encode("ascii")
            baseline_row = row(1, 84, 89)
            interval_row = row(2, 85, 89)
            thermal.write_bytes(header + baseline_row + interval_row)
            metadata = thermal.stat()
            baseline = screen_module.parse_thermal_sample(
                baseline_row, "test baseline")
            mark = {
                "dev": metadata.st_dev,
                "ino": metadata.st_ino,
                "offset": len(header) + len(baseline_row),
                "header": header,
                "baseline_row": baseline_row,
                "baseline": baseline,
            }
            with mock.patch.object(
                    screen_module, "validate_thermal_current",
                    return_value=None):
                separate = subject.thermal_finish(
                    thermal, mark, root / "separate.csv", limit_c=85.0,
                    dimm_limit_c=90.0)
                legacy_equal = subject.thermal_finish(
                    thermal, mark, root / "legacy-equal.csv", limit_c=85.0)
            self.assertEqual(separate["thermal_high_samples"], 1)
            self.assertEqual(
                separate["thermal_high_max_consecutive_samples"], 1)
            self.assertEqual(legacy_equal["thermal_high_samples"], 2)
            self.assertEqual(separate["thermal_limit_c"], 85.0)
            with self.assertRaises(ValueError):
                subject.thermal_finish(
                    thermal, mark, root / "invalid.csv", limit_c=85.0,
                    dimm_limit_c=float("nan"))

    def test_timing_host_defers_signals_until_mocked_restore_finishes(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        session = subject.TimingHostSession({}, timing, 0)
        previous = {
            signum: signal.getsignal(signum)
            for signum in (signal.SIGINT, signal.SIGTERM)
        }
        session._install_signal_handlers()
        try:
            with self.assertRaises(KeyboardInterrupt):
                session._handle_timing_signal(signal.SIGINT, None)
            with self.assertRaises(subject.CampaignError):
                session._handle_timing_signal(signal.SIGTERM, None)
            session.sibling_states = {64: "1"}
            session.slice_states = {"system.slice": ("", "0 64")}
            session.active = True

            def restart_fillers() -> None:
                handler = signal.getsignal(signal.SIGTERM)
                self.assertIs(
                    getattr(handler, "__self__", None), session)
                self.assertIs(
                    getattr(handler, "__func__", None),
                    subject.TimingHostSession._defer_timing_signal)
                handler(signal.SIGTERM, None)

            with ExitStack() as stack:
                read_text = stack.enter_context(mock.patch.object(
                    subject.Path, "read_text", side_effect=("0", "1")))
                write_online = stack.enter_context(mock.patch.object(
                    session, "_write_online"))
                set_allowed = stack.enter_context(mock.patch.object(
                    session, "_set_allowed"))
                effective = stack.enter_context(mock.patch.object(
                    session, "_systemctl_effective",
                    side_effect=("0 64", "0 64")))
                allowed = stack.enter_context(mock.patch.object(
                    session, "_systemctl_allowed", return_value=""))
                stack.enter_context(mock.patch.object(
                    session, "_restart_fillers", side_effect=restart_fillers))
                with self.assertRaisesRegex(
                        subject.CampaignError,
                        "interrupted during restoration by SIGTERM"):
                    session.close()
            self.assertEqual(read_text.call_count, 2)
            write_online.assert_called_once_with(64, "1")
            self.assertEqual(set_allowed.call_args_list, [
                mock.call("system.slice", "0,64"),
                mock.call("system.slice", ""),
            ])
            self.assertEqual(effective.call_count, 2)
            allowed.assert_called_once_with("system.slice")
            self.assertFalse(session.active)
            self.assertFalse(session.signals_installed)
            self.assertEqual(session.sibling_states, {})
            self.assertEqual(session.slice_states, {})
            for signum, handler in previous.items():
                self.assertEqual(signal.getsignal(signum), handler)
        finally:
            if session.signals_installed:
                session.pending_cleanup_signals.clear()
                session._restore_signal_handlers()

    def test_interrupted_turbostat_communication_is_reaped(self) -> None:
        probe = subject.LinuxTimingEvidenceProbe(
            mock.Mock(), mock.Mock(), Path("/thermal"), Path("/evidence"))
        for communicate_effect in (
                [KeyboardInterrupt()],
                [subprocess.TimeoutExpired("turbostat", 10),
                 KeyboardInterrupt()]):
            with self.subTest(calls=len(communicate_effect)):
                process = mock.Mock()
                process.pid = 12345
                process.poll.return_value = None
                process.communicate.side_effect = communicate_effect
                token = subject.TimingEvidenceToken(
                    {}, process, {}, Path("/thermal-out"),
                    Path("/performance-out"), ("turbostat",))
                with ExitStack() as stack:
                    stack.enter_context(mock.patch.object(subject.os, "killpg"))
                    reap = stack.enter_context(mock.patch.object(
                        subject.common, "stop_and_reap_process_group"))
                    with self.assertRaises(KeyboardInterrupt):
                        probe.finish(token)
                reap.assert_called_once_with(process)

    def test_interrupted_turbostat_startup_is_reaped(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-turbostat-startup-") as temporary:
            host = mock.Mock()
            host.isolation_is_live.return_value = True
            host.core = 0
            host.isolation.sibling_cpus = ()
            host.tool.side_effect = lambda name: f"/frozen/{name}"
            probe = subject.LinuxTimingEvidenceProbe(
                mock.Mock(), host, Path("/thermal"), Path(temporary))
            process = mock.Mock()
            process.pid = 12345
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "thermal_start", return_value={}))
                stack.enter_context(mock.patch.object(
                    subject.subprocess, "Popen", return_value=process))
                stack.enter_context(mock.patch.object(
                    subject.time, "sleep", side_effect=KeyboardInterrupt))
                reap = stack.enter_context(mock.patch.object(
                    subject.common, "stop_and_reap_process_group"))
                with self.assertRaises(KeyboardInterrupt):
                    probe.start(mock.Mock(), 1, None)
            reap.assert_called_once_with(process)

    def test_state_root_restart_and_manifest_mutation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-state-test-") as temporary:
            parent = Path(temporary)
            result = parent / "campaign"
            result.mkdir()
            self.required_h1_tree(result)
            repo, commit = self.fake_repo(parent)
            contract = self.fake_contract(repo, commit)
            prepare = {
                "prepared_utc_ns":
                    (subject.now_utc_ms() - 1000) * 1_000_000,
            }
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "verify_frozen_controller_runtime",
                    return_value=(prepare, contract),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "github_remote_url",
                    return_value="git@github.com:fake/wirehair.git",
                ))
                stack.enter_context(mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=os.environ.copy(),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "obtain_verified_latest_bound",
                    return_value=self.latest_bound(),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "verify_phase_ready_for_seal",
                    return_value=None,
                ))
                stack.enter_context(redirect_stdout(io.StringIO()))
                self.assertEqual(subject.seal_holdout(argparse.Namespace(
                    result_dir=result, phase="h1")), 0)
                state = subject.load_fixed_json(
                    result / "beacon/h1/state.json", "test state")
                _attempt, _manifest, seal, _seal_sha = subject.load_seal_attempt(
                    result, subject.HOLDOUT_PHASES["h1"], state)
                beacon = dict(subject.DRAND_KNOWN_ROUND)
                beacon["round"] = seal["round"]
                origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
                wave = {
                    "schema": "wirehair.wh2.drand_quicknet_wave.v1",
                    "status": "QUORUM",
                    "chain_hash": subject.DRAND_CHAIN_HASH,
                    "consensus_origins": origins[:3],
                    "observations": [
                        {"origin": origin, "ok": True} for origin in origins
                    ],
                    "beacon": beacon,
                }
                real_run = subject.subprocess.run

                def run(command: object, *args: object, **kwargs: object):
                    if (isinstance(command, tuple) and len(command) >= 4 and
                            command[-1] == str(seal["round"]) and
                            "wh2_drand_verify.cjs" in str(command[1])):
                        return subprocess.CompletedProcess(
                            command, 0, json.dumps(wave).encode("ascii"), b"")
                    return real_run(command, *args, **kwargs)

                def offline(
                    _contract: object, _result: object, value: dict[str, object],
                ) -> dict[str, object]:
                    return {
                        "beacon": value,
                        "canonical_sha256": hashlib.sha256(
                            subject.canonical_beacon_bytes(value)).hexdigest(),
                    }

                acquire = argparse.Namespace(
                    result_dir=result, phase="h1", no_wait=False,
                    max_wait_seconds=1.0,
                )
                with ExitStack() as acquire_stack:
                    acquire_stack.enter_context(mock.patch.object(
                        subject.subprocess, "run", side_effect=run))
                    acquire_stack.enter_context(mock.patch.object(
                        subject, "now_utc_ms",
                        return_value=seal["round_time_ms"] + 1,
                    ))
                    acquire_stack.enter_context(mock.patch.object(
                        subject, "offline_verify_beacon", side_effect=offline),
                    )
                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    root = subject.load_fixed_json(
                        result / "holdout/h1_roots.json", "test root")
                    self.assertEqual(root["status"], "ROOTED")
                    self.assertEqual(len(root["roots"]), subject.H1_ROOT_COUNT)
                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    self.assertEqual(len(list(
                        (result / "beacon/h1/waves").glob("wave-*.json"))), 1)
                    root_path = result / "holdout/h1_roots.json"
                    state_path = result / "beacon/h1/state.json"
                    root_bytes = root_path.read_bytes()
                    state_bytes = state_path.read_bytes()
                    for field, value in (
                            ("manifest_sha256", "f" * 64),
                            ("round", root["round"] + 1),
                            ("round_time_ms", root["round_time_ms"] + 1)):
                        mutated_root = dict(root)
                        mutated_root[field] = value
                        subject.atomic_state_json(root_path, mutated_root)
                        with self.assertRaises(subject.CampaignError):
                            subject.acquire_holdout(acquire)
                        subject.common.atomic_write(root_path, root_bytes)
                    rooted_state = subject.load_fixed_json(
                        state_path, "test rooted state")
                    mutated_state = dict(rooted_state)
                    mutated_state["manifest_sha256"] = "d" * 64
                    subject.atomic_state_json(state_path, mutated_state)
                    with self.assertRaises(subject.CampaignError):
                        subject.seal_holdout(argparse.Namespace(
                            result_dir=result, phase="h1"))
                    subject.common.atomic_write(state_path, state_bytes)
                    for field, value in (
                            ("root_file", "holdout/wrong.json"),
                            ("root_file_sha256", "e" * 64)):
                        mutated_state = dict(rooted_state)
                        mutated_state[field] = value
                        subject.atomic_state_json(state_path, mutated_state)
                        with self.assertRaises(subject.CampaignError):
                            subject.acquire_holdout(acquire)
                        with self.assertRaises(subject.CampaignError):
                            subject.seal_holdout(argparse.Namespace(
                                result_dir=result, phase="h2"))
                        subject.common.atomic_write(state_path, state_bytes)
                    (result / "development/frozen_table.json").write_text(
                        "mutated\n", encoding="ascii")
                    with self.assertRaises(subject.CampaignError):
                        subject.acquire_holdout(acquire)


if __name__ == "__main__":
    unittest.main()
