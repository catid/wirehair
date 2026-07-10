#!/usr/bin/env python3

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

try:
    from . import apply_regeneration_results as apply
    from . import regeneration_results as results
except ImportError:
    import apply_regeneration_results as apply
    import regeneration_results as results


class ApplyRegenerationResultsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source_path = Path(__file__).resolve().parents[1] / "WirehairTools.cpp"
        cls.source = cls.source_path.read_text(encoding="ascii")
        cls.tables = apply.SourceTables(cls.source)

    def test_source_parser_and_fingerprints(self):
        self.assertEqual(len(self.tables.tiny_counts), 65)
        self.assertEqual(len(self.tables.small_dense), 1983)
        self.assertEqual(len(self.tables.dense_seeds), 100)
        self.assertEqual(len(self.tables.peel_seeds), 2048)
        self.assertGreater(self.tables.small_hash(), 0)
        self.assertGreater(self.tables.dense_graph_hash(), 0)
        self.assertGreater(self.tables.peel_hash(), 0)

    def test_numeric_and_dense_point_replacement(self):
        seeds = list(self.tables.dense_seeds)
        seeds[13] ^= 1
        rendered = apply.replace_numeric_array(self.source, "kDenseSeeds", seeds)
        reparsed = apply.SourceTables(rendered)
        self.assertEqual(reparsed.dense_seeds, seeds)
        points = [(2048, 54), (64000, 346)]
        rendered = apply.replace_dense_points(self.source, points)
        reparsed = apply.SourceTables(rendered)
        self.assertEqual(reparsed.dense_points, points)
        table_generator = (
            self.source_path.parent / "tables" / "TableGenerator.cpp"
        ).read_text(encoding="ascii")
        rendered_generator = apply.replace_dense_points(table_generator, points)
        self.assertEqual(apply.parse_dense_points(rendered_generator)[0], points)

    def test_identity_small_application_verifies_current_values(self):
        rows = []
        source_hash = self.tables.small_hash()
        for n in range(2, 2048):
            count = self.tables.dense_count(n)
            dense = self.tables.dense_seed(n, count)
            peel = self.tables.peel_seed(n)
            rows.append({
                "N": str(n), "DenseCount": str(count),
                "DenseSeed": str(dense), "PeelSeed": str(peel),
                "DenseFailures": "0", "PeelScore": "0",
                "CurrentDenseCount": str(count),
                "CurrentDenseSeed": str(dense),
                "CurrentPeelSeed": str(peel),
                "SmallTablesHash": str(source_hash),
            })
        rendered, ledger = apply.apply_small(self.source, rows)
        self.assertEqual(len(ledger), 2046)
        self.assertEqual(apply.SourceTables(rendered).small_hash(), source_hash)

    def test_dense_seed_usage_is_recomputed_from_graph(self):
        expected = {}
        for n in range(2048, 64001):
            expected[self.tables.dense_count(n) // 4] = n
        rows = []
        source_hash = self.tables.dense_seed_hash()
        for index in range(2, 100):
            rows.append({
                "DenseIndex": str(index),
                "Used": str(int(index in expected)),
                "N": str(expected.get(index, 0)),
                "DenseSeed": str(self.tables.dense_seeds[index]),
                "CurrentDenseSeed": str(self.tables.dense_seeds[index]),
                "Failures": "0",
                "SourceTablesHash": str(source_hash),
            })
        rendered, ledger = apply.apply_dense_seeds(self.source, rows)
        self.assertEqual(len(ledger), 98)
        self.assertEqual(apply.SourceTables(rendered).dense_seed_hash(), source_hash)
        rows[0]["Used"] = "1"
        rows[0]["N"] = "2048"
        with self.assertRaisesRegex(apply.ApplyError, "usage/N"):
            apply.apply_dense_seeds(self.source, rows)

    def test_identity_peel_application(self):
        source_hash = self.tables.peel_hash()
        rows = [
            {
                "Subdivision": str(index),
                "PeelSeed": str(seed),
                "CurrentPeelSeed": str(seed),
                "Status": "searched",
                "Score": "0",
                "BaseTablesHash": str(source_hash),
            }
            for index, seed in enumerate(self.tables.peel_seeds)
        ]
        rendered, ledger = apply.apply_peel(self.source, rows)
        self.assertEqual(len(ledger), 2048)
        self.assertEqual(apply.SourceTables(rendered).peel_hash(), source_hash)

    def test_dense_candidate_uses_max_and_normalizes(self):
        source_hash = self.tables.dense_count_hash()
        rows = []
        for seed in (11, 12, 13, 14):
            for n in results.canonical_dense_count_ns():
                rows.append({
                    "N": str(n),
                    "DenseCount": str(self.tables.dense_count(n)),
                    "BaseSeed": str(seed),
                    "SourceTablesHash": str(source_hash),
                })
        rendered, ledger = apply.apply_dense_count(self.source, rows)
        candidate = apply.SourceTables(rendered)
        self.assertEqual(candidate.dense_points[0][0], 2048)
        self.assertEqual(candidate.dense_points[-1][0], 64000)
        self.assertEqual(len(ledger), len(results.canonical_dense_count_ns()))

    def test_atomic_write_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "output"
            apply.write_atomic(path, "first")
            with self.assertRaisesRegex(apply.ApplyError, "already exists"):
                apply.write_atomic(path, "second")
            self.assertEqual(path.read_text(encoding="ascii"), "first")

    def test_multi_output_rolls_back_second_and_third_publish_failures(self):
        real_link = os.link
        for failure_index in (2, 3):
            with self.subTest(failure_index=failure_index), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                entries = [(root / name, name) for name in ("one", "two", "three")]
                calls = 0

                def injected_link(source, destination, **kwargs):
                    nonlocal calls
                    calls += 1
                    if calls == failure_index:
                        raise OSError("injected publish failure")
                    return real_link(source, destination, **kwargs)

                with mock.patch.object(
                    apply.atomic_publish.os, "link", side_effect=injected_link
                ):
                    with self.assertRaisesRegex(OSError, "injected"):
                        apply.publish_outputs(entries)
                for output, _ in entries:
                    self.assertFalse(output.exists())
                    self.assertFalse(Path(str(output) + ".tmp").exists())

    def test_concurrent_destination_is_preserved_and_other_outputs_rollback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = [(root / name, name) for name in ("one", "two", "three")]
            real_link = os.link
            calls = 0

            def concurrent_link(source, destination, **kwargs):
                nonlocal calls
                calls += 1
                if calls == 2:
                    Path(destination).write_text("concurrent", encoding="ascii")
                return real_link(source, destination, **kwargs)

            with mock.patch.object(
                apply.atomic_publish.os, "link", side_effect=concurrent_link
            ):
                with self.assertRaises(FileExistsError):
                    apply.publish_outputs(entries)
            self.assertFalse(entries[0][0].exists())
            self.assertEqual(entries[1][0].read_text(encoding="ascii"), "concurrent")
            self.assertFalse(entries[2][0].exists())
            for output, _ in entries:
                self.assertFalse(Path(str(output) + ".tmp").exists())

    def test_post_link_directory_sync_failure_rolls_back_every_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = [(root / name, name) for name in ("one", "two", "three")]
            real_fsync = os.fsync
            calls = 0

            def fail_directory_sync(fd):
                nonlocal calls
                calls += 1
                if calls == 4:
                    raise OSError("injected directory sync failure")
                return real_fsync(fd)

            with mock.patch.object(
                apply.atomic_publish.os, "fsync", side_effect=fail_directory_sync
            ):
                with self.assertRaisesRegex(OSError, "directory sync"):
                    apply.publish_outputs(entries)
            for output, _ in entries:
                self.assertFalse(output.exists())
                self.assertFalse(Path(str(output) + ".tmp").exists())

    def test_swapped_temporary_cannot_publish_attacker_content(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            real_link = os.link

            def swap_then_link(source, destination, **kwargs):
                Path(source).unlink()
                Path(source).write_text("attacker", encoding="ascii")
                return real_link(source, destination, **kwargs)

            with mock.patch.object(
                apply.atomic_publish.os, "link", side_effect=swap_then_link
            ):
                with self.assertRaisesRegex(
                    apply.ApplyError,
                    "does not match staged",
                ):
                    apply.publish_outputs([(output, "expected")])
            self.assertFalse(output.exists())
            self.assertEqual(
                Path(str(output) + ".tmp").read_text(encoding="ascii"),
                "attacker",
            )

    def test_post_link_destination_replacement_is_never_deleted(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            real_link = os.link

            def replace_destination(source, destination, **kwargs):
                real_link(source, destination, **kwargs)
                Path(destination).unlink()
                Path(destination).write_text("concurrent", encoding="ascii")

            with mock.patch.object(
                apply.atomic_publish.os, "link", side_effect=replace_destination
            ):
                with self.assertRaisesRegex(apply.ApplyError, "does not match staged"):
                    apply.publish_outputs([(output, "expected")])
            self.assertEqual(output.read_text(encoding="ascii"), "concurrent")
            self.assertFalse(Path(str(output) + ".tmp").exists())

    @unittest.skipIf(os.name == "nt", "symlink creation is privilege-dependent")
    def test_swapped_temporary_symlink_is_rejected_without_following(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "output"
            target = root / "target"
            target.write_text("attacker", encoding="ascii")
            real_link = os.link

            def symlink_then_link(source, destination, **kwargs):
                Path(source).unlink()
                Path(source).symlink_to(target)
                return real_link(source, destination, **kwargs)

            with mock.patch.object(
                apply.atomic_publish.os, "link", side_effect=symlink_then_link
            ):
                with self.assertRaisesRegex(apply.ApplyError, "does not match staged"):
                    apply.publish_outputs([(output, "expected")])
            self.assertFalse(output.exists())
            self.assertTrue(Path(str(output) + ".tmp").is_symlink())
            self.assertEqual(target.read_text(encoding="ascii"), "attacker")

    def test_directory_sync_can_be_disabled_on_platforms_without_support(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            with mock.patch.object(
                apply.atomic_publish, "DIRECTORY_FSYNC_SUPPORTED", False
            ), mock.patch.object(
                apply.atomic_publish.os, "open",
                side_effect=AssertionError("directory open must be skipped"),
            ):
                apply.publish_outputs([(output, "payload")])
            self.assertEqual(output.read_text(encoding="ascii"), "payload")

    def test_rollback_continues_after_one_unlink_is_denied(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = [(root / name, name) for name in ("one", "two", "three")]
            real_unlink = apply.atomic_publish._unlink_if_identity
            real_fsync = os.fsync
            fsync_calls = 0

            def fail_directory_sync(fd):
                nonlocal fsync_calls
                fsync_calls += 1
                if fsync_calls >= 4:
                    raise OSError("injected sync failure")
                return real_fsync(fd)

            def deny_one_output(path, identity):
                if Path(path).name == "three":
                    raise PermissionError("injected unlink denial")
                return real_unlink(path, identity)

            with mock.patch.object(
                apply.atomic_publish.os, "fsync", side_effect=fail_directory_sync
            ), mock.patch.object(
                apply.atomic_publish, "_unlink_if_identity",
                side_effect=deny_one_output,
            ):
                with self.assertRaisesRegex(
                    apply.atomic_publish.AtomicPublishError,
                    "rollback had",
                ):
                    apply.atomic_publish.publish_files_no_replace(entries)
            self.assertFalse(entries[0][0].exists())
            self.assertFalse(entries[1][0].exists())
            self.assertTrue(entries[2][0].exists())

    def test_unreadable_existing_final_and_temp_are_never_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for suffix in ("", ".tmp"):
                with self.subTest(suffix=suffix):
                    output = root / f"output-{len(suffix)}"
                    sentinel = Path(str(output) + suffix)
                    sentinel.write_text("sentinel", encoding="ascii")
                    sentinel.chmod(0)
                    try:
                        with self.assertRaisesRegex(apply.ApplyError, "already exists"):
                            apply.write_atomic(output, "replacement")
                    finally:
                        sentinel.chmod(0o600)
                    self.assertEqual(sentinel.read_text(encoding="ascii"), "sentinel")
                    if suffix:
                        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
