#!/usr/bin/env python3
"""Synthetic portability-adapter tests; never run a native codec campaign."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2ProductionMix3RecoveryR0 as original
import Wh2ProductionMix3RecoveryBroadR0 as broad

ORIGINAL_STATE = (original.PROTOCOL, original.HOLDOUT_ROOTS, original.source_bytes,
                  original.profile_bytes, original.trace, original.validate_cell)
BROAD_STATE = (broad.c.PROTOCOL, broad.c.HOLDOUT_ROOTS, broad.c.source_bytes,
               broad.c.profile_bytes, broad.c.trace, broad.c.validate_cell)
import Wh2ProductionMix3RecoveryPortabilityR0 as p

c = p.c
fixture_spec = importlib.util.spec_from_file_location(
    "_wh2_portability_fixture_private", Path(__file__).with_name("test_Wh2ProductionMix3RecoveryR0.py"))
fixtures = importlib.util.module_from_spec(fixture_spec)
fixture_spec.loader.exec_module(fixtures)
fixtures.c = c
fixtures.LIBRARY_PATH = str(broad.PINNED_LIBRARY)
_original_cell = fixtures.cell


def width_cell(*args):
    row = _original_cell(*args)
    if row["packets_hex"] is not None:
        K, attempt = row["K"], row["attempt"]
        width = struct.unpack_from("<I", bytes.fromhex(row["profile_hex"]), 24)[0]
        source = c.source_bytes(K)
        row["packets_hex"] = b"".join(
            source[packet_id * width:(packet_id + 1) * width] if packet_id < K else
            bytes((K + attempt + packet_id + i) & 255 for i in range(width))
            for packet_id in row["ids"]).hex()
    return row


fixtures.cell = width_cell


def records(failures=None, selected=(5, 3), baseline=(0, 0)):
    return fixtures.transcript("holdout", selected, baseline, failures)


def parsed(rows):
    return fixtures.parse(rows, "holdout", [5, 3])


def synthetic_receipt():
    sources = {path: c.digest(path.encode("ascii")) for path in p.SOURCE_PATHS}
    basis = {"source_commit": broad.PRODUCTION_COMMIT, "library": str(broad.PINNED_LIBRARY),
             "library_sha256": broad.LIBRARY_SHA256,
             "source_files": {"bench/Wh2ProductionMix3RecoveryR0.py": sources["bench/Wh2ProductionMix3RecoveryR0.py"]}}
    return {"protocol": p.PROTOCOL, "source_commit": fixtures.SOURCE_COMMIT, "source_files": sources,
            "production_build": basis, "production_source_commit": broad.PRODUCTION_COMMIT,
            "library": str(broad.PINNED_LIBRARY), "library_sha256": broad.LIBRARY_SHA256,
            "worker": str(p.WORKER_SOURCE), "worker_sha256": "2" * 64,
            "interpreter": str(Path(sys.executable).resolve()), "interpreter_sha256": "3" * 64,
            "interpreter_version": sys.version, "interpreter_flags": ["-I", "-B", "-S"],
            "map_basis": {"protocol": broad.PROTOCOL, "selected_attempts": [5, 3], "train_sha256": "4" * 64}}


def synthetic_bundle(receipt):
    phases, files = [], {"freeze.json": p.frozen_protocol(receipt), "selection.json": receipt["map_basis"]}
    for width in p.WIDTHS:
        p.configure_width(width)
        rows = records()
        phases.append(parsed(rows))
        files["b%d.jsonl" % width] = fixtures.encoded(rows).encode("ascii")
        files["b%d.stderr" % width] = b""
    summary = p.summarize(phases)
    summary.update({"protocol": p.PROTOCOL, "elapsed_seconds": 0.25,
                    "promotion_claimed": False, "all_K_claimed": False, "timing_claimed": False})
    files["summary.json"] = summary
    return files


class PortabilityAdapterTest(unittest.TestCase):
    def setUp(self):
        p.configure_width(2)

    def tearDown(self):
        p.configure_width(2)

    def test_private_modules_and_exact_fresh_roots(self):
        self.assertIsNot(p.b, broad)
        self.assertIsNot(c, broad.c)
        self.assertIsNot(c, original)
        roots = tuple("0x" + hashlib.sha256((
            "wirehair.wh2.production-mix3-k3k6-portability-r0:holdout/" + str(i))
            .encode("ascii")).hexdigest()[:16] for i in range(32))
        self.assertEqual(c.HOLDOUT_ROOTS, roots)
        self.assertEqual(c.digest(c.canonical(roots)),
                         "a739277bae151083192f2cf9c9c56fe5c321467a639775b4821ba4f03293c53d")
        self.assertEqual(len(set(roots)), 32)
        reserved = {"0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3",
                    "0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b",
                    "0xefd20c982041a46b", "0x8827bc36ed906555", "0x86029f23d6132efa"}
        self.assertFalse(set(roots) & (reserved | set(original.TRAINING_ROOTS + original.HOLDOUT_ROOTS)
                                      | set(broad.TRAINING_ROOTS + broad.HOLDOUT_ROOTS)))
        for width in (64, 1280, 2):
            p.configure_width(width)
        self.assertEqual((original.PROTOCOL, original.HOLDOUT_ROOTS, original.source_bytes,
                          original.profile_bytes, original.trace, original.validate_cell), ORIGINAL_STATE)
        self.assertEqual((broad.c.PROTOCOL, broad.c.HOLDOUT_ROOTS, broad.c.source_bytes,
                          broad.c.profile_bytes, broad.c.trace, broad.c.validate_cell), BROAD_STATE)

    def test_width_profiles_sources_and_trace_seed_equivalence(self):
        for width in (2, 64, 1280):
            p.configure_width(width)
            for K in (3, 6):
                self.assertEqual(c.source_bytes(K), bytes((73 * i + 19 * K + 11) & 255
                                                         for i in range(K * width)))
                for attempt in (0, 3, 5, 255):
                    self.assertEqual(c.profile_bytes(K, attempt), b"WHV2" + struct.pack(
                        "<HHQQIB3s", 1, 32, 0x4b295bbb47f4f9c9, K * width, width, attempt, bytes(3)))
                for root in (c.HOLDOUT_ROOTS[0], c.HOLDOUT_ROOTS[-1], "0xffffffffffffffff"):
                    adjusted = (int(root, 16) ^ (width * 0xbf58476d1ce4e5b9)
                                ^ (2 * 0xbf58476d1ce4e5b9)) & ((1 << 64) - 1)
                    for schedule in c.SCHEDULES:
                        self.assertEqual(c.trace(K, root, schedule),
                                         original.trace(K, hex(adjusted), schedule))

    def test_every_width_accepts_complete_fixed_map_and_rejects_reselection(self):
        for width in (2, 64, 1280):
            p.configure_width(width)
            rows = records()
            result = parsed(rows)
            self.assertEqual(len(result["cells"]), 384)
            self.assertEqual(result["selected_attempts"], [5, 3])
            self.assertEqual(c.paired_counts(result, [5, 3])["cells_per_arm"], 192)
            for changed in (rows[:-2] + rows[-1:], records(selected=(5, 4))):
                with self.subTest(width=width), self.assertRaises(c.ValidationError):
                    parsed(changed)
            # A failed fixed candidate remains the same map, never starts training.
            failed = parsed(records(failures=lambda phase, arm, K, attempt, ri, schedule:
                                    (1, 0, 0, 0) if arm == "candidate" and ri == 0 and
                                    K == 6 and schedule == "burst" else None))
            self.assertEqual(failed["selected_attempts"], [5, 3])
            self.assertEqual(c.paired_counts(failed, [5, 3])["introductions"], [1, 0, 0, 0])

    def test_wide_packet_receipts_and_width_identity_are_not_two_byte_checks(self):
        for width in (64, 1280):
            p.configure_width(width)
            rows = records()
            row_index = next(i for i, row in enumerate(rows) if row.get("type") == "cell"
                             and any(packet_id < row["K"] for packet_id in row["ids"]))
            target = rows[row_index]
            for mutation in ("last_systematic_byte", "two_byte_packet", "wrong_width", "trace"):
                changed = copy.deepcopy(rows)
                row = changed[row_index]
                if mutation == "last_systematic_byte":
                    packet_index = next(i for i, packet_id in enumerate(row["ids"]) if packet_id < row["K"])
                    payload = bytearray.fromhex(row["packets_hex"])
                    payload[(packet_index + 1) * width - 1] ^= 1
                    row["packets_hex"] = payload.hex()
                elif mutation == "two_byte_packet":
                    row["packets_hex"] = row["packets_hex"][:4 * (row["K"] + 4)]
                elif mutation == "wrong_width":
                    row["profile_hex"] = original.profile_bytes(row["K"], row["attempt"]).hex()
                    row["profile_sha256"] = c.digest(bytes.fromhex(row["profile_hex"]))
                else:
                    row["ids"][0] ^= 1
                    row["trace_sha256"] = c.digest(b"".join(struct.pack("<I", x) for x in row["ids"]))
                with self.subTest(width=width, mutation=mutation), self.assertRaises(c.ValidationError):
                    parsed(changed)
            self.assertEqual(len(bytes.fromhex(target["packets_hex"])), (target["K"] + 4) * width)

    def test_native_signatures_and_encoder_lifetime_with_mocked_library(self):
        ct, lib = p.ct, mock.Mock()
        lib.wirehair_init_.return_value = 0
        handle_value = 0x123456789

        def assign(pointer, kind, value):
            ct.cast(pointer, ct.POINTER(kind))[0] = value

        def baseline(message, size, width, profile, capacity, written, handle):
            self.assertEqual((message.raw, size, width, capacity), (p.source_bytes(3, 64), 192, 64, 32))
            profile.raw = p.profile_bytes(3, 64, 7)
            assign(written, ct.c_uint32, 32)
            assign(handle, ct.c_void_p, handle_value)
            return 0

        def deserialize(profile, size, decoded):
            self.assertEqual(size, 32)
            ct.cast(decoded, ct.POINTER(p.Profile)).contents.seed_attempt = profile.raw[28]
            return 0

        def serialize(decoded, output, capacity, written):
            actual = ct.cast(decoded, ct.POINTER(p.Profile)).contents.seed_attempt
            output.raw = p.profile_bytes(3, 64, actual)
            assign(written, ct.c_uint32, 32)
            return 0

        def encode(handle, packet_id, output, capacity, written):
            self.assertEqual((handle.value, capacity), (handle_value, 64))
            output.raw = p.source_bytes(3, 64)[packet_id * 64:(packet_id + 1) * 64]
            assign(written, ct.c_uint32, 64)
            return 0

        lib.wirehair_v2_encoder_create.side_effect = baseline
        lib.wirehair_v2_profile_deserialize.side_effect = deserialize
        lib.wirehair_v2_profile_serialize.side_effect = serialize
        lib.wirehair_v2_encode.side_effect = encode
        with mock.patch.object(c, "file_digest", return_value=broad.LIBRARY_SHA256), \
                mock.patch.object(ct, "CDLL", return_value=lib) as loader:
            native = p.Native()
        loader.assert_called_once_with(str(broad.PINNED_LIBRARY))
        lib.wirehair_init_.assert_called_once_with(2)
        u32, u64, ptr = ct.c_uint32, ct.c_uint64, ct.c_void_p
        expected = {
            "wirehair_init_": [ct.c_int],
            "wirehair_v2_encoder_create": [ptr, u64, u32, ptr, u32, ct.POINTER(u32), ct.POINTER(ptr)],
            "wirehair_v2_encoder_create_profile": [ptr, ptr, u32, ct.POINTER(ptr)],
            "wirehair_v2_profile_deserialize": [ptr, u32, ct.POINTER(p.Profile)],
            "wirehair_v2_profile_serialize": [ct.POINTER(p.Profile), ptr, u32, ct.POINTER(u32)],
            "wirehair_v2_decoder_create": [ptr, u32, ct.POINTER(ptr)],
            "wirehair_v2_encode": [ptr, u32, ptr, u32, ct.POINTER(u32)],
            "wirehair_v2_decode": [ptr, u32, ptr, u32],
            "wirehair_v2_recover": [ptr, ptr, u64, ct.POINTER(u64)], "wirehair_v2_free": [ptr]}
        for name, arguments in expected.items():
            self.assertEqual(getattr(lib, name).argtypes, arguments)
            self.assertIs(getattr(lib, name).restype, None if name.endswith("free") else ct.c_int)
        with self.assertRaisesRegex(RuntimeError, "consumer failure"):
            with native.encoder(3, 64, None) as (handle, profile, message, source, bad, actual):
                self.assertEqual((handle.value, actual, bad), (handle_value, 7, False))
                self.assertEqual((profile, message.raw), (p.profile_bytes(3, 64, 7), source))
                raise RuntimeError("consumer failure")
        self.assertEqual(lib.wirehair_v2_free.call_args.args[0].value, handle_value)
        lib.wirehair_v2_encoder_create_profile.return_value = 10
        with native.encoder(3, 64, 5) as (handle, profile, message, source, bad, actual):
            self.assertEqual((handle.value, actual, bad), (None, 5, True))
        self.assertIsNone(lib.wirehair_v2_free.call_args.args[0].value)
        self.assertEqual(lib.wirehair_v2_encode.call_count, 3)  # BadSeed never encodes.

    def test_native_decode_offsets_byte_oracle_and_cleanup_are_mocked(self):
        ct, native = p.ct, object.__new__(p.Native)
        native.lib = lib = mock.Mock()
        source, packets = p.source_bytes(3, 64), b"A" * 64 + b"B" * 64 + b"C" * 64

        def create(profile, size, handle):
            ct.cast(handle, ct.POINTER(ct.c_void_p))[0] = 0x123456789
            return 0

        def decode(handle, packet_id, data, size):
            self.assertEqual((handle.value, size), (0x123456789, 64))
            expected = b"A" if packet_id == 0xffffffff else b"B"
            self.assertEqual(ct.string_at(data, size), expected * 64)
            return 1 if packet_id == 0xffffffff else 0

        def recover(handle, output, capacity, written):
            self.assertEqual(capacity, len(source))
            output.raw = source
            ct.cast(written, ct.POINTER(ct.c_uint64))[0] = len(source)
            return 0

        lib.wirehair_v2_decoder_create.side_effect = create
        lib.wirehair_v2_decode.side_effect = decode
        lib.wirehair_v2_recover.side_effect = recover
        args = (p.profile_bytes(3, 64, 5), source, 64, [0xffffffff, 7, 9], packets, 3)
        self.assertEqual(native.decode(*args), 0)
        self.assertEqual(lib.wirehair_v2_decode.call_count, 2)  # Stops after success.
        lib.wirehair_v2_recover.side_effect = lambda *unused: 0
        with self.assertRaisesRegex(c.ValidationError, "recovered bytes"):
            native.decode(*args)
        self.assertEqual(lib.wirehair_v2_free.call_count, 2)
        self.assertTrue(all(call.args[0].value == 0x123456789 for call in lib.wirehair_v2_free.call_args_list))

    def test_receipt_pins_and_coherent_inherited_source_drift(self):
        receipt = synthetic_receipt()
        live = {receipt[key]: receipt[key + "_sha256"] for key in ("worker", "library", "interpreter")}
        map_hash = c.digest(c.canonical(receipt["map_basis"]) + b"\n")
        live[str(p.MAP_PATH)] = map_hash
        with mock.patch.object(p.b, "BASIS_BUILD_SHA256", c.digest(c.canonical(receipt["production_build"]))), \
                mock.patch.object(p, "MAP_SHA256", map_hash), \
                mock.patch.object(c, "command", return_value=fixtures.SOURCE_COMMIT), \
                mock.patch.object(p, "source_receipt", return_value=receipt["source_files"]), \
                mock.patch.object(c, "file_digest", side_effect=lambda path: live[str(path)]):
            p.check_build(receipt)
            for path in live:
                saved, live[path] = live[path], "f" * 64
                with self.subTest(artifact=path), self.assertRaises(c.ValidationError):
                    p.check_build(receipt)
                live[path] = saved
            forged = copy.deepcopy(receipt)
            path = "bench/Wh2ProductionMix3RecoveryR0.py"
            forged["source_files"][path] = "f" * 64
            with mock.patch.object(p, "source_receipt", return_value=forged["source_files"]), \
                    self.assertRaisesRegex(c.ValidationError, "inherited source drift"):
                p.check_build(forged)
            for key in ("production_source_commit", "map_basis"):
                forged = copy.deepcopy(receipt)
                forged[key] = "f" * 40 if key == "production_source_commit" else {"selected_attempts": [5, 4]}
                with self.subTest(identity=key), self.assertRaises(c.ValidationError):
                    p.frozen_protocol(forged)

    def test_offline_replay_rehashed_mutations_and_publication_deadline(self):
        receipt = synthetic_receipt()
        with mock.patch.object(p.b, "BASIS_BUILD_SHA256", c.digest(c.canonical(receipt["production_build"]))), \
                mock.patch.object(p, "MAP_SHA256", c.digest(c.canonical(receipt["map_basis"]) + b"\n")), \
                mock.patch.object(c, "command", side_effect=AssertionError("live replay tool")), \
                mock.patch.object(p, "Native", side_effect=AssertionError("native replay")), \
                tempfile.TemporaryDirectory() as directory:
            bundle, files = Path(directory), synthetic_bundle(receipt)
            fixtures.publish_bundle(bundle, files)
            self.assertEqual(p.replay(bundle)["outcome"], "PASS")
            self.assertEqual(files["freeze.json"]["max_prefix_decodes"], 4608)
            self.assertIs(files["freeze.json"]["training"], False)
            for mutation in ("width_order", "map", "worker_stderr", "summary", "deadline", "scope"):
                changed = copy.deepcopy(files)
                if mutation == "width_order":
                    changed["freeze.json"]["widths"] = [2, 1280, 64]
                elif mutation == "map":
                    changed["selection.json"]["selected_attempts"] = [5, 4]
                elif mutation == "worker_stderr":
                    changed["b64.stderr"] = b"unexpected diagnostic\n"
                elif mutation == "summary":
                    changed["summary.json"]["widths"][1]["paired"]["candidate_failures"][0] = 1
                elif mutation == "deadline":
                    changed["summary.json"]["elapsed_seconds"] = 70
                else:
                    changed["summary.json"]["known_prior_counterexample_resolved"] = True
                fixtures.publish_bundle(bundle, changed)
                with self.subTest(mutation=mutation), self.assertRaises(c.ValidationError):
                    p.replay(bundle)
            # Hashing can consume the last outer-budget time; only INVALID may seal.
            fixtures.publish_bundle(bundle, files)
            (bundle / "COMPLETE").unlink()
            (bundle / "summary.json").unlink()
            now, real_digest = [0.0], c.file_digest

            def delayed_digest(path):
                now[0] = 70.0
                return real_digest(path)

            with mock.patch.object(c.time, "monotonic", side_effect=lambda: now[0]), \
                    mock.patch.object(c, "file_digest", side_effect=delayed_digest):
                result = c.publish_result(bundle, files["summary.json"], 0.0)
            self.assertEqual(result["outcome"], "INVALID")
            self.assertEqual(p.replay(bundle), result)


if __name__ == "__main__":
    unittest.main()
