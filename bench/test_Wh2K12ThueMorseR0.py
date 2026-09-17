"""Independent synthetic checks for the finite K12 structural screen."""
import hashlib
import itertools
import tempfile
import unittest
from unittest import mock

import Wh2K12ThueMorseR0 as M


class Tests(unittest.TestCase):
    def test_field_and_rank_oracles(self):
        identity = M.identity()
        self.assertEqual(M.matrix_rank(identity), 12)
        self.assertEqual(M.matrix_rank(list(identity[:-1]) + [(0,) * 12]), 11)
        for index in range(256):
            raw = hashlib.sha512(('k12-neutral/' + str(index)).encode()).digest() * 3
            rows = [raw[offset:offset + 12] for offset in range(0, 144, 12)]
            self.assertEqual(M.matrix_rank(rows), 12)
        with self.assertRaises(ValueError): M.matrix_rank([[True] + [0] * 11])
        with self.assertRaises(ValueError): M.matrix_rank([[256] + [0] * 11])

    def test_feedback_companion_and_local_roster(self):
        feedback = M.fixed_feedback()
        self.assertEqual(len(feedback), 12)
        pair = (M.companion(feedback), M.companion((feedback[0] ^ 1,) + feedback[1:]))
        self.assertNotEqual(M.matrix_multiply(pair[0], pair[1]), M.matrix_multiply(pair[1], pair[0]))
        for word in M.WORDS:
            columns = M.local_columns(pair, word)
            self.assertEqual(len(columns), 16)
            self.assertEqual(M.matrix_rank(columns[:12]), 12)
        self.assertEqual(len(M.MINORS), 1820)
        self.assertEqual(len(set(M.MINORS)), 1820)

    def test_mapper_sequential_and_seams(self):
        feedback = M.fixed_feedback()
        pair = (M.companion(feedback), M.companion((feedback[0] ^ 1,) + feedback[1:]))
        mapper = M.Mapper(pair, M.Budget(20))
        self.assertEqual(len(mapper.payload), 135168)
        product = M.identity()
        for index in range(2050):
            self.assertEqual(mapper.row(index), tuple(row[0] for row in product))
            product = M.matrix_multiply(product, pair[M.parity(index)])
        for exponent in range(3, 32):
            for offset in (-1, 0, 1):
                packet_id = (1 << exponent) + offset
                self.assertEqual(mapper.row(packet_id), mapper.reference_row(packet_id))
        self.assertEqual(mapper.row(M.MAX_ID), mapper.reference_row(M.MAX_ID))
        for bad in (-1, 1.0, True, 1 << 32):
            with self.assertRaises(ValueError): mapper.row(bad)

    def test_trace_reconstruction_and_fresh_roots(self):
        root = '0x123456789abcdef0'
        def reference(width, schedule):
            state = (int(root, 16) ^ 12 * 0x9e3779b97f4a7c15 ^ width * 0xbf58476d1ce4e5b9) & M.MASK64
            if schedule != 'iid': state ^= 0x10fade
            threshold = (1 / 9 if schedule == 'burst' else .1 if schedule == 'iid' else .5) * 2**53
            result, burst = [], 0
            for candidate in range(65536):
                if schedule == 'burst' and burst:
                    burst -= 1
                    continue
                state = (state + 0x9e3779b97f4a7c15) & M.MASK64
                value = ((state ^ (state >> 30)) * 0xbf58476d1ce4e5b9) & M.MASK64
                value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & M.MASK64
                if ((value ^ (value >> 31)) >> 11) < threshold:
                    if schedule == 'burst': burst = 7
                    continue
                result.append(M.MAX_ID - 2 * candidate if schedule == 'adversarial' else
                              12 + candidate if schedule == 'repair-only' else candidate)
                if len(result) == 16: return result
            raise AssertionError('trace bound')
        for width, schedule in itertools.product(M.WIDTHS, M.SCHEDULES):
            self.assertEqual(M.trace(width, root, schedule), reference(width, schedule))
        roots = M.fresh_roots({'roots': {'0x' + M.digest(('old/' + str(index)).encode())[:16]
                                         for index in range(64)}})
        self.assertEqual(len(roots), 512)
        with mock.patch.object(M, 'digest', return_value='0' * 64):
            with self.assertRaisesRegex(ValueError, 'fresh root collision'):
                M.fresh_roots({'roots': set()})

    def test_fresh_gate_and_no_candidate_reselection(self):
        rows = [dict(b=width, schedule=schedule, root=str(index), ranks=[12] * 5)
                for width, schedule in itertools.product(M.WIDTHS, M.SCHEDULES)
                for index in range(8)]
        self.assertTrue(M.summarize_fresh(rows, 8)['fresh_pass'])
        rows[0]['ranks'] = [11] * 5
        result = M.summarize_fresh(rows, 8)
        self.assertFalse(result['fresh_pass'])
        self.assertEqual(result['cells'][0]['first_success'], [7, 0, 0, 0, 0, 1])
        with self.assertRaises(ValueError): M.summarize_fresh(rows[:-1], 8)

    def test_compact_controller_claim_authentication(self):
        receipt = {'protocol': M.PROTOCOL, 'sentinel': True}
        claim = {'protocol': M.PROTOCOL, 'receipt': receipt,
                 'receipt_sha256': hashlib.sha256((M.json.dumps(
                     receipt, sort_keys=True, separators=(',', ':'),
                     ensure_ascii=True, allow_nan=False)).encode('ascii')).hexdigest()}
        with tempfile.TemporaryDirectory() as directory:
            output = M.Path(directory)
            raw = (M.json.dumps(claim, sort_keys=True, separators=(',', ':'),
                                ensure_ascii=True, allow_nan=False)).encode('ascii')
            (output / 'CLAIM.json').write_bytes(raw)
            with mock.patch.object(M, 'OUTPUT', output):
                self.assertEqual(M.claimed_inputs(), M.digest(raw))
            (output / 'CLAIM.json').write_bytes(M.canonical(claim))
            with mock.patch.object(M, 'OUTPUT', output):
                with self.assertRaisesRegex(ValueError, 'claim identity'):
                    M.claimed_inputs()

    def test_screen_stops_before_fresh_on_history_failure(self):
        pair = (M.companion((1,) + (0,) * 11),) * 2
        class FakeMapper:
            payload = bytes(135168)
            def __init__(self, *unused): self.pair = pair
            def row(self, packet_id): return M.identity()[packet_id % 12][0:]
        history = dict(origins=[], prefixes=[], roots=[])
        with mock.patch.object(M, 'history_inputs', return_value=history), \
             mock.patch.object(M, 'fixed_feedback', return_value=(1,) + (0,) * 11), \
             mock.patch.object(M, 'choose_pair', return_value=(pair, [])), \
             mock.patch.object(M, 'Mapper', FakeMapper), \
             mock.patch.object(M, 'check_mapper', return_value=[]), \
             mock.patch.object(M, 'check_history', side_effect=ValueError('history')):
            with self.assertRaisesRegex(ValueError, 'history'): M.run_screen('0' * 64)


if __name__ == '__main__': unittest.main()
