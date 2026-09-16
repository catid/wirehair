import copy
import io
import unittest
import Analyze as A


def observation():
    common = dict.fromkeys(A.HEADER, 0)
    common.update(k=3, width=2, tail=2, dso_base=0x800000)
    result = []
    pointers = (0x1000, 0x2000, 0x3000, 0x3000, 0x2000, 0x1000)
    for i in range(192):
        slot = i % 6
        r = dict(common)
        r.update(event=i, allocate=int(slot < 3), array=int(slot in (2, 3)),
                 bytes=(296, 104, 8, 0, 0, 0)[slot], pointer=pointers[slot],
                 caller_owner=int(slot == 5), caller_offset=0x100+slot, stack=0x4000+slot,
                 phase=0 if slot < 3 else 3)
        result.append(r)
    for i in range(192):
        cycle, slot = divmod(i, 6)
        phase = 0 if slot == 0 else 1 if slot <= 3 else 2 if slot == 4 else 3
        r = dict(common)
        first = 0x5000 if phase == 0 else 0x8000+(slot-1)*2 if phase == 1 else 0x9000+cycle*18 if phase == 2 else 0
        r.update(kind=1, event=i, phase=phase, pointer=0 if phase == 0 else 0x1000,
                 bytes=(32, 2, 2, 2, 6, 0)[slot], arg1=first,
                 arg2=0x6000 if phase in (0, 2) else 0, packet=slot+2 if phase == 1 else 0)
        result.append(r)
    return result


class Test(unittest.TestCase):
    def test_complete_observation(self):
        e, c, feeds = A.validate(observation(), (0, 0, 0, 0, 0))
        self.assertEqual((len(e), len(c), feeds), (192, 192, 3))

    def test_allocations_and_frees(self):
        for row, column, value in ((0, 'bytes', 295), (1, 'pointer', 0x1000),
                                   (3, 'pointer', 0x2000), (3, 'phase', 1),
                                   (0, 'caller_owner', 1), (5, 'array', 1)):
            rows = observation()
            rows[row][column] = value
            with self.assertRaises(ValueError):
                A.validate(rows, (0, 0, 0, 0, 0))

    def test_arguments(self):
        for row, column, value in ((192, 'arg1', 0x1000), (192, 'arg2', 0x2000),
                                   (193, 'packet', 4), (193, 'pointer', 0x3000),
                                   (196, 'bytes', 5), (197, 'arg1', 0x5000),
                                   (196, 'arg1', 2**64-3), (192, 'arg2', 0x5000),
                                   (196, 'arg2', 0x9000), (194, 'arg1', 0x8003),
                                   (202, 'arg1', 0x9000), (198, 'arg1', 0x7000)):
            rows = observation()
            rows[row][column] = value
            with self.assertRaises(ValueError):
                A.validate(rows, (0, 0, 0, 0, 0))

    def test_roster_and_extent(self):
        for rows in ([], observation()[:192], observation()[:-1], observation()*4):
            with self.assertRaises(ValueError):
                A.validate(rows, (0, 0, 0, 0, 0))
        with self.assertRaises(ValueError):
            A.validate(observation(), (0, 1, 0, 0, 0))

    def test_comparisons(self):
        a = A.validate(observation(), (0, 0, 0, 0, 0))
        b = copy.deepcopy(a)
        self.assertTrue(all(A.comparison(a, b).values()))
        b[0][0]['pointer'] += 1
        self.assertEqual(A.comparison(a, b)['same_private_addresses'], 0)
        b[0][0]['stack'] += 1
        self.assertEqual(A.comparison(a, b)['same_allocator_stack'], 0)
        b[1][0]['arg2'] += 1
        self.assertEqual(A.comparison(a, b)['same_public_arguments'], 0)

    def test_csv(self):
        header = ','.join(A.HEADER)+'\n'
        row = ','.join(str(observation()[0][n]) for n in A.HEADER)+'\n'
        self.assertEqual(list(A.read_rows([header, row])), observation()[:1])
        for bad in (row.rstrip(), row.replace('296', '-1'), row.replace('296', str(2**64)), 'x'*1025+'\n'):
            with self.assertRaises(ValueError):
                list(A.read_rows([header, bad]))
        with self.assertRaises(ValueError):
            list(A.read_rows(['wrong header\n']))

    def test_incomplete_campaign(self):
        with self.assertRaises(ValueError):
            A.analyze(io.StringIO(','.join(A.HEADER)+'\n'))

    def test_shapes(self):
        self.assertEqual(A.shape(0), (3, 2, 2, 0))
        self.assertEqual(A.shape(41), (8, 1280, 1, 1))
        for index in (-1, 42):
            with self.assertRaises(ValueError):
                A.shape(index)


if __name__ == '__main__':
    unittest.main()
