"""Neutral ELF mutation tests; no codec executable or scientific mode is run."""
import copy
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest

import Wh2OrderedRhsReleaseIsolation as I


class IsolationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = I.QUALIFIED/'native/CMakeFiles/release_candidate.dir/WirehairV2SolveOrdered.cpp.o'
        cls.before = I.Elf(cls.source.read_bytes())
        cls.mapping = I.mapping_for([cls.before])
        cls.temporary = tempfile.TemporaryDirectory(prefix='wh2-isolation-test-')
        directory = Path(cls.temporary.name)
        cls.mapfile = directory/'symbols.map'
        cls.mapfile.write_text(''.join(a+' '+b+'\n' for a, b in cls.mapping.items()))
        cls.destination = directory/'mapped.o'
        subprocess.check_call(['/usr/bin/objcopy', '--redefine-syms='+str(cls.mapfile),
                               str(cls.source), str(cls.destination)], timeout=10)
        cls.after = I.Elf(cls.destination.read_bytes())

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_exact_mapping_and_shared_groups(self):
        result = I.audit(self.before, self.after, self.mapping)
        self.assertGreater(result['mapped'], 0)
        self.assertGreater(I.audit_shared_groups(self.before, self.mapping), 0)
        signatures = [self.before.symbols[s[7]] for s in self.before.sections if s[1] == 17]
        self.assertTrue(any(s[0] in self.mapping and s[1] >> 4 == 0 for s in signatures))

    def test_local_signature_omission_rejected(self):
        partial = {s[0]: self.mapping[s[0]] for s in self.before.symbols
                   if s[0] in self.mapping and s[1] >> 4 != 0}
        self.assertLess(len(partial), len(self.mapping))
        with self.assertRaisesRegex(ValueError, 'symbol metadata'):
            I.audit(self.before, self.after, partial)

    def test_code_byte_change_rejected(self):
        section = next(s for s in self.after.sections if s[2] & 4 and s[5] > 0)
        data = bytearray(self.after.data)
        data[section[4]] ^= 1
        with self.assertRaisesRegex(ValueError, 'code/data/debug'):
            I.audit(self.before, I.Elf(bytes(data)), self.mapping)

    def test_relocation_change_rejected(self):
        section = next(s for s in self.after.sections if s[1] == 4 and s[5] > 0)
        data = bytearray(self.after.data)
        data[section[4]+8] ^= 1
        with self.assertRaisesRegex(ValueError, 'relocation semantics'):
            I.audit(self.before, I.Elf(bytes(data)), self.mapping)

    def test_symbol_binding_change_rejected(self):
        obj = copy.deepcopy(self.after)
        i = next(i for i, s in enumerate(obj.symbols) if s[0].startswith(I.PREFIX) and s[1] >> 4 == 1)
        s = obj.symbols[i]
        obj.symbols[i] = (s[0], (s[1] & 15) | 32) + s[2:]
        with self.assertRaisesRegex(ValueError, 'symbol metadata'):
            I.audit(self.before, obj, self.mapping)

    def test_bss_size_change_rejected(self):
        obj = copy.deepcopy(self.after)
        i = next(i for i, s in enumerate(obj.sections) if s[1] == 8 and s[5] > 0)
        s = obj.sections[i]
        obj.sections[i] = s[:5] + (s[5]+8,) + s[6:]
        with self.assertRaisesRegex(ValueError, 'section size'):
            I.audit(self.before, obj, self.mapping)

    def test_group_signature_cross_binding_rejected(self):
        obj = copy.deepcopy(self.after)
        i = next(i for i, s in enumerate(obj.sections) if s[1] == 17 and
                 obj.symbols[s[7]][0].startswith(I.PREFIX))
        replacement = next(s[7] for s in obj.sections if s[1] == 17 and
                           not obj.symbols[s[7]][0].startswith(I.PREFIX))
        s = obj.sections[i]
        obj.sections[i] = s[:7] + (replacement,) + s[8:]
        with self.assertRaisesRegex(ValueError, 'COMDAT identity'):
            I.audit(self.before, obj, self.mapping)

    def test_private_reference_from_shared_group_rejected(self):
        obj = copy.deepcopy(self.before)
        name = next(s[0] for s in obj.symbols if s[3] == 0 and s[0] == '_ZdlPv')
        mapping = dict(self.mapping, **{name: I.PREFIX+name})
        with self.assertRaisesRegex(ValueError, 'shared COMDAT references V2'):
            I.audit_shared_groups(obj, mapping)

    def test_writable_literal_rejected(self):
        obj = copy.deepcopy(self.before)
        i = next(i for i, s in enumerate(obj.sections) if s[2] == 0x32)
        s = obj.sections[i]
        obj.sections[i] = s[:2] + (s[2] | 1,) + s[3:]
        with self.assertRaisesRegex(ValueError, 'private allocated section'):
            I.audit_shared_groups(obj, self.mapping)

    def test_sanitizer_groups_isolated(self):
        path = I.QUALIFIED/'asan/CMakeFiles/release_candidate.dir/WirehairV2SolveOrdered.cpp.o'
        obj = I.Elf(path.read_bytes())
        mapping = I.mapping_for([obj], sanitizer=True)
        self.assertEqual(I.audit_shared_groups(obj, mapping), 0)
        self.assertGreater(len(mapping), len(I.mapping_for([obj])))

    def test_bad_elf_rejected(self):
        for data in (b'', b'\0'*64, self.before.data[:63]):
            with self.assertRaises(ValueError):
                I.Elf(data)
        data = bytearray(self.before.data)
        struct.pack_into('<Q', data, 40, len(data)+1)
        with self.assertRaisesRegex(ValueError, 'section table'):
            I.Elf(bytes(data))

    def test_double_transform_rejected(self):
        with self.assertRaisesRegex(ValueError, 'already transformed'):
            I.mapping_for([self.after])


if __name__ == '__main__':
    unittest.main()
