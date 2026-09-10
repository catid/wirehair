"""Synthetic coverage and evidence-reader tests; no codec calls."""
import copy
import _ctypes
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import Qualify as Q


class QualificationTest(unittest.TestCase):
    def suite(self):
        suite = ET.Element('testsuite', tests='31', failures='0', disabled='0', skipped='0')
        for name in sorted(Q.test_names()):
            ET.SubElement(suite, 'testcase', name=name, status='run')
        return suite

    def test_complete_success(self):
        self.assertEqual(len(Q.consumer_names()), 30)
        self.assertEqual(len(Q.test_names()), 31)
        Q.check_test_results(ET.tostring(self.suite()))

    def test_missing_duplicate_failed_skipped_disabled_reject(self):
        original = self.suite()
        for change in ('missing', 'duplicate', 'failed', 'skipped', 'disabled', 'not-run', 'count'):
            suite = copy.deepcopy(original)
            if change == 'missing': suite.remove(suite[0])
            elif change == 'duplicate': suite[0].set('name', suite[1].attrib['name'])
            elif change == 'failed': ET.SubElement(suite[0], 'failure')
            elif change == 'skipped': ET.SubElement(suite[0], 'skipped')
            elif change == 'disabled': suite.set('disabled', '1')
            elif change == 'not-run': suite[0].set('status', 'notrun')
            else: suite.set('tests', '30')
            with self.assertRaises(ValueError): Q.check_test_results(ET.tostring(suite))

    def test_transitive_helpers_and_ctypes_are_included(self):
        paths = Q.helper_files()
        self.assertIn(Path(Q.P.A.__file__).resolve(), paths)
        self.assertIn(Path(Q.L.R.N.__file__).resolve(), paths)
        self.assertIn(Path(Q.L.R.N.A.__file__).resolve(), paths)
        native_path = getattr(_ctypes, '__file__', None)
        if native_path:
            self.assertIn(Path(native_path).resolve(), paths)
        else:
            # This Python 3.8 distribution links _ctypes into its separately
            # pinned interpreter rather than loading an extension file.
            self.assertEqual(_ctypes.__spec__.origin, 'built-in')

    def test_nonexecutable_shared_runtime_detected(self):
        with tempfile.TemporaryDirectory(prefix='wh2-shared-elf-test.') as directory:
            path = Path(directory)/'runtime.so'
            data = bytearray(64)
            data[:4] = b'\x7fELF'
            data[5] = 1
            data[16] = 3
            path.write_bytes(data)
            path.chmod(0o600)
            self.assertTrue(Q.runtime_elf(path))
            data[16] = 1
            path.write_bytes(data)
            self.assertFalse(Q.runtime_elf(path))


if __name__ == '__main__': unittest.main()
