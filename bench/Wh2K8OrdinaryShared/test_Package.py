"""Synthetic package evidence tests; no codec or scientific workload."""
import copy
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import Package as K


class PackageTest(unittest.TestCase):
    def suite(self, arm):
        names = K.names(arm) | {'package_plugin_round_trip'}
        suite = ET.Element('testsuite', tests=str(len(names)), failures='0', disabled='0')
        for name in sorted(names): ET.SubElement(suite, 'testcase', name=name, status='run')
        return suite

    def test_exact_rosters(self):
        for arm, count in (('baseline', 9), ('candidate', 10)):
            self.assertEqual(len(self.suite(arm)), count)
            K.check_tests(ET.tostring(self.suite(arm)), arm)
        with self.assertRaises(ValueError): K.names('unknown')

    def test_incomplete_or_failed_results_reject(self):
        for arm in ('baseline', 'candidate'):
            for mutation in ('omit', 'duplicate', 'fail', 'skip', 'disabled', 'error', 'notrun'):
                suite = self.suite(arm)
                if mutation == 'omit': suite.remove(suite[0])
                elif mutation == 'duplicate': suite[0].set('name', suite[1].attrib['name'])
                elif mutation == 'fail': ET.SubElement(suite[0], 'failure')
                elif mutation == 'skip': ET.SubElement(suite[0], 'skipped')
                elif mutation == 'disabled': suite.set('disabled', '1')
                elif mutation == 'error': ET.SubElement(suite[0], 'error')
                else: suite[0].set('status', 'notrun')
                with self.assertRaises(ValueError): K.check_tests(ET.tostring(suite), arm)
        with self.assertRaises(ValueError): K.check_tests(ET.tostring(self.suite('baseline')), 'candidate')

    def test_only_dso_bytes_may_change(self):
        original = {str(K.LIBRARY): dict(kind='file', bytes=10, sha256='a'*64),
                    'include/api.h': dict(kind='file', bytes=5, sha256='b'*64),
                    'lib.so': dict(kind='symlink', target='lib.so.2')}
        K.check_overlay(original, copy.deepcopy(original))
        candidate = copy.deepcopy(original)
        candidate[str(K.LIBRARY)]['sha256'] = 'c'*64
        candidate[str(K.LIBRARY)]['bytes'] = 11
        candidate_pin = dict(bytes=11, sha256='c'*64)
        K.check_overlay(original, candidate, candidate_pin)
        for mutation in ('header', 'link', 'extra', 'missing', 'library'):
            bad = copy.deepcopy(candidate)
            if mutation == 'header': bad['include/api.h']['sha256'] = 'd'*64
            elif mutation == 'link': bad['lib.so']['target'] = 'elsewhere'
            elif mutation == 'extra': bad['extra'] = bad['include/api.h']
            elif mutation == 'missing': del bad['include/api.h']
            else: bad[str(K.LIBRARY)]['sha256'] = 'd'*64
            with self.assertRaises(ValueError): K.check_overlay(original, bad, candidate_pin)

    def test_relative_links_and_escape_rejection(self):
        with tempfile.TemporaryDirectory(prefix='wh2-package-inventory-test.') as directory:
            parent = Path(directory)
            root = parent/'prefix'
            root.mkdir()
            (root/'library').write_bytes(b'content')
            (root/'alias').symlink_to('library')
            result = K.inventory(root)
            self.assertEqual(result['alias'], dict(kind='symlink', target='library'))
            (parent/'outside').write_bytes(b'outside')
            (root/'escape').symlink_to('../outside')
            with self.assertRaises(ValueError): K.inventory(root)

    def test_test_helpers_are_not_codec_symbols(self):
        K.check_defined_symbols(b'000 T wirehair_package_round_trip\n'
                                b'001 t wirehair_package_v2_round_trip.constprop.0\n', 'package_c_consumer')
        K.check_defined_symbols(b'000 T wirehair_plugin_round_trip\n', 'package_plugin')
        for target in ('package_plugin', 'package_c_consumer', 'package_ordinary_cpp'):
            for name in ('wirehair_v2_decode', 'wirehair_init_', 'gf256_mul_mem', 'GF256Ctx'):
                with self.assertRaises(ValueError): K.check_defined_symbols(('000 T '+name+'\n').encode(), target)
        with self.assertRaises(ValueError):
            K.check_defined_symbols(b'000 T wirehair_package_round_trip\n', 'package_ordinary_cpp')
        with self.assertRaises(ValueError):
            K.check_defined_symbols(b'000 T wirehair_plugin_round_trip\n', 'package_c_consumer')

    def test_symlinked_directory_rejects(self):
        with tempfile.TemporaryDirectory(prefix='wh2-package-directory-test.') as directory:
            root = Path(directory)
            (root/'data').mkdir()
            (root/'data/file').write_bytes(b'content')
            (root/'alias').symlink_to('data', target_is_directory=True)
            with self.assertRaises(ValueError): K.inventory(root)


if __name__ == '__main__': unittest.main()
