"""Read-only/synthetic current-library closure checks; no codec workload."""
import copy
import json
from pathlib import Path
import unittest

import Wh2K8PublicRecoveryBuildR0 as B


class BuildTest(unittest.TestCase):
    def test_sanitizer_startup_is_frozen_before_link(self):
        self.assertEqual(len(B.LINK_INPUTS),len(set(B.LINK_INPUTS)))
        for name in ('libasan_preinit.o','libasan.so','libubsan.so'):
            self.assertIn(name,B.LINK_INPUTS)
            path=Path(B.command(['c++','-print-file-name='+name]).decode().strip())
            self.assertTrue(path.is_absolute() and path.is_file())

    def test_all_current_production_recipes(self):
        for mode in B.MODES:
            database=json.loads((B.PRODUCTION/mode/'compile_commands.json').read_bytes())
            recipes=B.producer_recipes(database,mode)
            self.assertEqual(len(recipes),19)
            self.assertEqual([str(s.relative_to(B.ROOT)) for s,_,_ in recipes],list(B.PRODUCERS))
            self.assertTrue(all('WIREHAIR_TESTING' not in ' '.join(flags) for _,_,flags in recipes))
            for change in ('drop','duplicate','flag','directory','source','output'):
                altered=copy.deepcopy(database)
                prefix='CMakeFiles/wirehair_objects.dir/' if mode=='native' else 'CMakeFiles/wirehair.dir/'
                index=next(i for i,r in enumerate(altered) if r['output'].startswith(prefix))
                if change=='drop':altered.pop(index)
                elif change=='duplicate':altered.append(copy.deepcopy(altered[index]))
                elif change=='flag':altered[index]['command']+=' -DWIREHAIR_TESTING=1'
                elif change=='directory':altered[index]['directory']='/different/build'
                elif change=='source':altered[index]['file']='/different/source.cpp'
                else:altered[index]['output']+='other'
                with self.assertRaises(ValueError):B.producer_recipes(altered,mode)

    def test_pin_identity_and_conflicting_duplicate(self):
        record=dict(path='/synthetic/file',bytes=12,sha256='a'*64)
        self.assertEqual(B.pin_map([record,record]),{'/synthetic/file':record})
        other=dict(record,sha256='b'*64)
        with self.assertRaises(ValueError):B.pin_map([record,other])
        for changed in (dict(record,path='relative'),dict(record,path='/a/../file'),
                        dict(record,bytes=True),dict(record,sha256='g'*64)):
            with self.assertRaises(ValueError):B.pin_map([changed])

    def test_single_private_context(self):
        row=b'00001000 00022810 B GF256Ctx\n'
        self.assertEqual(B.context_size(row),141328)
        for raw in (b'',row+row,row.replace(b'B GF',b'T GF')):
            with self.assertRaises(ValueError):B.context_size(raw)

    def test_dependency_target_and_resolution(self):
        target=Path('/synthetic/observer.o')
        source=Path(__file__).resolve()
        raw=(str(target)+': '+str(source)+' \\\n '+str(source)+'\n').encode()
        self.assertEqual(B.preprocessor_dependencies(raw,target),{source})
        for altered in (raw.replace(b'observer.o:',b'wrong.o:'),b'/synthetic/observer.o: relative.h\n'):
            with self.assertRaises(ValueError):B.preprocessor_dependencies(altered,target)


if __name__=='__main__':unittest.main()
