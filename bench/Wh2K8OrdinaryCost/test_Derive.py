"""Derivation/selector/linkage checks; no scientific worker is launched."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import Candidate as C
import Derive as D


class DeriveTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = D.templates()
        cls.output = D.derive(cls.inputs)

    def test_exact_closed_templates_and_anchors(self):
        for name in self.inputs:
            bad = dict(self.inputs)
            bad[name] += '\n'
            with self.assertRaises(ValueError):
                D.derive(bad)
        for source in ('missing', 'anchor anchor'):
            with self.assertRaises(ValueError):
                D.replace(source, 'anchor', 'changed')
        self.assertEqual(D.replace('anchor', 'anchor', 'changed'), 'changed')

    def test_constructor_routes_and_unchanged_work(self):
        old = self.inputs['Wh2K8PublicCostR0.cpp']
        new = self.output['Wh2K8OrdinaryCostR0.cpp']
        control = new[new.index('int PCreate('):new.index('int PEncode(')]
        candidate = new[new.index('int SCreate('):new.index('int SIndependent(')]
        self.assertIn('wirehair_v2_encoder_create_profile_id(', control)
        self.assertIn('wirehair_v2_encoder_create_profile_id_with_options(', control)
        self.assertEqual(control.count('WIREHAIR_V2_PROFILE_CERTIFIED_2026_07'), 2)
        self.assertNotIn('wirehair_v2_encoder_create(', control)
        self.assertIn('wirehair_v2_encoder_create(source', candidate)
        self.assertIn('wirehair_v2_encoder_create_with_options(source', candidate)
        self.assertNotIn('encoder_create_profile', candidate)
        for start, end in (('const Api apis[6]', 'void Initialize()'),
                           ('int Worker(', 'struct FakeReader')):
            self.assertEqual(old[old.index(start):old.index(end)], new[new.index(start):new.index(end)])
        old = self.inputs['Wh2K8PublicCostR0.py']
        new = self.output['Wh2K8OrdinaryCostR0.py']
        self.assertEqual(old[old.index('def roster('):old.index('def clocks(')],
                         new[new.index('def roster('):new.index('def clocks(')])

    def test_generated_binding_and_tamper_rejection(self):
        with tempfile.TemporaryDirectory(prefix='wh2-ordinary-derive-test.') as tmp:
            directory = D.prepare(Path(tmp)/'generated')
            self.assertEqual(D.verify(directory), directory)
            with self.assertRaises(ValueError):
                D.prepare(directory)
            for name, raw in D.encoded_outputs().items():
                path = directory/name
                path.write_bytes(raw+b'\n')
                with self.subTest(name=name), self.assertRaises(ValueError):
                    D.verify(directory)
                path.write_bytes(raw)
            receipt = json.loads((directory/'DERIVATION.json').read_bytes())
            receipt['outputs'].pop('Wh2K8OrdinaryCostR0.cpp')
            (directory/'DERIVATION.json').write_text(json.dumps(receipt))
            with self.assertRaises(ValueError):
                D.verify(directory)

    def test_fresh_import_routes(self):
        with tempfile.TemporaryDirectory(prefix='wh2-ordinary-import-test.') as tmp:
            directory = D.prepare(Path(tmp)/'generated')
            spec = importlib.util.spec_from_file_location('ordinary_cost_under_test', directory/'Wh2K8OrdinaryCostR0.py')
            controller = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(controller)
            self.assertEqual(controller.ROOT, D.ROOT)
            self.assertEqual(controller.PROTOCOL, D.PROTOCOL)
            self.assertEqual(controller.R.PROTOCOL, D.PROTOCOL)
            self.assertEqual(Path(controller.R.__file__).parent, directory)
            self.assertEqual(controller.R.OBSERVERS[0], directory/'Wh2K8OrdinaryCostR0.cpp')
            self.assertEqual(controller.OUTPUT, Path('/var/tmp/wh2-k8-ordinary-cost-r0'))
            self.assertEqual(controller.CALLBACKS, 77760)
            self.assertEqual(controller.PAIRS, ((0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,1),(2,1),(3,4),(5,4)))
            self.assertEqual(controller.ARM_NAMES, ('certified_independent','ordinary_k8_independent','wh1_owned',
                                                   'certified_borrowed','ordinary_k8_borrowed','wh1_borrowed'))
            assertions = [node for node in ast.walk(ast.parse(self.output['test_Wh2K8OrdinaryCostR0.py']))
                          if isinstance(node, ast.Call) and len(node.args) == 2 and
                          isinstance(node.args[0], ast.Attribute) and node.args[0].attr == 'ARM_NAMES']
            self.assertEqual(len(assertions), 1)
            self.assertEqual(ast.literal_eval(assertions[0].args[1]), controller.ARM_NAMES)
            self.assertEqual(C.verify_derivation(directory), {directory/n for n in D.encoded_outputs()})

    def test_selected_source_and_real_recipes(self):
        C.source_identity((D.ROOT/'codec/WirehairV2Profile.cpp').read_bytes(),
                          (C.QUALIFIED/C.SOURCE_NAME).read_bytes())
        with self.assertRaises(ValueError):
            C.source_identity(b'wrong', b'wrong')
        database = json.loads((C.QUALIFIED/'compile_commands.json').read_bytes())
        for mode in C.OBJECTS:
            row, flags = C.recipe(database, mode)
            self.assertIn('-O3', flags)
            self.assertIn(C.OBJECTS[mode][0], row['output'])
            for change in ('drop', 'duplicate', 'hook', 'directory', 'source', 'backend'):
                bad = copy.deepcopy(database)
                index = next(i for i, x in enumerate(bad) if x['output'] == row['output'])
                if change == 'drop':
                    bad.pop(index)
                elif change == 'duplicate':
                    bad.append(copy.deepcopy(bad[index]))
                elif change == 'hook':
                    bad[index]['command'] = bad[index]['command'].replace(' -o ', ' -DWIREHAIR_TESTING=1 -o ')
                elif change == 'directory':
                    bad[index]['directory'] = '/wrong'
                elif change == 'source':
                    bad[index]['file'] = '/wrong'
                else:
                    bad[index]['command'] = bad[index]['command'].replace(' -o ', ' -DANDROID=1 -o ') if mode != 'scalar' else bad[index]['command'].replace('-DANDROID=1 ', '')
                with self.subTest(mode=mode, change=change), self.assertRaises(ValueError):
                    C.recipe(bad, mode)

    def test_exact_member_extraction(self):
        archive = Path('/synthetic/libwirehair.a')
        candidate = Path('/synthetic/candidate.o')
        members = ['wirehair.cpp.o', 'gf256.cpp.o', 'WirehairSmall.cpp.o', 'WirehairSmallK8.cpp.o',
                   'WirehairV2Profile.cpp.o']
        load = 'LOAD '+str(candidate)+'\nLOAD '+str(archive)+'\n'
        selected = ''.join(str(archive)+'('+name+')\n' for name in members[:-1])
        self.assertEqual(C.check_link(load+selected, archive, candidate, members), sorted(members[:-1]))
        for bad in (selected, load+selected+str(archive)+'(WirehairV2Profile.cpp.o)\n',
                    load+selected+str(archive)+'(unknown.o)\n', load+selected.replace('gf256.cpp.o','other.o'),
                    load+selected+'LOAD '+str(candidate)+'\n',
                    'LOAD '+str(archive)+'\nLOAD '+str(candidate)+'\n'+selected):
            with self.assertRaises(ValueError):
                C.check_link(bad, archive, candidate, members)


if __name__ == '__main__':
    unittest.main()
