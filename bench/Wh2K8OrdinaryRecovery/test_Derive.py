"""Closed derivation and installed-parity mutations; never launch codec science."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import Derive as D
import Retained as I


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class DeriveTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = D.templates()
        cls.output = D.derive(cls.inputs)
        cls.temporary = tempfile.TemporaryDirectory(prefix='wh2-ordinary-recovery-test.')
        cls.directory = D.prepare(Path(cls.temporary.name)/'generated')
        cls.controller = module('ordinary_recovery_test', cls.directory/'Wh2K8OrdinaryRecoveryR0.py')

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_templates_and_anchor_rejection(self):
        for name in self.inputs:
            bad = dict(self.inputs); bad[name] += '\n'
            with self.assertRaises(ValueError): D.derive(bad)
        for text in ('missing', 'anchor anchor'):
            with self.assertRaises(ValueError): D.replace(text, 'anchor', 'changed')
        with self.assertRaises(ValueError): D.derive({})

    def test_exact_constructor_routes_and_unchanged_worker(self):
        old = self.inputs['Wh2K8PublicRecoveryR0.cpp']
        new = self.output['Wh2K8OrdinaryRecoveryR0.cpp']
        control = new[new.index('int PCreate('):new.index('int PEncode(')]
        candidate = new[new.index('int SCreate('):new.index('int SIndependent(')]
        self.assertIn('wirehair_v2_encoder_create_profile_id(', control)
        self.assertIn('wirehair_v2_encoder_create_profile_id_with_options(', control)
        self.assertEqual(control.count('WIREHAIR_V2_PROFILE_CERTIFIED_2026_07'), 2)
        self.assertNotIn('wirehair_v2_encoder_create(', control)
        self.assertIn('wirehair_v2_encoder_create(source', candidate)
        self.assertIn('wirehair_v2_encoder_create_with_options(source', candidate)
        self.assertNotIn('encoder_create_profile', candidate)
        start = 'const Api apis[6]'
        self.assertEqual(old[old.index(start):].replace('ordinary selects certified profile',
                         'explicit control selects certified profile'), new[new.index(start):])

    def test_unchanged_roster_math_decisions_and_prototype_parity(self):
        names = ('retained', 'roster', 'products', 'powers', 'coefficient', 'message',
                 'packet_hash', 'first_success', 'candidate_first', 'summarize',
                 'retained_parity', 'capture')
        def bodies(text):
            return {node.name: ast.dump(node) for node in ast.parse(text).body
                    if isinstance(node, ast.FunctionDef) and node.name in names}
        self.assertEqual(bodies(self.inputs['Wh2K8PublicRecoveryR0.py']),
                         bodies(self.output['Wh2K8OrdinaryRecoveryR0.py']))
        self.assertEqual(len(bodies(self.output['Wh2K8OrdinaryRecoveryR0.py'])), len(names))

    def test_fresh_imports_scope_and_bound_installed_evidence(self):
        C = self.controller
        self.assertEqual(C.ROOT, D.ROOT)
        self.assertEqual(C.PROTOCOL, D.PROTOCOL)
        self.assertEqual(C.U.PROTOCOL, D.PROTOCOL)
        self.assertEqual(Path(C.U.__file__).parent, self.directory)
        self.assertEqual(C.U.G.__file__, str(D.HERE/'Wh2K8OrdinaryRecovery/Derive.py'))
        self.assertEqual(C.OUTPUT, Path('/var/tmp/wh2-k8-ordinary-recovery-r0'))
        self.assertEqual(C.PAIRED, ((1,0),(1,2),(4,3),(4,5)))
        self.assertEqual(C.ARMS, ('certified_independent','ordinary_k8_independent','wh1_owned',
                                 'certified_borrowed','ordinary_k8_borrowed','wh1_borrowed'))
        self.assertEqual(C.SOURCES, D.SOURCES)
        self.assertEqual(C.I.PINS, I.PINS)
        self.assertEqual(len(C.roster()), 6260)
        self.assertEqual(sum(len(e['ids']) for e in C.roster()), 74946)
        builder = self.output['Wh2K8OrdinaryRecoveryBuildR0.py']
        self.assertIn('G.verify(GENERATED)', builder)
        self.assertNotIn('C.verify_derivation', builder)
        self.assertIn('dependencies.update(reader.I.PINS)', builder)
        self.assertIn('reader.I.neutral(fixture_records, mode, A)', builder)
        self.assertIn('archives = [candidate["object"], original_archive]', builder)
        self.assertIn("'encoder_create_profile_id', 'encoder_create_profile_id_with_options'", builder)

    def test_generation_tamper_and_fresh_directory_checks(self):
        with tempfile.TemporaryDirectory(prefix='wh2-ordinary-recovery-tamper.') as tmp:
            directory = D.prepare(Path(tmp)/'generated')
            self.assertEqual(D.verify(directory), directory)
            with self.assertRaises(ValueError): D.prepare(directory)
            with self.assertRaises(ValueError): D.prepare(D.ROOT/'not-created')
            for name, raw in D.encoded_outputs().items():
                path = directory/name; path.write_bytes(raw+b'\n')
                with self.subTest(name=name), self.assertRaises(ValueError): D.verify(directory)
                path.write_bytes(raw)
            receipt = json.loads((directory/'DERIVATION.json').read_bytes())
            receipt['outputs'].pop('Wh2K8OrdinaryRecoveryR0.cpp')
            (directory/'DERIVATION.json').write_text(json.dumps(receipt))
            with self.assertRaises(ValueError): D.verify(directory)

    def test_qualified_selector_objects_and_recipes(self):
        C = self.controller.U.C
        C.source_identity((D.ROOT/'codec/WirehairV2Profile.cpp').read_bytes(),
                          (C.QUALIFIED/C.SOURCE_NAME).read_bytes())
        database = json.loads((C.QUALIFIED/'compile_commands.json').read_bytes())
        for mode in C.OBJECTS:
            row, flags = C.recipe(database, mode)
            self.assertIn('-O3', flags)
            self.assertEqual(D.sha((C.QUALIFIED/row['output']).read_bytes()), C.OBJECTS[mode][1])

    def test_actual_installed_neutral_evidence(self):
        A = self.controller.A
        reference = None
        for mode in ('native', 'scalar', 'asan'):
            raw = I.read(I.QUALIFIED/mode/'fixtures.jsonl', A)
            records = I.records(raw, 48, mode, 'neutral', A)
            I.neutral(records, mode, A)
            if reference is None: reference = records
            else: self.assertEqual(records, reference)
        with self.assertRaises(ValueError): I.neutral([], 'unknown', A)

    def test_installed_parity_exact_and_nonmutating(self):
        A = self.controller.A
        path = I.BUNDLE/'native.raw.jsonl'
        record = dict(type='record', group=0, index=0, B=2, tail=2, ids=list(range(12)),
                      arms=[dict(profile='profile', packets=['packet'], rows='rows', feed=[1,0],
                                 first=8, recoveries=2, counts=[9,108,1,8,2,10], checked=True)
                            for _ in range(6)])
        expected = [record]*6260
        header = dict(type='header', protocol=I.PROTOCOL, backend='native', scope='retained')
        footer = dict(type='footer', records=6260, checked=True)
        raw = b''.join(A.canonical(r) for r in [header]+expected+[footer])
        pin = dict(path=str(path), bytes=len(raw), sha256=A.sha(raw))
        complete = A.canonical(dict(protocol=I.PROTOCOL, outcome='PASS', files=[pin]))
        pins = {path:A.sha(raw), I.BUNDLE/'COMPLETE.json':A.sha(complete)}
        with patch.object(I, 'PINS', pins), patch.object(A, 'read_regular',
                side_effect=lambda p, cap:raw if p==path else complete):
            before = copy.deepcopy(expected)
            I.verify(expected, A, lambda p:pin)
            self.assertEqual(expected, before)
            changes = [lambda r:r.pop(), lambda r:r[0].__setitem__('index',1),
                       lambda r:r[0]['ids'].reverse()]
            for arm in range(6):
                for field in record['arms'][arm]:
                    changes.append(lambda r, a=arm, f=field:r[0]['arms'][a].__setitem__(f,None))
            for change in changes:
                changed = list(expected); changed[0] = copy.deepcopy(record)
                change(changed)
                with self.assertRaisesRegex(ValueError, 'parity'): I.verify(changed, A, lambda p:pin)
            with self.assertRaisesRegex(ValueError, 'COMPLETE binding'):
                I.verify(expected, A, lambda p:dict(pin,bytes=1))
            for altered_path in pins:
                with patch.dict(I.PINS, {altered_path:'0'*64}):
                    with self.assertRaisesRegex(ValueError, 'evidence bytes'):
                        I.verify(expected, A, lambda p:pin)


if __name__ == '__main__':
    unittest.main()
