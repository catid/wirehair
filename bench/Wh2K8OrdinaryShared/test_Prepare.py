"""Read-only and synthetic shared-build recipe checks; no codec workload."""
import copy
import json
from pathlib import Path
import shlex
import unittest

import Prepare as P


class PreparationTest(unittest.TestCase):
    def recipes(self):
        database = json.loads((P.CURRENT/'compile_commands.json').read_bytes())
        return P.B.producer_recipes(database, 'native')

    def test_original_shared_recipe(self):
        recipes = self.recipes()
        raw = P.B.command(['/usr/bin/ninja', '-C', P.CURRENT, '-t', 'commands', 'libwirehair.so.2.0.0'])
        recipe = P.link_recipe(raw, recipes)
        self.assertEqual(recipe[9:-1], [str(obj.relative_to(P.CURRENT)) for _, obj, _ in recipes])
        self.assertEqual(recipe[-1], '-lm')
        self.assertFalse(any('symbolic' in arg or 'semantic-interposition' in arg for arg in recipe))
        for before, after in ((b'-shared', b'-shared -Wl,-Bsymbolic'),
                              (b'-O3', b'-O2'),
                              (b'-soname,libwirehair.so.2', b'-soname,libcandidate.so.2'),
                              (b'-Wl,--version-script=', b'-Wl,--bad-script=')):
            with self.assertRaises(ValueError):
                P.link_recipe(raw.replace(before, after), recipes)
        with self.assertRaises(ValueError):
            P.link_recipe(b'\n'.join(raw.splitlines()[1:])+b'\n', recipes)
        reversed_rows = raw.splitlines()
        reversed_rows[0], reversed_rows[1] = reversed_rows[1], reversed_rows[0]
        with self.assertRaises(ValueError):
            P.link_recipe(b'\n'.join(reversed_rows)+b'\n', recipes)

    def test_replacement_preserves_object_position(self):
        objects = [Path('/fresh')/obj.name for _, obj, _ in self.recipes()]
        candidate = Path('/fresh/selector.o')
        selected = P.replacement_objects(objects, candidate)
        self.assertEqual(selected[P.PROFILE_INDEX], candidate)
        self.assertEqual(selected[:P.PROFILE_INDEX], objects[:P.PROFILE_INDEX])
        self.assertEqual(selected[P.PROFILE_INDEX+1:], objects[P.PROFILE_INDEX+1:])
        for bad in (objects[:-1], list(reversed(objects)), objects+objects[:1]):
            with self.assertRaises(ValueError):
                P.replacement_objects(bad, candidate)
        with self.assertRaises(ValueError):
            P.replacement_objects(objects, objects[P.PROFILE_INDEX])

    def test_selector_is_exactly_the_qualified_six_lines(self):
        original = (P.ROOT/'codec/WirehairV2Profile.cpp').read_bytes()
        candidate = (P.C.QUALIFIED/P.C.SOURCE_NAME).read_bytes()
        P.C.source_identity(original, candidate)
        with self.assertRaises(ValueError):
            P.C.source_identity(original, candidate+b'\n')
        database = json.loads((P.C.QUALIFIED/'compile_commands.json').read_bytes())
        row, flags = P.C.recipe(database, 'native')
        self.assertNotIn('-DWIREHAIR_TESTING=1', flags)
        altered = copy.deepcopy(database)
        next(r for r in altered if r['output'] == row['output'])['command'] += ' -DWIREHAIR_TESTING=1'
        with self.assertRaises(ValueError):
            P.C.recipe(altered, 'native')

    def test_link_inputs_cover_shared_startup(self):
        self.assertIn('crtbeginS.o', P.LINK_INPUTS)
        self.assertIn('crtendS.o', P.LINK_INPUTS)
        self.assertEqual(len(P.LINK_INPUTS), len(set(P.LINK_INPUTS)))
        self.assertEqual(P.PROFILE_INDEX, 16)
        self.assertEqual(shlex.split('-shared -Wl,-soname,libwirehair.so.2'),
                         ['-shared', '-Wl,-soname,libwirehair.so.2'])

    def test_complete_command_graph(self):
        recipes = self.recipes()
        raw = P.B.command(['/usr/bin/ninja', '-C', P.CURRENT, '-t', 'commands', 'libwirehair.so.2.0.0'])
        original = P.link_recipe(raw, recipes)
        commands = P.expected_commands(Path('/fresh'), recipes, original)
        self.assertEqual(len(commands), 62)
        self.assertEqual(commands[0]['argv'][0], '/usr/bin/ninja')
        self.assertEqual(commands[20]['cwd'], str(P.C.QUALIFIED))
        self.assertEqual(commands[59]['cwd'], str(P.C.QUALIFIED))
        for index in range(19):
            self.assertEqual(commands[22+2*index]['argv'][:2], ['/usr/bin/ar', 'p'])
        self.assertEqual(commands[60]['argv'][9+P.PROFILE_INDEX], '/fresh/WirehairV2Profile.cpp.o')
        self.assertEqual(commands[61]['argv'][9+P.PROFILE_INDEX], '/fresh/'+P.C.SOURCE_NAME+'.o')

    def test_dependencies_cannot_omit_headers_or_candidate_inputs(self):
        files = [Path(__file__).resolve(), P.ROOT/'include/wirehair/wirehair.h']
        declared = [P.B.pin(path) for path in files]
        pins = P.B.pin_map(declared)
        target = Path('/synthetic/producer.o')
        raw = (str(target)+': '+' '.join(map(str, files))+'\n').encode()
        P.check_dependencies(raw, target, declared, pins)
        for bad in ([], declared[:1], declared+declared[:1]):
            with self.assertRaises(ValueError): P.check_dependencies(raw, target, bad, pins)
        altered = copy.deepcopy(declared)
        altered[0]['sha256'] = '0'*64
        with self.assertRaises(ValueError): P.check_dependencies(raw, target, altered, pins)

    def test_shape_rejects_empty_or_redirected_proof(self):
        # In-memory synthetic upgrade of the retained first neutral proof. It
        # deliberately does not claim the initial proof passes today's verifier.
        output = Path('/tmp/wh2-k8-ordinary-shared-prep.NdDkOaWn/prepared')
        proof = json.loads((output/'PREPARED.json').read_bytes())
        proof['schema'] = 3
        proof['interpreter'] = next(row for row in proof['files'] if row['path'] == '/usr/bin/python3.12')
        proof['captures'] = [P.B.pin(output/('command-%03d.%s'%(i, suffix))) for i in range(62)
                             for suffix in ('request.json', 'stdout', 'stderr', 'result.json')]
        proof['captures'].extend(P.B.pin(path) for path in P.dependency_paths(output))
        proof['files'].append(P.B.pin(P.CURRENT/'CMakeFiles/rules.ninja'))
        P.proof_shape(proof, output)
        for key in ('files', 'commands', 'original_link', 'captures', 'members', 'links'):
            bad = copy.deepcopy(proof)
            bad[key] = []
            with self.assertRaises((ValueError, KeyError, IndexError)):
                P.proof_shape(bad, output)
        bad = copy.deepcopy(proof)
        bad['links'][1]['dso']['path'] = str(output/'baseline/libwirehair.so.2.0.0')
        with self.assertRaises(ValueError): P.proof_shape(bad, output)
        bad = copy.deepcopy(proof)
        bad['files'] = [row for row in bad['files'] if not row['path'].endswith('/cc1plus')]
        with self.assertRaises(ValueError): P.proof_shape(bad, output)
        bad = copy.deepcopy(proof)
        bad['captures'].append(bad['captures'][0])
        with self.assertRaises(ValueError): P.proof_shape(bad, output)


if __name__ == '__main__':
    unittest.main()
