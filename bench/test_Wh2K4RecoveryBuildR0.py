"""Parser/provenance-helper tests only; no builds or codec programs."""
import copy
import os
from pathlib import Path
import shlex
import unittest
from unittest.mock import patch

import Wh2K4RecoveryBuildR0 as B


class BuildParsingTest(unittest.TestCase):
    def production_database(self, mode):
        prefix = 'CMakeFiles/'+('wirehair_objects' if mode == 'native' else 'wirehair')+'.dir/'
        flags = ['-DWIREHAIR_BUILDING=1']+([] if mode == 'native' else ['-DWIREHAIR_STATIC=1'])
        flags += ['-I'+str(B.ROOT/'include')]
        flags += dict(native=['-O3','-DNDEBUG'], scalar=['-DANDROID','-O3','-DNDEBUG'],
                      asan=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-g'])[mode]
        flags += ['-std=gnu++11','-fPIC','-Wall','-Wextra','-Wpedantic','-Werror']
        if mode == 'asan': flags.append('-march=native')
        result = []
        for name in B.PRODUCERS:
            source, output = str(B.ROOT/name), prefix+name+'.o'
            argv = ['/usr/bin/c++']+flags+['-o',output,'-c',source]
            result.append(dict(file=source, output=output, directory=str(B.PRODUCTION/mode),
                               command=' '.join(map(shlex.quote,argv))))
        return result

    def test_nineteen_production_recipes_are_data_not_executed(self):
        with patch.object(B, 'command', side_effect=AssertionError('must not execute recipe')):
            for mode in B.MODES:
                recipes = B.producer_recipes(self.production_database(mode), mode)
                self.assertEqual(len(recipes), 19)
                self.assertEqual([s for s,_,_ in recipes], [B.ROOT/name for name in B.PRODUCERS])
                if mode == 'asan':
                    self.assertTrue(all('-O1' not in flags and '-march=native' in flags for _,_,flags in recipes))
                for mutation in (lambda d:d.pop(), lambda d:d.append(copy.deepcopy(d[0])),
                                 lambda d:d[0].__setitem__('file',str(B.ROOT/'wrong.cpp')),
                                 lambda d:d[0].__setitem__('output',d[0]['output']+'other'),
                                 lambda d:d[0].__setitem__('directory','/different'),
                                 lambda d:d[0].__setitem__('command',d[0]['command']+' -O1')):
                    database = self.production_database(mode); mutation(database)
                    with self.assertRaises(ValueError): B.producer_recipes(database, mode)

    def test_rejects_changed_original_anchor_before_inspections(self):
        with patch.object(B, 'pin', return_value=dict(path=str(B.PRIOR/'COMPLETE.json'), bytes=1, sha256='0'*64)), \
                patch.object(B, 'command', side_effect=AssertionError('unauthenticated inspection')):
            with self.assertRaisesRegex(ValueError, 'original library COMPLETE'):
                B.qualified_inputs('native', Path('/tmp/synthetic-k4/native'))

    def test_read_only_actual_three_backend_provenance(self):
        # New reader of immutable data, not an invocation of an old verifier.
        # No output files, compiler or codec workloads are created here.
        for mode in B.MODES:
            archives, deps, proof, snapshots = B.qualified_inputs(mode, Path('/tmp/synthetic-k4')/mode)
            self.assertEqual(archives, [B.SMALL/mode/'k4/libwh2_small_serialized.a',
                                        B.PRODUCTION/mode/'libwirehair.a'])
            self.assertEqual(len(proof['historical_production']['members']),19)
            self.assertEqual(len(snapshots),1)
            original = proof['historical_snapshots'][0]
            self.assertEqual(original['original']['sha256'], B.CORE_OLD_SHA)
            self.assertEqual(B.A.sha(next(iter(snapshots.values()))), B.CORE_OLD_SHA)
            self.assertIn(B.ROOT/'codec/WirehairSmallCore.h', deps)
            self.assertNotIn('producing_source_closure',proof)

    def database(self, mode):
        build = B.SMALL/mode/'k4'
        target = 'CMakeFiles/wh2_small_serialized.dir/home/catid/wirehair/bench/Wh2SmallSerialized.cpp.o'
        flags = ['-DWH2_SMALL_CODEC_K=4']
        if mode == 'scalar': flags = ['-DANDROID']+flags+['-DWH2_SMALL_EXPECT_PORTABLE=1']
        flags += ['-I'+str(B.ROOT), '-I'+str(B.ROOT/'include'), '-I'+str(build)]
        flags += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-march=native', '-g'] if mode == 'asan' else ['-O3', '-DNDEBUG']
        flags += ['-std=c++11', '-fPIC', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-strict-aliasing', '-fno-lto']
        argv = ['/usr/bin/c++']+flags+['-o', target, '-c', str(B.HERE/'Wh2SmallSerialized.cpp')]
        entry = dict(directory=str(build), file=str(B.HERE/'Wh2SmallSerialized.cpp'), output=target,
                     command=' '.join(map(shlex.quote, argv)))
        return build, [entry]+[dict(output='CMakeFiles/test%d.dir/other.o'%i) for i in range(7)]

    def test_boundary_exact_flags_and_single_archive_member(self):
        for mode in B.MODES:
            build, database = self.database(mode)
            obj, argv = B.boundary_recipe(database, build, mode)
            self.assertEqual(obj.name, 'Wh2SmallSerialized.cpp.o')
            self.assertEqual('-march=native' in argv, mode == 'asan')
            self.assertEqual('-DANDROID' in argv, mode == 'scalar')
            self.assertIn('-fno-lto', argv)

    def test_boundary_rejects_duplicate_selector_and_backend_drift(self):
        build, original = self.database('native')
        for mutation in (lambda d:d.pop(), lambda d:d.__setitem__(1, copy.deepcopy(d[0])),
                         lambda d:d[0].__setitem__('file', str(B.HERE/'Wh2SmallNativeTest.cpp')),
                         lambda d:d[0].__setitem__('command', d[0]['command'].replace('-fno-lto', '-flto')),
                         lambda d:d[0].__setitem__('command', d[0]['command'].replace('=4', '=3', 1))):
            database = copy.deepcopy(original); mutation(database)
            with self.assertRaises(ValueError): B.boundary_recipe(database, build, 'native')

    def test_pin_overlap_only_identical(self):
        record = dict(path='/tmp/example', bytes=7, sha256='a'*64)
        self.assertEqual(B.pin_map([record, copy.deepcopy(record)]), {record['path']:record})
        for key, value in [('bytes', 8), ('sha256', 'b'*64)]:
            other = dict(record, **{key:value})
            with self.assertRaisesRegex(ValueError, 'conflicting duplicate pin'): B.pin_map([record, other])
        for bad in [dict(record, path='relative'), dict(record, path='/tmp/../bad'),
                    dict(record, bytes=True), dict(record, sha256='A'*64), dict(record, extra=0)]:
            with self.assertRaises(ValueError): B.pin_map([bad])

    def test_dependency_target_normalization_and_duplicates(self):
        source = B.HERE/'Wh2SmallSerialized.cpp'
        lexical = B.HERE/'../codec/WirehairSmallCore.h'
        raw = ('a.o: \\\n '+str(source)+' '+str(lexical)+' '+str(source)+'\n').encode()
        self.assertEqual(B.preprocessor_dependencies(raw, 'a.o'), {source, lexical.resolve()})
        for bad in [raw.replace(b'a.o:', b'b.o:', 1), b'a.o: relative.h\n', b'a.o: \n']:
            with self.assertRaises(ValueError): B.preprocessor_dependencies(bad, 'a.o')

    def test_context_symbol_is_unique_sized_and_known(self):
        for size in B.GF_BYTES.values():
            raw = ('00000000 %08x B GF256Ctx\n'%size).encode()
            self.assertEqual(B.context_size(raw), size)
            with self.assertRaises(ValueError): B.context_size(raw+raw)
        for bad in [b'', b'0 B GF256Ctx\n', b'0 00000000 B GF256Ctx\n', b'0 00022810 T GF256Ctx\n']:
            with self.assertRaises(ValueError): B.context_size(bad)

    def test_freeze_rejects_later_mutation(self):
        path = Path('/tmp/input'); first = dict(path=str(path), bytes=3, sha256='a'*64)
        with patch.object(B, 'pin', return_value=first):
            frozen = {}; B.freeze_inputs({path}, frozen); B.freeze_inputs({path}, frozen)
        with patch.object(B, 'pin', return_value=dict(first, sha256='b'*64)):
            with self.assertRaisesRegex(ValueError, 'input changed'): B.freeze_inputs({path}, frozen)

    def test_environment_excludes_implicit_compiler_and_loader_inputs(self):
        with patch.dict(os.environ, {'PATH':'/tmp/hostile', 'CPATH':'/tmp/includes', 'LIBRARY_PATH':'/tmp/libs',
                                     'LD_PRELOAD':'/tmp/preload', 'ASAN_OPTIONS':'checked', 'UBSAN_OPTIONS':'checked'}, clear=True):
            self.assertEqual(B.process_environment(), dict(PATH='/usr/bin:/bin', LANG='C', LC_ALL='C', TZ='UTC',
                                                           ASAN_OPTIONS='checked', UBSAN_OPTIONS='checked'))

    def test_qualified_inputs_rejects_scope_before_commands(self):
        with patch.object(B, 'command', side_effect=AssertionError('unexpected command')):
            for mode, path in [('wrong', Path('/tmp/native')), ('native', B.ROOT),
                               ('native', B.ROOT/'build'), ('native', Path('relative'))]:
                with self.assertRaises(ValueError): B.qualified_inputs(mode, path)
            for root in (B.PRODUCTION, B.PROOFS, B.PRIOR, B.SMALL, B.SEALED, B.BOUNDARY_AUDIT.parent,
                         Path('/var/tmp/other-science')):
                with self.assertRaises(ValueError): B.qualified_inputs('native',root/'new-native')

    def test_historical_header_not_rebound(self):
        self.assertNotEqual(B.CORE_OLD_SHA, B.CORE_SHA)
        self.assertEqual(B.CORE_OLD_SHA, '26167117230258275cc0d522c3039564f5ec5abe8009cbc48faae0e7b9f04d51')
        self.assertEqual(B.SOURCE_HEAD, '5a300505da3754dbdbb8cbfcba4b8a79bad1b0cb')
        self.assertEqual(B.CORE_SHA, '5b0acdd096d24b76351bacd1718c44ca5b37d4df587fe7334822a9f61f1e0b8c')


if __name__ == '__main__': unittest.main()
