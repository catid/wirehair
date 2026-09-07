#!/usr/bin/env python3
"""Neutral-only isolation of the qualified .78 objects; never runs a timer."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]
QUALIFIED = Path('/tmp/wh2-ordered-release-build.VIQirG')
PREFIX = 'wh2_ordered_release_'
MEMBERS = tuple('WirehairV2' + n + '.cpp.o' for n in
                ('Codec', 'Peel', 'Plan', 'Policy', 'Precode', 'PrecodeDecode',
                 'PrecodeEncode', 'Profile', 'Seeds', 'Solve'))


def require(ok, why):
    if not ok:
        raise ValueError(why)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run(*args):
    return subprocess.check_output(list(map(str, args)), timeout=60)


class Elf:
    """Small, deliberately strict reader for these ELF64-LE x86-64 objects."""
    def __init__(self, data):
        self.data = data
        require(64 <= len(data) <= 64 * 1024**2, 'ELF size')
        h = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
        require(h[0][:7] == b'\x7fELF\x02\x01\x01' and h[1:4] == (1, 62, 1), 'ELF identity')
        require(h[8] == 64 and h[10] == 0 and h[11] == 64 and 0 < h[12] < 4096,
                'ordinary relocatable ELF headers')
        require(h[6] + 64*h[12] <= len(data) and 0 < h[13] < h[12], 'section table')
        self.sections = [struct.unpack_from('<IIQQQQIIQQ', data, h[6]+64*i) for i in range(h[12])]
        self.shstrings = h[13]
        for s in self.sections:
            require(s[1] == 8 or s[4]+s[5] <= len(data), 'section bounds')
        self.names = [self.string(self.body(h[13]), s[0]) for s in self.sections]
        tables = [i for i, s in enumerate(self.sections) if s[1] == 2]
        require(len(tables) == 1, 'single symtab')
        self.table = tables[0]
        symsection = self.sections[self.table]
        require(symsection[9] == 24 and symsection[5] % 24 == 0, 'symbol records')
        self.strings = symsection[6]
        require(0 < self.strings < len(self.sections) and self.sections[self.strings][1] == 3,
                'symbol strings')
        strings = self.body(self.strings)
        self.symbols = []
        for name, info, other, section, value, size in struct.iter_unpack('<IBBHQQ', self.body(self.table)):
            require(section < len(self.sections) or section in (0xfff1, 0xfff2), 'symbol section')
            self.symbols.append((self.string(strings, name), info, other, section, value, size))
        require(not any(s[1] >> 4 == 10 for s in self.symbols), 'no GNU unique symbols')

    @staticmethod
    def string(strings, start):
        require(start < len(strings), 'string start')
        end = strings.find(b'\0', start)
        require(end >= start, 'string terminator')
        return strings[start:end].decode('ascii')

    def body(self, index):
        s = self.sections[index]
        return b'' if s[1] == 8 else self.data[s[4]:s[4]+s[5]]

    def symbol(self, index, mapping):
        require(index < len(self.symbols), 'symbol index')
        s = self.symbols[index]
        return (mapping.get(s[0], s[0]),) + s[1:]

    def relocations(self, index, mapping):
        s = self.sections[index]
        require(s[6] == self.table and s[7] < len(self.sections), 'relocation linkage')
        width = 24 if s[1] == 4 else 16
        require(s[9] == width and s[5] % width == 0, 'relocation records')
        result = []
        for row in struct.iter_unpack('<QQq' if width == 24 else '<QQ', self.body(index)):
            offset, info = row[:2]
            result.append((offset, info & 0xffffffff, self.symbol(info >> 32, mapping)) + row[2:])
        return result

    def group(self, index, mapping):
        s = self.sections[index]
        require(s[6] == self.table and s[5] >= 8 and s[5] % 4 == 0, 'group linkage')
        words = struct.unpack('<' + 'I'*(s[5]//4), self.body(index))
        require(words[0] == 1 and all(0 < i < len(self.sections) for i in words[1:]), 'COMDAT members')
        return self.symbol(s[7], mapping), words


def mapping_for(objects, sanitizer=False):
    names = {s[0] for obj in objects for s in obj.symbols}
    require(not any(n.startswith(PREFIX) for n in names), 'already transformed input')
    selected = {n for n in names if 'wirehair_v2' in n}
    if sanitizer:
        # Instrumented std COMDAT code references private ASAN/UBSAN metadata.
        # Keep those groups separate in the neutral build, not in timing code.
        selected.update(s[0] for obj in objects for s in obj.symbols
                        if s[0] and s[1] >> 4 == 2 and s[3] != 0)
        selected.update(obj.symbols[s[7]][0] for obj in objects for s in obj.sections if s[1] == 17)
    return {n: PREFIX+n for n in sorted(selected)}


def audit_shared_groups(obj, mapping):
    """A deduplicable standard-library group must not borrow arm-private state."""
    shared = 0
    for i, section in enumerate(obj.sections):
        if section[1] != 17:
            continue
        signature, members = obj.group(i, {})
        if signature[0] in mapping:
            continue
        shared += 1
        for member in members[1:]:
            if obj.sections[member][1] not in (4, 9):
                continue
            for relocation in obj.relocations(member, {}):
                symbol = relocation[2]
                require(symbol[0] not in mapping, 'shared COMDAT references V2 symbol')
                if symbol[1] >> 4 == 0 and 0 < symbol[3] < len(obj.sections):
                    target = obj.sections[symbol[3]]
                    # GCC emits length_error message literals outside COMDAT.
                    # Mergeable, read-only strings with no relocations are not state.
                    literal = (target[1] == 1 and target[2] == 0x32 and target[9] == 1 and
                               not any(s[1] in (4, 9) and s[7] == symbol[3] for s in obj.sections))
                    require(not target[2] & 2 or symbol[3] in members[1:] or literal,
                            'shared COMDAT references private allocated section')
    return shared


def audit(before, after, mapping):
    require(before.names == after.names and before.table == after.table and
            before.strings == after.strings and before.shstrings == after.shstrings,
            'unchanged section identities')
    require(Counter(before.symbol(i, mapping) for i in range(len(before.symbols))) ==
            Counter(after.symbols), 'exact mapped symbol metadata')
    for i, (a, b) in enumerate(zip(before.sections, after.sections)):
        # File offsets and string-table sizes may change; allocated layout may not.
        require(a[1:4] == b[1:4] and a[6] == b[6] and a[8:] == b[8:], 'section metadata')
        if i not in (before.strings, before.shstrings):
            require(a[5] == b[5], 'section size')
        if a[1] == 17:
            require(before.group(i, mapping) == after.group(i, {}), 'isolated COMDAT identity')
        elif a[1] in (4, 9):
            require(a[7] == b[7] and before.relocations(i, mapping) == after.relocations(i, {}),
                    'exact relocation semantics')
        else:
            require(a[7] == b[7], 'section info')
            if i not in (before.table, before.strings, before.shstrings):
                require(before.body(i) == after.body(i), 'unchanged code/data/debug bytes: '+before.names[i])
    return dict(symbols=len(before.symbols), sections=len(before.sections),
                mapped=sum(s[0] in mapping for s in before.symbols),
                groups=sum(s[1] == 17 for s in before.sections),
                allocated_bytes=sum(s[5] for s in before.sections if s[2] & 2))


def pin(path):
    require(path.is_file() and not path.is_symlink(), 'regular artifact: '+str(path))
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path.read_bytes()))


def write(path, data):
    with path.open('xb') as stream:
        stream.write(data)


def build(mode, output):
    require(output.is_absolute(), 'absolute output')
    output = output.parent.resolve(strict=True)/output.name
    require(ROOT not in output.parents and output != ROOT, 'external output')
    require(not output.exists() and not output.is_symlink(), 'fresh output')
    require(mode in ('native', 'asan', 'scalar'), 'qualified mode')
    qualified = QUALIFIED / mode
    # Authenticate the retained scientific input closure, deliberately not its old HEAD.
    receipt_bytes = (QUALIFIED / 'receipt.json').read_bytes()
    require(sha(receipt_bytes) == '3057037621a58d4a00708bb9426119865edffaf53a8e4b0d67f75fc1c3cf2a89',
            'qualified receipt')
    spec = importlib.util.spec_from_file_location('release_prior', ROOT/'bench/Wh2OrderedRhsReleaseR0.py')
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    for p in json.loads(receipt_bytes)['pins']:
        require(prior.pin_input(Path(p['path'])) == p, 'qualified input changed: '+p['path'])
    require(all(prior.os.environ.get(k) is None for k in prior.ENVIRONMENT_KEYS), 'allocator environment')
    require(sha((qualified/'WirehairV2SolveOrdered.cpp').read_bytes()) ==
            '0c0b447c6273edd861176d10be378312d9fcb525d0f0bc3a4fe76a3df06fd4ac', 'five-swap source')
    archive = qualified/'production/libwirehair.a'
    paths = [qualified/'production/CMakeFiles/wirehair.dir/codec'/n for n in MEMBERS]
    # The nonsolver member instructions are exactly those of the baseline archive.
    for path in paths:
        require(run('/usr/bin/ar', 'p', archive, path.name) == path.read_bytes(), 'archive member parity')
    paths[-1] = qualified/'CMakeFiles/release_candidate.dir/WirehairV2SolveOrdered.cpp.o'
    original = [Elf(p.read_bytes()) for p in paths]
    mapping = mapping_for(original, sanitizer=mode == 'asan')
    require(mapping, 'nonempty symbol map')
    shared_groups = [audit_shared_groups(obj, mapping) for obj in original]
    definitions = {s[0] for obj in original for s in obj.symbols if s[3] != 0 and s[1] >> 4 in (1, 2)}
    references = {s[0] for obj in original for s in obj.symbols if s[3] == 0 and s[0] in mapping}
    require(references <= definitions, 'closed candidate symbol references')
    public = sorted({s[0] for obj in original for s in obj.symbols
                     if s[0].startswith('wirehair_v2_') and s[1] >> 4 == 1 and s[3] != 0})
    require(len(public) == 16 and all(re.fullmatch(r'wirehair_v2_[a-z_]+', n) for n in public),
            'exact public exports')
    # No non-V2 strong definitions may collide with the unchanged archive.
    require(not any(s[0] not in mapping and s[1] >> 4 == 1 and s[3] != 0
                    for obj in original for s in obj.symbols), 'no shared strong definitions')
    output.mkdir(mode=0o700)
    mapfile = output/'symbols.map'
    write(mapfile, ''.join(a+' '+b+'\n' for a, b in mapping.items()).encode())
    transformed = []
    audits = []
    for name, path, obj in zip(MEMBERS, paths, original):
        destination = output/name
        run('/usr/bin/objcopy', '--redefine-syms='+str(mapfile), path, destination)
        audits.append(audit(obj, Elf(destination.read_bytes()), mapping))
        transformed.append(destination)
    header = '#pragma once\n#include "wirehair/wirehair.h"\nextern "C" {\n'
    header += ''.join('extern decltype('+n+') '+mapping[n]+';\n' for n in public) + '}\n'
    write(output/'OrderedReleaseApi.h', header.encode())
    source = ROOT/'bench/Wh2OrderedRhsReleaseIsolationTest.cpp'
    executable = output/'isolation_test'
    command = ['/usr/bin/c++', '-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror',
               '-fno-lto', '-fPIC', '-no-pie', '-DWIREHAIR_STATIC=1',
               '-I'+str(ROOT/'include'), '-I'+str(output)]
    if mode == 'asan':
        command += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer']
    else:
        command += ['-O3', '-g1']
    if mode == 'scalar':
        command += ['-DANDROID=1']
    command += [str(source)] + list(map(str, transformed)) + [str(archive), '-pthread',
               '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    run(*command)
    symbols = run('/usr/bin/nm', '-g', '--defined-only', executable).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    require(names.count('GF256Ctx') == 1 and names.count('gf256_init_') == 1, 'one GF context/init')
    for n in public:
        require(names.count(n) == names.count(mapping[n]) == 1, 'distinct complete C APIs')
    require('libwirehair.a(WirehairV2Solve.cpp.o)' in (output/'link.map').read_text(), 'real baseline solver')
    for p in json.loads(receipt_bytes)['pins']:
        require(prior.pin_input(Path(p['path'])) == p, 'qualified input changed during build: '+p['path'])
    manifest = dict(schema='wirehair.wh2.ordered-release-isolation.v1', mode=mode,
                    scientific_worker=False, command=command, audits=audits,
                    shared_groups=shared_groups, closed_references=len(references),
                    inputs=[pin(p) for p in paths+[archive, source, Path(__file__).resolve()]],
                    artifacts=[pin(p) for p in sorted(output.iterdir())])
    write(output/'manifest.json', (json.dumps(manifest, sort_keys=True)+'\n').encode())
    print(json.dumps(dict(mode=mode, mapped_names=len(mapping), objects=len(audits),
                         executable=str(executable), scientific_worker=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('native', 'asan', 'scalar'))
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    build(args.mode, args.output)
