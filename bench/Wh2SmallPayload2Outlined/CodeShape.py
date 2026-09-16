#!/usr/bin/env python3
"""Strict native x86-64 code-shape diagnostic, not a speed qualification."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import subprocess

HELPER = 'wh2_small_payload2_outlined::Apply(void*, void const* const*, unsigned char const*, int, int)'
ENCODER = '(anonymous namespace)::EncodeSmall((anonymous namespace)::PublicCodec*, unsigned int, void*, unsigned int)'
EXPECTED = {ENCODER: 3, 'wirehair_small_encode': 1}
PROTECTED = (ENCODER,
    '(anonymous namespace)::DecodeSmall((anonymous namespace)::PublicCodec*, unsigned int, void const*, unsigned int)',
    '(anonymous namespace)::RecoverSmall((anonymous namespace)::PublicCodec*, void*)',
    'wirehair_v2_encode', 'wirehair_v2_decode', 'wirehair_v2_recover')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def text_section(path):
    header = subprocess.check_output(['readelf', '-SW', str(path)], text=True)
    matches = re.findall(r'^\s*\[\s*\d+\]\s+\.text\s+PROGBITS\s+([0-9a-f]+)\s+([0-9a-f]+)\s+([0-9a-f]+)\s', header, re.M)
    require(len(matches) == 1, 'one .text section')
    address, offset, length = (int(x, 16) for x in matches[0])
    data = path.read_bytes()[offset:offset + length]
    require(len(data) == length, 'complete .text')
    return address, data


def symbols(path):
    listing = subprocess.check_output(['nm', '-S', '-C', '--defined-only', str(path)], text=True)
    result = {}
    for line in listing.splitlines():
        match = re.fullmatch(r'([0-9a-f]+) ([0-9a-f]+) [tT] (.*)', line)
        if match and match[3] in PROTECTED + (HELPER, 'gf256_mulset_multi_mem'):
            require(match[3] not in result, 'unique local function names')
            result[match[3]] = (int(match[1], 16), int(match[2], 16))
    return result


def redirects(path):
    listing = subprocess.check_output(['objdump', '-d', '-w', '-C', '-j', '.text', str(path)], text=True)
    result, owner = [], None
    for line in listing.splitlines():
        symbol = re.fullmatch(r'[0-9a-f]+ <(.*)>:', line)
        if symbol:
            owner = symbol[1]
        call = re.fullmatch(r'\s*([0-9a-f]+):\s+((?:[0-9a-f]{2}\s+)+)call\s+([0-9a-f]+) <(.*)>', line)
        if call and call[4] == HELPER:
            result.append((int(call[1], 16), int(call[3], 16), owner))
    return result


def compare(address, baseline, candidate, calls, original_target, protected):
    require(len(candidate) > len(baseline), 'new helper must append text')
    require(collections.Counter(owner for _, _, owner in calls) == EXPECTED, 'exact four redirected callers')
    require(len({pc for pc, _, _ in calls}) == 4, 'unique call sites')
    normalized = bytearray(candidate[:len(baseline)])
    for pc, target, _ in calls:
        offset = pc - address
        require(0 <= offset <= len(baseline) - 5, 'call site in baseline text')
        before, after = baseline[offset:offset + 5], candidate[offset:offset + 5]
        require(before[0] == after[0] == 0xe8, 'direct rel32 call opcode')
        require(pc + 5 + int.from_bytes(before[1:], 'little', signed=True) == original_target, 'original kernel target')
        require(pc + 5 + int.from_bytes(after[1:], 'little', signed=True) == target, 'new helper target')
        require(address + len(baseline) <= target < address + len(candidate), 'appended helper target')
        normalized[offset:offset + 5] = before
    for start, length in protected:
        offset = start - address
        require(length > 0 and 0 <= offset <= len(baseline) - length, 'protected function bounds')
        require(normalized[offset:offset + length] == baseline[offset:offset + length], 'protected function changed')
    # Do not silently normalize other relocations, function reordering, or
    # register allocation. Record residual differences even if the core matches.
    return [address + i for i, (a, b) in enumerate(zip(baseline, normalized)) if a != b]


def audit(baseline, candidate):
    address, before = text_section(baseline)
    candidate_address, after = text_section(candidate)
    require(candidate_address == address, 'same .text base')
    old_symbols, new_symbols = symbols(baseline), symbols(candidate)
    require('gf256_mulset_multi_mem' in old_symbols, 'original kernel symbol')
    require(HELPER in new_symbols, 'new helper symbol')
    require(all(name in old_symbols and name in new_symbols and
                old_symbols[name] == new_symbols[name] for name in PROTECTED), 'protected function address/size')
    calls = redirects(candidate)
    require(all(target == new_symbols[HELPER][0] for _, target, _ in calls), 'same named helper target')
    differences = compare(address, before, after, calls, old_symbols['gf256_mulset_multi_mem'][0],
                          [old_symbols[name] for name in PROTECTED])
    return dict(outcome='CORE_IDENTICAL_OTHER_TEXT_DIFFERS' if differences else 'SAME_TEXT_EXCEPT_FOUR_CALL_TARGETS',
                baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
                candidate_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest(),
                baseline_text_bytes=len(before), appended_text_bytes=len(after) - len(before),
                protected_functions={name:dict(address=hex(old_symbols[name][0]), bytes=old_symbols[name][1]) for name in PROTECTED},
                other_changed_bytes=[hex(pc) for pc in differences],
                redirects=[dict(pc=hex(pc), target=hex(target), caller=owner) for pc, target, owner in calls],
                speed_claimed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('candidate', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.baseline, args.candidate), indent=2, sort_keys=True))
