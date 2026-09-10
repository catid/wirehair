"""Explicit hash-checked derivation for the changed serializer binding graph.

No historical module or file is modified. The private derived module retains
the original receipt/current/run/replay globals and changes only metadata and
binding-header generation. Its exact source is an external preparation artifact.
"""
import hashlib
from pathlib import Path
import types

HERE = Path(__file__).resolve().parent
PREPARED = Path('/tmp/wh2-validated-serialization-inputs.zloiJHko')
BASELINE = Path('/tmp/wh2-v2-k8-admission.YNN2XmAx/native/libwirehair.so.2.0.0')
BASELINE_SHA = 'ef684abac606667e30e3b0de1204b0897aa0fa93a4f5ec2268d1ba951dce03b2'
ORIGINAL_SHA = '167e08af526f7527e6afa12276c465bbb9d844b126efd6bbdbfb78e923921dcb'
SLOTS = ('wirehair_recover_block_ex', 'wirehair_v2_encoder_create',
         'wirehair_v2_encoder_create_profile', 'wirehair_v2_encoder_create_profile_id',
         'wirehair_v2_profile_deserialize', 'wirehair_v2_profile_serialize')


def replace(text, old, new, count=1):
    if text.count(old) != count:
        raise ValueError('missing or ambiguous serialization runtime anchor: '+old)
    return text.replace(old, new)


def runtime_source():
    raw = (HERE/'Wh2AdmissionRegressionCostR0.py').read_bytes()
    if hashlib.sha256(raw).hexdigest() != ORIGINAL_SHA:
        raise ValueError('unchanged historical runtime source')
    text = raw.decode()
    text = replace(text, 'Path(__file__).with_name(filename)', 'Path('+repr(str(HERE))+') / filename')
    text = replace(text, 'def metadata(libraries=None):',
        'def validate_slot_roster(elf, index):\n'
        '    A.require(type(index) is int and index in (0,1), "explicit library slot index")\n'
        '    expected_slots = '+repr((SLOTS,SLOTS[:-1]))+'\n'
        '    A.exact(tuple(n for n,_ in sorted(elf.slots)), expected_slots[index], "exact per-library internal GOT names")\n'
        '    A.exact(len({p for _,p in elf.slots}),len(elf.slots), "unique internal GOT offsets")\n'
        '    A.require("wirehair_v2_profile_serialize" in elf.exports, "serializer remains publicly exported")\n'
        '\n\ndef metadata(libraries=None):')
    text = replace(text, 'def metadata(libraries=None):\n    result = []',
        'def metadata(libraries=None):\n'
        '    A.require(libraries is not None and len(libraries)==2, "explicit two-library serialization roster")\n'
        '    A.exact(tuple(str(p) for p,_ in libraries), '+
        repr((str(BASELINE),str(PREPARED/'libwirehair.so.2.0.0')))+', "exact serialization library paths and order")\n'
        '    A.exact(libraries[0][1], '+repr(BASELINE_SHA)+', "exact current library digest")\n'
        '    result = []')
    text = replace(text, "        A.exact(len(elf.slots),6,'exact internal public GOT roster')",
        '        validate_slot_roster(elf, len(result))')
    text = replace(text, "    A.exact([s['name'] for s in result[0]['runtime_slots']],",
        '    A.exact([s["name"] for s in result[0]["exports"]],\n'
        '            [s["name"] for s in result[1]["exports"]], "same exact public export names")\n'
        "    A.exact([s['name'] for s in result[0]['runtime_slots']],")
    text = replace(text, "entries(lib['exports']),entries(lib['slots']),entries(lib['runtime_slots'])",
        "entries(lib['exports']),entries(lib['slots']),str(len(lib['slots'])),entries(lib['runtime_slots'])")
    return text.encode()


def runtime(path):
    if path != PREPARED/'Wh2ValidatedSerializationDerivedR0.py':
        raise ValueError('exact derived runtime path')
    raw = runtime_source()
    if path.exists() and path.read_bytes() != raw:
        raise ValueError('unchanged published derived runtime')
    module = types.ModuleType('_validated_serialization_private_runtime')
    module.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), vars(module))
    if module.current.__globals__ is not vars(module):
        raise ValueError('receipt recomputation owns adapted runtime globals')
    return module


def cpp_source(text):
    text = replace(text, 'std::array<SymbolSpec,53> exports; std::array<SymbolSpec,6> slots;',
        'std::array<SymbolSpec,53> exports; std::array<SymbolSpec,6> slots; unsigned slot_count;')
    text = replace(text, '        for(const auto& s:spec.slots) {\n',
        '        Check(spec.slot_count==(index==0?6u:5u) && spec.slot_count<=spec.slots.size(), "exact internal GOT count");\n'
        '        for(unsigned i=0;i<spec.slot_count;++i) { const auto& s=spec.slots[i];\n')
    text = replace(text, "        for(const auto& s:spec.slots) { if(comma) putchar(','); comma=true; uintptr_t v=0;",
        '        Check(spec.slot_count==(a==0?6u:5u) && spec.slot_count<=spec.slots.size(), "exact published GOT count");\n'
        "        for(unsigned i=0;i<spec.slot_count;++i) { const auto& s=spec.slots[i]; if(comma) putchar(','); comma=true; uintptr_t v=0;")
    return text
