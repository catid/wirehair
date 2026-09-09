#!/usr/bin/env python3
"""Bounded, neutral-only pre/post DSO binding and payload audit; never a timer.

Run each load order in a fresh process. The unmodified libraries remain loaded
until process exit; no codec handle or function pointer outlives its provider.
This is not a static-library or shared-library performance qualification.
"""
import argparse
import ctypes as T
import importlib.util
import os
from pathlib import Path
import resource
import struct

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('_admission_neutral_io', HERE/'Wh2AlignedIntermediateCostR0.py')
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
LIBRARIES = (
    (Path('/tmp/wh2-small-production.BFlfvz/native/libwirehair.so.2.0.0'),
     '4eb4f53c6c7eba5edf56a377a80cb620bbb455971b0f2f50efa6ca2859f655d1'),
    (Path('/tmp/wh2-v2-k3-admission.TECbfU/native/libwirehair.so.2.0.0'),
     '0a45828773c67823221ea2bada230f3175cd9f7e3fc728d4d0d8b519564fd98e'))
ENV_KEYS = ('LD_PRELOAD', 'LD_LIBRARY_PATH', 'LD_AUDIT', 'LD_DEBUG',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES')
PROVIDERS = ('malloc', 'free', 'memcpy', '_Znwm', '_ZdlPv')
V, U32, U64, SIZE, INT = T.c_void_p, T.c_uint32, T.c_uint64, T.c_size_t, T.c_int


class Created(T.Structure):
    _fields_ = [('status', INT), ('codec', V)]


class Result(T.Structure):
    _fields_ = [('status', INT), ('required', U64), ('written', U64)]


class Options(T.Structure):
    _fields_ = [('bytes', U32), ('version', U32), ('policy', U32), ('reserved', U32)]


class DlInfo(T.Structure):
    _fields_ = [('path', T.c_char_p), ('base', V), ('name', T.c_char_p), ('symbol', V)]


class Elf:
    """Strict ELF64-LE x86-64 DSO metadata for authenticated local inputs."""
    def __init__(self, data):
        A.require(64 <= len(data) <= 4*1024**2, 'ELF size')
        h = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
        A.require(h[0][:7] == b'\x7fELF\x02\x01\x01' and h[1:4] == (3,62,1), 'ELF DSO identity')
        A.require(h[8] == 64 and h[11] == 64 and 0 < h[12] < 4096 and
                  h[6]+64*h[12] <= len(data), 'ELF section table')
        self.sections = [struct.unpack_from('<IIQQQQIIQQ', data, h[6]+64*i) for i in range(h[12])]
        self.data = data
        for s in self.sections:
            A.require(s[1] == 8 or s[4]+s[5] <= len(data), 'ELF section bounds')
        tables = {}
        for i, s in enumerate(self.sections):
            if s[1] not in (2,11):
                continue
            A.require(s[9] == 24 and s[5] % 24 == 0 and s[6] < len(self.sections), 'ELF symbol table')
            A.require(self.sections[s[6]][1] == 3, 'ELF symbol string table')
            strings = self.body(s[6]); symbols = []
            for name, info, other, section, value, size in struct.iter_unpack('<IBBHQQ', self.body(i)):
                end = strings.find(b'\0', name)
                A.require(name < len(strings) and end >= name, 'ELF symbol name')
                symbols.append((strings[name:end].decode('ascii'), info, section, value, size))
            tables[i] = symbols
        full = [tables[i] for i,s in enumerate(self.sections) if s[1] == 2]
        dynamic = [tables[i] for i,s in enumerate(self.sections) if s[1] == 11]
        A.require(len(full) == len(dynamic) == 1, 'one full/dynamic symbol table')
        self.defined = {}
        for symbol in full[0]:
            name, info, section, value, size = symbol
            if name and 0 < section < len(self.sections) and info & 15 in (1,2):
                self.defined.setdefault(name, []).append(symbol)
        self.exports = [s[0] for s in dynamic[0] if s[0].startswith('wirehair_') and s[2] != 0]
        A.require(len(self.exports) == len(set(self.exports)) and len(self.exports) == 53,
                  'exact public export roster')
        self.slots = []
        for i,s in enumerate(self.sections):
            if s[1] != 4:
                continue
            A.require(s[9] == 24 and s[5] % 24 == 0 and s[6] in tables, 'ELF RELA table')
            for offset, info, addend in struct.iter_unpack('<QQq', self.body(i)):
                symbols = tables[s[6]]; index = info >> 32
                A.require(index < len(symbols), 'ELF relocation symbol')
                name = symbols[index][0]
                if (info & 0xffffffff) == 7 and name.startswith('wirehair_'):
                    A.require(addend == 0 and offset % 8 == 0 and self.allocated(offset, 8, 3), 'wirehair GOT slot')
                    self.slots.append((name, offset))
        A.require(self.slots, 'nonempty internal public binding audit')

    def body(self, i):
        s = self.sections[i]
        return b'' if s[1] == 8 else self.data[s[4]:s[4]+s[5]]

    def allocated(self, offset, size, flags):
        return any(s[2] & flags == flags and s[3] <= offset and offset+size <= s[3]+s[5]
                   for s in self.sections)

    def symbol(self, name, kind):
        found = self.defined.get(name, [])
        A.require(len(found) == 1, 'unique private/public symbol: '+name)
        _, info, section, value, size = found[0]
        A.require(info & 15 == kind and size > 0 and self.allocated(value,size,6 if kind == 2 else 3),
                  'allocated symbol type/range')
        return value, size


def address(function):
    return T.cast(function, V).value


class Library:
    def __init__(self, index, libraries=None):
        self.path, digest = (LIBRARIES if libraries is None else libraries)[index]
        raw = A.read_regular(self.path, 4*1024**2)
        A.exact(A.sha(raw), digest, 'exact qualified DSO')
        self.elf = Elf(raw)
        self.global_scope = T.CDLL(None)
        A.require(not hasattr(self.global_scope, 'wirehair_init_'), 'no globally linked Wirehair')
        self.dso = T.CDLL(str(self.path), mode=os.RTLD_NOW | os.RTLD_LOCAL)
        dladdr = self.global_scope.dladdr
        dladdr.argtypes, dladdr.restype = [V,T.POINTER(DlInfo)], INT
        info = DlInfo()
        A.require(dladdr(address(self.dso.wirehair_init_),T.byref(info)) == 1 and
                  Path(os.fsdecode(info.path)).resolve() == self.path, 'owned DSO base')
        self.base = info.base
        self.report = dict(path=str(self.path), sha256=digest, base=self.base, exports={}, slots=[], providers={})
        for name in self.elf.exports:
            pointer = address(getattr(self.dso,name)); offset, size = self.elf.symbol(name,2)
            A.require(pointer == self.base+offset, 'public export provider: '+name)
            self.report['exports'][name] = dict(address=pointer, offset=offset, size=size)
        for name,offset in self.elf.slots:
            target = V.from_address(self.base+offset).value
            A.require(target == address(getattr(self.dso,name)), 'internal public symbol preemption: '+name)
            self.report['slots'].append(dict(name=name, offset=offset, target=target))
        for name in PROVIDERS:
            self.report['providers'][name] = address(getattr(self.dso,name))
        init = self.call('wirehair_init_', INT, INT)
        A.exact(init(2), 0, 'library initialization')
        ctx, ctx_size = self.elf.symbol('GF256Ctx',1)
        getter, _ = self.elf.symbol('gf256_get_active_x86_cpu_features',2)
        features = (INT*4)()
        T.CFUNCTYPE(None, T.POINTER(INT))(self.base+getter)(features)
        A.exact(list(features), [1]*4, 'active native arithmetic')
        self.report.update(context=self.base+ctx, context_bytes=ctx_size, features=list(features))
        self.check_bindings()

    def call(self, name, result, *arguments):
        A.require(name in self.report['exports'], 'declared public function')
        return T.CFUNCTYPE(result, *arguments)(self.report['exports'][name]['address'])

    def check_bindings(self):
        A.exact(A.sha(A.read_regular(self.path,4*1024**2)),self.report['sha256'],'unchanged DSO bytes')
        A.require(not hasattr(self.global_scope,'wirehair_init_'), 'Wirehair remains local')
        for name,r in self.report['exports'].items():
            A.exact(address(getattr(self.dso,name)),r['address'],'stable public export')
        for r in self.report['slots']:
            A.exact(V.from_address(self.base+r['offset']).value,r['target'],'stable internal GOT target')


def cases():
    for family in ('certified','wh1','small','k6'):
        ks = (3,) if family == 'small' else (6,) if family == 'k6' else (2,3,4,6,8,16,64,128)
        for k in ks:
            for b in (2,64,1280):
                for tail in (1,b):
                    for policy in (1,2):
                        yield family,k,b,tail,policy


def exercise(lib, case):
    family,k,b,tail,policy = case
    m = (k-1)*b+tail
    source_bytes = bytes((37*i+i//11)%256 for i in range(m))
    source = T.create_string_buffer(source_bytes,m)
    profile = (T.c_ubyte*32)(); encoder=V(); decoder=V()
    small = family in ('small','k6')
    prefix = 'wirehair_' + ('v2_' if family == 'certified' else family+'_' if small else '')
    free = lib.call(prefix+'free',None,V)
    packets=[]; statuses=[]
    high = [0xffffffff-2*j for j in range(6)]
    ids = list(range(k+8))+high
    try:
        if small:
            create = lib.call(prefix+'encoder_create',Created,V,U64,U32,U32,V,SIZE)
            r=create(source,m,b,policy,profile,32); encoder.value=r.codec
            A.exact(r.status,0,'small encoder create')
        elif family == 'certified':
            options=Options(16,1,policy,0); written=U32()
            create=lib.call(prefix+'encoder_create_profile_id_with_options',INT,U64,V,U64,U32,V,V,U32,V,V)
            A.exact(create(0x4b295bbb47f4f9c9,source,m,b,T.byref(options),profile,32,
                           T.byref(written),T.byref(encoder)),0,'certified encoder create')
            A.exact(written.value,32,'descriptor length')
        else:
            create=lib.call(prefix+('encoder_create_owned_ex' if policy == 1 else 'encoder_create_ex'),INT,V,V,U64,U32,V)
            A.exact(create(None,source,m,b,T.byref(encoder)),0,'WH1 encoder create')
        A.require(encoder.value is not None,'live encoder')
        if family == 'certified':
            A.exact(bytes(profile)[:16],struct.pack('<4sHHQ',b'WHV2',1,32,0x4b295bbb47f4f9c9),'explicit certified descriptor')
            A.exact(bytes(profile)[16:28],struct.pack('<QI',m,b),'certified dimensions')
            A.exact(bytes(profile)[29:],bytes(3),'certified reserved')
        elif small:
            magic,identity = (b'WHK3',0x5748324b33544d31) if family == 'small' else (b'WHK6',0x5748324b36544d31)
            A.exact(bytes(profile),struct.pack('<4sHHQQII',magic,1,32,identity,m,b,0),'opt-in descriptor')
        if policy == 1:
            T.memset(source,0xcc,m)
        for packet_id in ids:
            output=(T.c_ubyte*(b+2))(*([0xa5]*(b+2))); written=U32()
            if small:
                r=lib.call(prefix+'encode',Result,V,U32,V,SIZE)(encoder,packet_id,T.byref(output,1),b)
                A.exact(r.status,0,'small encode'); count=r.written
                A.exact(r.required,count,'small packet length')
            else:
                A.exact(lib.call(prefix+'encode',INT,V,U32,V,U32,V)(encoder,packet_id,T.byref(output,1),b,T.byref(written)),0,'encode')
                count=written.value
            A.exact(count,tail if packet_id == k-1 else b,'packet length')
            A.require(output[0] == 0xa5 and all(v == 0xa5 for v in output[1+count:]),'packet guards')
            packets.append(bytes(output[1:1+count]))
        free(encoder); encoder.value=None
        if small:
            r=lib.call(prefix+'decoder_create',Created,V,SIZE)(profile,32);decoder.value=r.codec
            A.exact(r.status,0,'small decoder create')
        elif family == 'certified':
            A.exact(lib.call(prefix+'decoder_create',INT,V,U32,V)(profile,32,T.byref(decoder)),0,'certified decoder create')
        else:
            A.exact(lib.call(prefix+'decoder_create_ex',INT,V,U64,U32,V)(None,m,b,T.byref(decoder)),0,'WH1 decoder create')
        A.require(decoder.value is not None,'live decoder')
        feed=lib.call(prefix+'decode',INT,V,U32,V,SIZE if small else U32)
        for slot in list(range(k,k+8))+list(range(k+8,k+14))+list(range(k)):
            packet=packets[slot]
            status=feed(decoder,ids[slot],packet,len(packet));A.require(status in (0,1),'feed status')
            statuses.append(status)
            if status == 0:
                break
        A.require(statuses[-1] == 0 and len(statuses) >= k,'successful endpoint')
        for _ in range(2):
            output=(T.c_ubyte*(m+2))(*([0xa5]*(m+2)));written=U64()
            if small:
                r=lib.call(prefix+'recover',Result,V,V,SIZE)(decoder,T.byref(output,1),m)
                A.exact((r.status,r.required,r.written),(0,m,m),'small recover')
            elif family == 'certified':
                A.exact(lib.call(prefix+'recover',INT,V,V,U64,V)(decoder,T.byref(output,1),m,T.byref(written)),0,'certified recover')
                A.exact(written.value,m,'recover length')
            else:
                A.exact(lib.call(prefix+'recover',INT,V,V,U64)(decoder,T.byref(output,1),m),0,'WH1 recover')
            A.require(output[0] == output[-1] == 0xa5 and bytes(output[1:-1]) == source_bytes,'recovered payload/guards')
        A.exact(source.raw,bytes([0xcc])*m if policy == 1 else source_bytes,'source unchanged')
        return dict(case=list(case),profile=bytes(profile).hex(),ids=ids,
                    packets=[A.sha(p) for p in packets],feed=statuses,first=len(statuses),recoveries=2)
    finally:
        if decoder.value:free(decoder)
        if encoder.value:free(encoder)


def probe(order):
    A.require(order in ('old-new','new-old'),'explicit load order')
    A.require(all(os.environ.get(k) is None for k in ENV_KEYS),'clean loader/allocator environment')
    A.exact((T.sizeof(V),T.sizeof(SIZE),T.sizeof(Created),T.sizeof(Result),T.sizeof(Options)),
            (8,8,16,24,16),'native C ABI sizes')
    libs={i:Library(i) for i in ((0,1) if order == 'old-new' else (1,0))}
    A.require(libs[0].report['context']+libs[0].report['context_bytes'] <= libs[1].report['context'] or
              libs[1].report['context']+libs[1].report['context_bytes'] <= libs[0].report['context'],
              'distinct private GF contexts')
    A.exact(libs[0].report['providers'],libs[1].report['providers'],'shared C/C++ allocator and memcpy providers')
    records=[]
    for case in cases():
        old,new=exercise(libs[0],case),exercise(libs[1],case)
        A.exact(old,new,'paired unchanged equations/payload/receive endpoint')
        records.append(old)
    for lib in libs.values():lib.check_bindings()
    A.exact(len(records),216,'whole neutral roster')
    return dict(protocol='wirehair.wh2.admission-regression-neutral.v1',order=order,
                libraries=[libs[i].report for i in (0,1)],records=records,
                measured=False,speed_claimed=False,recovery_rate_claimed=False,all_K_claimed=False)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order',required=True,choices=('old-new','new-old'))
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,30))
    resource.setrlimit(resource.RLIMIT_AS,(256*1024**2,256*1024**2))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    A.require(args.output.is_absolute(),'absolute output')
    output=args.output.parent.resolve(strict=True)/args.output.name
    A.require(output.parent.is_dir() and not output.exists() and not output.is_symlink(),'fresh external result')
    A.require(A.ROOT not in output.parents,'result outside repository')
    result=probe(args.order)
    A.publish(output,A.canonical(result))
    print('PASS neutral216 paired cases, both unmodified DSOs, public/GOT ownership; no timing')
