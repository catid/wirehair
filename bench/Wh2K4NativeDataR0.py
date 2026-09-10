#!/usr/bin/env python3
"""Authenticated retained K4 fixtures; no selection, new traces or scoring."""
import argparse
import importlib.util
from pathlib import Path
import sys

SPEC=importlib.util.spec_from_file_location('_k4_native_reader',Path(__file__).with_name('Wh2K3NativeDataR0.py'))
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)
C,require,integer,ids=R.C,R.require,R.integer,R.ids
ROOT=Path('/var/tmp/wh2-k4-thue-morse-r0')
PROTOCOL='wirehair.wh2.k4-thue-morse-r0'
MANIFEST_SHA='8f1c9357b125edffc282df4cde29592130e27566f86946a7a85b0e92f175d281'
RAW_SHA='1e6422eeb1a62631fb9f0e8df86dac286b39ec7021cc5bbcb1d3f2335a351ce2'
LOOKUP_SHA='f4c3bc7d9fffac727acefc217d0481206b4c2d608af051993f34fb6a19f47c69'
FEEDBACK=(64,120,54,15)
PAIR=[[[int(r==c+1) if c<3 else FEEDBACK[r]^(phase if r==0 else 0)
        for c in range(4)] for r in range(4)] for phase in range(2)]
WIDTHS,CAP=R.WIDTHS,R.CAP


def load_report():
    raw=C.read_regular(ROOT/'COMPLETE.json',65536)
    require(C.sha(raw)==MANIFEST_SHA,'manifest identity')
    manifest=C.strict_json(raw)
    require(raw==C.canonical(manifest) and manifest['protocol']==PROTOCOL and
        manifest['outcome']=='PASS' and set(manifest['files'])==set(R.MEMBERS),'manifest schema')
    require({p.name for p in ROOT.iterdir()}==set(R.MEMBERS)|{'COMPLETE.json'},'bundle roster')
    content={};total=len(raw)
    for name in R.MEMBERS:
        path=ROOT/name;data=C.read_regular(path,CAP);total+=len(data)
        require(path.stat().st_mode&0o777==0o400 and
            manifest['files'][name]==dict(bytes=len(data),sha256=C.sha(data)),'member identity')
        content[name]=data
    require(total<=8*1024**2 and not content['stderr.txt'] and
        C.sha(content['raw.json'])==RAW_SHA,'raw/stderr identity')
    report=C.strict_json(content['raw.json'])
    require(content['raw.json']==C.canonical(report)+b'\n' and report['protocol']==PROTOCOL and
        report['outcome']=='PASS','report encoding/outcome')
    require(report['evidence']['lookup_bytes']==20480 and
        report['evidence']['lookup_sha256']==LOOKUP_SHA,'recorded lookup identity')
    return report


def extract(report):
    """Copy all retained ranks and original history widths/tails, without scoring."""
    require(report['pair']==PAIR,'sealed companion pair')
    require(len(report['fresh'])==6144 and len(report['hard'])==72,'trace roster')
    traces=[]
    for row in report['fresh']+report['hard']:
        require(type(row['B']) is int and row['B'] in WIDTHS and len(row['ranks'])==5,'trace shape')
        ranks=[integer(v,0,4) for v in row['ranks']]
        require(all(0<=b-a<=1 for a,b in zip(ranks,ranks[1:])),'nested prefix ranks')
        traces.append(dict(B=row['B'],ids=ids(row['ids'],8),ranks=ranks))
    require([sum(t['ranks'][oh]<4 for t in traces[:6144]) for oh in range(5)]==[9,0,0,0,0] and
        all(t['ranks']==[4]*5 for t in traces[6144:]),'retained fresh/hard rank counts')
    origins=report['inputs']['origins'];prefixes=report['inputs']['prefixes']
    require(len(origins)==38 and len(prefixes)==37 and len(report['history'])==37,'history roster')
    require(C.sha(C.canonical(origins))=='2efbce06e3db3928635764c51c74dd0f79966b7050b4751f592c60f14d75b578',
        'origin projection identity')
    require(C.sha(C.canonical(prefixes))=='96cb80f4458adca9f5cc6ed773f9cca27ba9e0b5a048dc0cd8c9f1ca11be42ff',
        'prefix projection identity')
    for recorded,prefix in zip(report['history'],prefixes):
        require(recorded==dict(prefix,rank=4) and type(recorded['rank']) is int,'history concordance')
    history=[]
    for origin in origins:
        require(type(origin['b']) is int and origin['b'] in WIDTHS,'original width')
        tail=integer(origin['tail'],1,origin['b'])
        packet_ids=ids(origin['ids'],integer(len(origin['ids']),4,8))
        require(any(p['ids']==packet_ids and origin['b'] in p['original_widths'] for p in prefixes),'origin coverage')
        history.append(dict(B=origin['b'],tail=tail,ids=packet_ids))
    require(sum(len(p['ids']) for p in history)==163,'history packet accounting')
    require(len(report['seams'])==31,'window roster')
    windows=[]
    for window,start in zip(report['seams'],[(1<<e)-4 for e in range(2,32)]+[2**32-8]):
        require(window['deficient']==[] and window['ids']==list(range(start,start+8)),'window certificate')
        windows.append(ids(window['ids'],8))
    rows=report['evidence']['unique_rows']
    require(len(rows)==2252 and [r['id'] for r in rows]==sorted({r['id'] for r in rows}),'row roster')
    for row in rows:
        integer(row['id'],0,2**32-1);require(len(row['row'])==4,'row shape')
        for value in row['row']:integer(value,0,255)
    require(len(traces)+len(history)+len(windows)*70==8424,'case accounting')
    require(len(traces)*8+sum(len(p['ids']) for p in history)+len(windows)*70*4==58571,'packet accounting')
    return dict(traces=traces,history=history,windows=windows,rows=rows)


def build_lookup():
    """Reconstruct only the audited fixed pair using polynomial multiplication."""
    table=[bytes(R.multiply(a,b) for b in range(256)) for a in range(256)]
    def product(a,b):
        output=[[0]*4 for _ in range(4)]
        for r in range(4):
            for c in range(4):
                for k in range(4):output[r][c]^=table[a[r][k]][b[k][c]]
        return output
    blocks=[PAIR]
    for level in range(31):
        a,b=blocks[level];blocks.append([product(a,b),product(b,a)])
    output=bytearray()
    for start,width,phase,vectors in ((0,10,0,True),(0,10,1,True),(10,7,0,False),
            (10,7,1,False),(17,7,0,False),(17,7,1,False),(24,8,0,False)):
        p=[[int(r==c) for c in range(4)] for r in range(4)]
        for index in range(1<<width):
            output.extend([r[0] for r in p] if vectors else [v for r in p for v in r])
            p=product(p,blocks[start][phase^(bin(index).count('1')&1)])
    require(len(output)==20480 and C.sha(output)==LOOKUP_SHA,'lookup bytes/hash')
    return bytes(output)


def render(data,lookup):
    require(len(lookup)==20480 and C.sha(lookup)==LOOKUP_SHA,'lookup bytes/hash')
    lines=['// Generated only from authenticated sealed K4 evidence.','#include <cstdint>',
        'namespace wh2_k4_data {','struct Trace { unsigned B; std::uint32_t ids[8]; unsigned ranks[5]; };',
        'struct Origin { unsigned B, tail, count; std::uint32_t ids[8]; };',
        'struct Row { std::uint32_t id; std::uint8_t values[4]; };',
        'static const char kRawSha[] = "'+RAW_SHA+'";','static const char kLookupSha[] = "'+LOOKUP_SHA+'";',
        'alignas(64) static const std::uint8_t kLookup[20480] = {']
    for offset in range(0,len(lookup),24):lines.append(','.join(str(v) for v in lookup[offset:offset+24])+',')
    lines.append('};\nstatic const Trace kTraces[] = {')
    for trace in data['traces']:
        lines.append('{%d,{%s},{%s}},'%(trace['B'],','.join(str(i)+'u' for i in trace['ids']),','.join(str(v) for v in trace['ranks'])))
    lines.append('};\nstatic const Origin kHistory[] = {')
    for prefix in data['history']:
        lines.append('{%d,%d,%d,{%s}},'%(prefix['B'],prefix['tail'],len(prefix['ids']),','.join(str(i)+'u' for i in prefix['ids'])))
    lines.append('};\nstatic const std::uint32_t kWindows[][8] = {')
    lines.extend('{'+','.join(str(i)+'u' for i in window)+'},' for window in data['windows'])
    lines.append('};\nstatic const Row kRows[] = {')
    lines.extend('{%du,{%s}},'%(row['id'],','.join(str(v) for v in row['row'])) for row in data['rows'])
    lines.append('};\n} // namespace wh2_k4_data\n')
    raw='\n'.join(lines).encode('ascii');require(len(raw)<=CAP,'generated header cap')
    return raw


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--header',required=True,type=Path)
    args=parser.parse_args();raw=render(extract(load_report()),build_lookup())
    R.write_header(args.header,raw)
    print('K4_DATA bytes=%d sha256=%s cases=8424 packets=58571'%(len(raw),C.sha(raw)))


if __name__=='__main__':
    try:main()
    except Exception as error:
        print(type(error).__name__+': '+str(error)[:1000],file=sys.stderr);sys.exit(1)
