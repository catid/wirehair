#!/usr/bin/env python3
"""Authenticated retained K2 fixtures; no selection, new traces or scoring."""
import argparse
import importlib.util
from pathlib import Path
import sys

SPEC=importlib.util.spec_from_file_location('_k2_native_reader',Path(__file__).with_name('Wh2K3NativeDataR0.py'))
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)
C,require,integer,ids=R.C,R.require,R.integer,R.ids
ROOT=Path('/var/tmp/wh2-k2-thue-morse-r0')
PROTOCOL='wirehair.wh2.k2-thue-morse-r0'
MANIFEST_SHA='7b76bc8584ccc175d45bce9e1b1dbc36ec6c849f9ec8356f0a137b2854766a75'
RAW_SHA='984d0836e5d297c72ffa568527d0d23cea3da5588d554ec715fc69f331aa5e73'
LOOKUP_SHA='c0529964214c5032bfc3d5b7d407ebc318bf23232ad60f694d9af1449cbdfd7a'
PAIR=(((0,2),(1,3)),((0,3),(1,3)))
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
        p=ROOT/name;data=C.read_regular(p,CAP);total+=len(data)
        require(p.stat().st_mode&0o777==0o400 and manifest['files'][name]==dict(bytes=len(data),sha256=C.sha(data)),'member identity')
        content[name]=data
    require(total<=8*1024**2 and not content['stderr.txt'] and C.sha(content['raw.json'])==RAW_SHA,'raw/stderr identity')
    report=C.strict_json(content['raw.json'])
    require(content['raw.json']==C.canonical(report)+b'\n' and report['protocol']==PROTOCOL and report['outcome']=='PASS','report encoding/outcome')
    require(report['evidence']['lookup_bytes']==7168 and report['evidence']['lookup_sha256']==LOOKUP_SHA,'recorded lookup identity')
    return report


def extract(report):
    """Copy original widths/tails and recorded ranks, including deficient pairs."""
    require(report['pair']==[[list(row) for row in m] for m in PAIR],'sealed companion pair')
    require(len(report['fresh'])==6144 and len(report['hard'])==72,'trace roster')
    traces=[]
    for row in report['fresh']+report['hard']:
        require(type(row['B']) is int and row['B'] in WIDTHS and row['ranks']==[2]*5 and
            all(type(v) is int for v in row['ranks']),'retained full-rank trace')
        traces.append(dict(B=row['B'],ids=ids(row['ids'],6),ranks=list(row['ranks'])))
    origins=report['inputs']['origins'];prefixes=report['inputs']['prefixes']
    require(len(origins)==56 and len(prefixes)==13 and len(report['history'])==13,'history roster')
    require(C.sha(C.canonical(origins))=='74361b5ac2588afd50a58e9a0ccf4c6b10f5f464e577f7c3d3c9daaeaa354516','origin projection identity')
    require(C.sha(C.canonical(prefixes))=='2e0968c8ffb53546547cb4034c6a9c4c0d247981f6a18825f44b258a2d926a8b','prefix projection identity')
    for recorded,p in zip(report['history'],prefixes):
        require(recorded==dict(p,rank=2) and type(recorded['rank']) is int,'history concordance')
    history=[]
    for o in origins:
        require(type(o['b']) is int and o['b'] in WIDTHS,'original width')
        tail=integer(o['tail'],1,o['b']);packet_ids=ids(o['ids'],2)
        require(any(p['ids']==packet_ids and o['b'] in p['original_widths'] for p in prefixes),'origin coverage')
        history.append(dict(B=o['b'],tail=tail,ids=packet_ids))
    require(len(report['seams'])==30,'window roster')
    windows=[]
    for w in report['seams']:
        require(w['deficient']==[],'window certificate');windows.append(ids(w['ids'],6))
    require(len(report['legacy'])==3 and len(report['strides'])==3,'pair roster')
    pairs=[]
    for p,expected in zip(report['legacy'],([265,270],[2,966],[1056,1313])):
        require(p['ids']==expected and integer(p['rank'],2,2)==2,'legacy witness')
        pairs.append(dict(ids=ids(p['ids'],2),rank=2))
    for entry,stride,failures in zip(report['strides'],(255,257,65537),(1,3,1)):
        require(entry['stride']==stride and entry['failures']==failures and entry['passed'] is True and
            len(entry['pairs'])==512,'stride cell')
        cell=[]
        for p in entry['pairs']:
            packet_ids=ids(p['ids'],2);rank=integer(p['rank'],1,2);determinant=integer(p['determinant'],0,255)
            require(packet_ids[1]-packet_ids[0]==stride and (rank==2)==(determinant!=0),'stride rank/determinant')
            cell.append(dict(ids=packet_ids,rank=rank))
        require(len({tuple(p['ids']) for p in cell})==512 and sum(p['rank']==1 for p in cell)==failures,'stride coverage')
        pairs.extend(cell)
    rows=report['evidence']['unique_rows']
    require(len(rows)==5274 and [r['id'] for r in rows]==sorted({r['id'] for r in rows}),'row roster')
    for row in rows:
        integer(row['id'],0,2**32-1);require(len(row['row'])==2,'row shape')
        for value in row['row']:integer(value,0,255)
    require(len(traces)+len(history)+len(windows)*15+len(pairs)*6==15956,'case accounting')
    require(len(traces)*6+len(history)*2+len(windows)*30+len(pairs)*12==56776,'packet accounting')
    return dict(traces=traces,history=history,windows=windows,pairs=pairs,rows=rows)


def build_lookup():
    table=[bytes(R.multiply(a,b) for b in range(256)) for a in range(256)]
    def product(a,b):
        return tuple(tuple(table[a[r][0]][b[0][c]]^table[a[r][1]][b[1][c]] for c in range(2)) for r in range(2))
    blocks=[PAIR]
    for level in range(31):
        a,b=blocks[level];blocks.append((product(a,b),product(b,a)))
    output=bytearray()
    for start,width,phase,vectors in ((0,10,0,True),(0,10,1,True),(10,7,0,False),
            (10,7,1,False),(17,7,0,False),(17,7,1,False),(24,8,0,False)):
        p=((1,0),(0,1))
        for index in range(1<<width):
            output.extend([r[0] for r in p] if vectors else [v for r in p for v in r])
            p=product(p,blocks[start][phase^(bin(index).count('1')&1)])
    require(len(output)==7168 and C.sha(output)==LOOKUP_SHA,'lookup bytes/hash')
    return bytes(output)


def render(data,lookup):
    require(len(lookup)==7168 and C.sha(lookup)==LOOKUP_SHA,'lookup bytes/hash')
    lines=['// Generated only from authenticated sealed K2 evidence.','#include <cstdint>',
        'namespace wh2_k2_data {','struct Trace { unsigned B; std::uint32_t ids[6]; unsigned ranks[5]; };',
        'struct Origin { unsigned B, tail, count; std::uint32_t ids[2]; };',
        'struct Pair { std::uint32_t ids[2]; unsigned rank; };',
        'struct Row { std::uint32_t id; std::uint8_t values[2]; };',
        'static const char kRawSha[] = "'+RAW_SHA+'";','static const char kLookupSha[] = "'+LOOKUP_SHA+'";',
        'alignas(64) static const std::uint8_t kLookup[7168] = {']
    for offset in range(0,len(lookup),24):lines.append(','.join(str(v) for v in lookup[offset:offset+24])+',')
    lines.append('};\nstatic const Trace kTraces[] = {')
    for t in data['traces']:
        lines.append('{%d,{%s},{%s}},'%(t['B'],','.join(str(i)+'u' for i in t['ids']),','.join(str(v) for v in t['ranks'])))
    lines.append('};\nstatic const Origin kHistory[] = {')
    for p in data['history']:
        lines.append('{%d,%d,%d,{%s}},'%(p['B'],p['tail'],len(p['ids']),','.join(str(i)+'u' for i in p['ids'])))
    lines.append('};\nstatic const std::uint32_t kWindows[][6] = {')
    lines.extend('{'+','.join(str(i)+'u' for i in w)+'},' for w in data['windows'])
    lines.append('};\nstatic const Pair kPairs[] = {')
    lines.extend('{{%s},%d},'%(','.join(str(i)+'u' for i in p['ids']),p['rank']) for p in data['pairs'])
    lines.append('};\nstatic const Row kRows[] = {')
    lines.extend('{%du,{%s}},'%(r['id'],','.join(str(v) for v in r['row'])) for r in data['rows'])
    lines.append('};\n} // namespace wh2_k2_data\n')
    raw='\n'.join(lines).encode('ascii');require(len(raw)<=CAP,'generated header cap')
    return raw


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--header',required=True,type=Path)
    args=parser.parse_args();raw=render(extract(load_report()),build_lookup())
    R.write_header(args.header,raw)
    print('K2_DATA bytes=%d sha256=%s cases=15956 packets=56776'%(len(raw),C.sha(raw)))


if __name__=='__main__':
    try:main()
    except Exception as error:
        print(type(error).__name__+': '+str(error)[:1000],file=sys.stderr);sys.exit(1)
