#!/usr/bin/env python3
"""Describe sealed gate-A records. Never executes a codec or changes its verdict."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import statistics
import struct

ROOT = Path('/home/catid/wirehair')
BUNDLE = Path('/var/tmp/wh2-k8-ordinary-shared-retention-r0')
AUDIT = Path('/tmp/wh2-k8-retention-auditor.jY1sGYtx/POSTRUN.json')
AUDIT_SHA = 'f030a6f8e27005dacdf818fd0d9e402f03f872c298263a6ae94ea5185ff8e4ec'
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
COUNTERS = ('minor_faults','major_faults','voluntary_switches','involuntary_switches')
METRICS = ('work_ns','thread_bracket_ns','capture_wall_ns','wait_wall_ns',
           'wait_thread_ns','wait_overshoot_ns','gap_wall_ns','gap_thread_ns')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def decode(raw):
    def pairs(items):
        result = {}
        for key,value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def constant(value):
        raise ValueError('nonfinite JSON constant: '+value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def canonical(value):
    return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


def identity(info):
    return (info.st_dev,info.st_ino,info.st_mode,info.st_nlink,info.st_size,
            info.st_mtime_ns,info.st_ctime_ns)


def sealed_lines(pin, cap=192*1024**2):
    """Check the full byte stream and stable regular-file identity at EOF."""
    path = Path(pin['path'])
    fd = os.open(path, os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and stat.S_IMODE(before.st_mode)==0o400
                and before.st_nlink==1 and before.st_size==pin['bytes']
                and before.st_size<=cap, 'sealed bounded regular input')
        digest = hashlib.sha256(); size = 0
        with os.fdopen(fd,'rb',closefd=False) as stream:
            for line in iter(lambda:stream.readline(4*1024**2+1),b''):
                require(len(line)<=4*1024**2 and line.endswith(b'\n'), 'bounded complete JSON line')
                size += len(line); digest.update(line)
                require(size<=cap, 'bounded stream')
                yield line
        require(size==pin['bytes'] and digest.hexdigest()==pin['sha256'], 'exact retained file hash')
        require(identity(before)==identity(os.fstat(fd)), 'stable retained file')
    finally:
        os.close(fd)


def sealed_json(pin):
    return decode(b''.join(sealed_lines(pin,2*1024**2)))


def annotation(row, previous):
    c = row['coordinate']; o = row['observation']; clocks = o['clocks']; w = row['wait']
    require(len(c)==9 and len(clocks)==6 and len(w)==4, 'named record schema')
    require(row['complete'] is True and row['checked'] is True, 'checked complete work')
    require(previous['clocks'][5]<=clocks[0]<=clocks[2]<clocks[3]<=clocks[5]
            and previous['clocks'][4]<=clocks[1]<=clocks[4], 'clock chronology')
    require(row['ready']<=w[0]<=w[2]<=clocks[0] and row['target']<=w[2]
            and w[1]<=w[3]<=clocks[1], 'wait chronology')
    delta = [a-b for a,b in zip(o['after'],o['before'])]
    gap = [a-b for a,b in zip(o['before'],previous['after'])]
    require(len(delta)==4 and len(gap)==4 and min(delta+gap)>=0, 'counter chronology')
    addresses = row['addresses'][:row['address_count']]
    require(addresses and all(type(a) is int and 0<a<2**64 for a in addresses), 'retained addresses')
    require(all(a==0 for a in row['addresses'][row['address_count']:]), 'unused addresses')
    packed = struct.pack('<'+'Q'*len(addresses),*addresses)
    return dict(coordinate=c, logical_side=SIDES[c[6]]^c[2], prelude=c[6]<2,
        work_start_ns=clocks[2], work_ns=clocks[3]-clocks[2],
        thread_bracket_ns=clocks[4]-clocks[1],capture_wall_ns=clocks[5]-clocks[0],
        wait_wall_ns=w[2]-w[0],wait_thread_ns=w[3]-w[1],wait_overshoot_ns=w[2]-row['target'],
        gap_wall_ns=clocks[0]-previous['clocks'][5],gap_thread_ns=clocks[1]-previous['clocks'][4],
        capture_counters=delta,gap_counters=gap,
        address_count=len(addresses),first_address=addresses[0],distinct_batch_addresses=len(set(addresses)),
        address_vector_sha256=hashlib.sha256(packed).hexdigest(),
        all_handle_mod64=dict(Counter(a%64 for a in addresses)))


def distribution(values):
    require(bool(values), 'nonempty descriptive series')
    return dict(count=len(values),minimum=min(values),median=statistics.median(values),
                maximum=max(values),total=sum(values))


def side_summary(rows):
    require(bool(rows), 'nonempty side')
    offsets = Counter()
    for r in rows: offsets.update(r['all_handle_mod64'])
    return dict(records=len(rows),preludes=sum(r['prelude'] for r in rows),
        metrics={name:distribution([r[name] for r in rows]) for name in METRICS},
        capture_counters=dict(zip(COUNTERS,[sum(r['capture_counters'][i] for r in rows) for i in range(4)])),
        gap_counters=dict(zip(COUNTERS,[sum(r['gap_counters'][i] for r in rows) for i in range(4)])),
        captures_with_counter_change=sum(any(r['capture_counters']) for r in rows),
        unique_address_vectors=len({r['address_vector_sha256'] for r in rows}),
        first_addresses=sorted({r['first_address'] for r in rows}),
        first_address_mod4096=dict(Counter(r['first_address']%4096 for r in rows)),
        all_handle_mod64=dict(offsets),
        distinct_batch_addresses=distribution([r['distinct_batch_addresses'] for r in rows]))


def describe(rows, stored, expected_count=98496):
    """Visit every record. Extra focus records never replace the all-cell report."""
    rows = iter(rows); header = next(rows)
    require(header['type']=='header', 'header')
    cases = [f['case'] for f in header['fixtures']]
    expected = {}
    for s in stored['statistics']:
        key = (cases.index(s['case']),s['metric'],s['comparison'],s['order'])
        require(key not in expected, 'unique stored cell')
        expected[key] = s
    failed = {tuple(k) for k in stored['failed_controls']}
    require(failed<=set(expected), 'every failed control is a declared cell')
    groups = defaultdict(list); previous = header['prelude']; total = 0
    for index in range(expected_count):
        row = next(rows)
        require(row['type']=='record' and row['coordinate'][0]==index, 'complete sequential raw chronology')
        a = annotation(row,previous); previous = row['observation']; total += a['work_ns']
        c = a['coordinate']; key = (c[3],c[4],c[5],c[2])
        require(key in expected, 'declared cell')
        groups[key].append(a)
    footer = next(rows)
    require(footer==dict(type='footer',complete=True,records=expected_count,work_ns=total), 'exact footer')
    sentinel = object(); require(next(rows,sentinel) is sentinel, 'no extra raw data')
    require(set(groups)==set(expected), 'all stored cells covered')
    cells = []; focused = []
    for key,group in sorted(groups.items()):
        require(len(group)==216, 'all twelve eighteen-position panels')
        by_rep = defaultdict(list)
        for r in group: by_rep[r['coordinate'][1]].append(r)
        require(set(by_rep)==set(range(12)), 'all twelve replicates')
        logs = []; panels = []
        for rep,panel in sorted(by_rep.items()):
            require([r['coordinate'][6] for r in panel]==list(range(18)), 'full panel positions')
            contrasts = []
            for p in range(2,18,2):
                pair = panel[p:p+2]
                require(pair[0]['coordinate'][8]==pair[1]['coordinate'][8], 'same delay in pair')
                values = {r['logical_side']:r['work_ns'] for r in pair}
                require(set(values)=={0,1}, 'opposite logical sides')
                contrasts.append(math.log(values[1])-math.log(values[0]))
            logs.append(math.fsum(contrasts)/8)
            panels.append(dict(replicate=rep,paired_log=logs[-1],records=panel))
        require(logs==expected[key]['replicate_logs'], 'unchanged original paired statistics')
        sides = [[r for r in group if r['logical_side']==side] for side in range(2)]
        require(all(len(side)==108 for side in sides), 'full logical sides')
        cells.append(dict(key=list(key),case=cases[key[0]],
            original_statistics=expected[key],sides=[side_summary(side) for side in sides],
            same_address_vector_pairs=sum(group[p]['address_vector_sha256']==group[p+1]['address_vector_sha256']
                                         for p in range(0,len(group),2)),
            measured_pairs=96,all_pairs_including_preludes=108))
        if key in failed:
            focused.append(dict(key=list(key),case=cases[key[0]],panels=panels))
    return dict(original_outcome=stored['outcome'],records=expected_count,work_ns=total,
        binding_bases=[b['base'] for b in header['bindings']],
        gf_context_addresses=[b['context'] for b in header['bindings']],
        cells=cells,failed_control_panels=focused)


def analyze():
    audit = sealed_json(dict(path=str(AUDIT),bytes=4792,sha256=AUDIT_SHA))
    require(audit['audit_valid'] is True and audit['scientific_outcome']=='CONTROL_FAIL', 'accepted terminal audit')
    pins = {Path(p['path']).name:p for p in audit['scientific_files']}
    require(len(pins)==8 and all(Path(p['path']).parent==BUNDLE for p in pins.values()), 'fixed bundle members')
    complete = sealed_json(pins['COMPLETE.json'])
    require(complete['outcome']=='CONTROL_FAIL' and len(complete['files'])==7, 'unchanged terminal decision')
    require({p['path']:p for p in complete['files']}==
            {p['path']:p for name,p in pins.items() if name!='COMPLETE.json'}, 'complete audit member agreement')
    stored = sealed_json(pins['analysis.json'])
    for name in ('CLAIM.json','processes.json'):
        sealed_json(pins[name])
    result = []
    for index,name in enumerate(('old-new','new-old')):
        require(b''.join(sealed_lines(pins[name+'.stderr.txt']))==b'', 'empty successful stderr')
        result.append(describe((decode(line) for line in sealed_lines(pins[name+'.raw.jsonl'])),
                               stored['load_orders'][index]))
    return dict(schema=1,diagnostic='retained gate-A clocks, addresses and controls',
        scientific_outcome_unchanged='CONTROL_FAIL',codec_calls=0,new_timing_observations=0,
        rescore_performed=False,retention_qualified=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        accepted_audit_sha256=AUDIT_SHA,scientific_files=list(pins.values()),load_orders=result,
        caveats=['Thread brackets include clock/capture work and are not exclusive codec CPU cost.',
                 'Counter deltas cover capture brackets, not just the inner WORK timestamps.',
                 'Inter-capture gaps include output checking, preparation and waiting; gap and wait totals are nested, not additive.',
                 'Gap counter deltas cannot attribute a fault or switch specifically to the wait.',
                 'Virtual handle addresses omit internal allocations and do not identify physical pages.',
                 'All records, preludes and cells are included; focused panels only expand failed controls.',
                 'Descriptive distributions and associations do not change the frozen verdict or prove causation.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    output = args.output.parent.resolve(strict=True)/args.output.name
    require(output.is_absolute() and ROOT not in output.parents and BUNDLE not in output.parents
            and output!=ROOT and output!=BUNDLE, 'external diagnostic output only')
    require(not output.exists() and not output.is_symlink(), 'fresh diagnostic output')
    result = analyze(); raw = canonical(result)
    fd = os.open(output,os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600)
    try:
        with os.fdopen(fd,'wb',closefd=False) as stream:
            require(stream.write(raw)==len(raw), 'complete diagnostic publication')
            stream.flush(); os.fsync(fd)
        os.fchmod(fd,0o400)
    finally:
        os.close(fd)
    print(json.dumps(dict(output=str(output),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
                          records=sum(r['records'] for r in result['load_orders']),outcome_unchanged='CONTROL_FAIL')))


if __name__=='__main__': main()
