#!/usr/bin/env python3
"""Read-only K4 history reconstruction, not producer or speed qualification.

No project imports, codec execution, historical receipt rebinding, resampling,
or selection of a favorable run. Every fixed retained namespace is inspected.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = 'wirehair.wh2.k4-serialized-cost-r0'
BUNDLES = (
    ('R0', Path('/var/tmp/wh2-k4-serialized-cost-r0'), 'f1504ab7422cb9a64164f45b0d46ef74684d34c0c76f185a1558107646a8d113'),
    ('R17', Path('/var/tmp/wh2-k4-serialized-cost-r1.R17/science'), '570a347bbaf26e7bd05e291ad894691557c8fd89498096cd6c722b1ee9adbad8'),
    ('R18', Path('/var/tmp/wh2-k4-serialized-cost-r1.R18/science'), 'c3834e908db596709f98b8374a414ebeff47fe055609716f30044cf5502be9b8'),
    ('R22', Path('/var/tmp/wh2-k4-serialized-cost-r1.R22/science'), '9f8a524662cf99e4485ad1ed2a5cd05ab0d538ec4a883666f68bd87ac2a2bdf7'),
    ('R23', Path('/var/tmp/wh2-k4-serialized-cost-r1.R23/science'), '0bf22c6c4d2dfc84f74ccc7707eb0fb6d50c80e993dac568d7eed4559d28e9b4'),
    ('R24', Path('/var/tmp/wh2-k4-serialized-cost-r1.R24/science'), 'be20c9fb660d488ec172899302421a473e3748b49411344382cc6f404954e09c'),
)
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
PAIRS = ((0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,1),(2,1),(3,4),(5,4))
SHAPES = ((2,2),(2,1),(64,64),(64,1),(1280,1280),(1280,1))
ARMS = ('ordinary_independent','k4_independent','wh1_owned','ordinary_borrowed','k4_borrowed','wh1_borrowed')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def decode(raw):
    def unique(pairs):
        result = {}
        for key,value in pairs:
            require(key not in result,'duplicate JSON key')
            result[key] = value
        return result
    def invalid_constant(value):
        raise ValueError('nonfinite JSON constant: '+value)
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid_constant)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


def confidence(values):
    require(len(values) == 12 and all(math.isfinite(v) for v in values),'12 finite replicate logs')
    mean = math.fsum(values)/12
    radius = 2.200985160082949*math.sqrt(math.fsum((v-mean)**2 for v in values)/11/12)
    return dict(replicates=12,mean_log=mean,lower95_log=mean-radius,upper95_log=mean+radius,
                ratio=math.exp(mean),lower95=math.exp(mean-radius),upper95=math.exp(mean+radius))


def compare(a, b):
    if type(a) is float:
        require(type(b) is float and math.isfinite(b) and abs(a-b) < 1e-12,'finite statistical agreement')
    elif type(a) is dict:
        require(type(b) is dict and set(a) == set(b),'exact statistical schema')
        for key in a: compare(a[key],b[key])
    elif type(a) is list:
        require(type(b) is list and len(a) == len(b),'exact list extent')
        for x,y in zip(a,b): compare(x,y)
    else:
        require(type(a) is type(b) and a == b,'exact retained value')


def coordinates():
    index = 0
    for rep in range(12):
        for sweep in range(2):
            order = (rep+sweep)%2
            for width_slot in range(6):
                shape = (rep+sweep+width_slot)%6
                for metric_slot in range(3):
                    metric = (rep+sweep+width_slot+metric_slot)%3
                    for pair_slot in range(10):
                        pair = (2*rep+sweep+width_slot+metric+pair_slot)%10
                        for pos in range(18):
                            phase = rep+12*(rep%4) if pos < 2 else (rep+6*((pos-2)//8))%12+12*(((pos-2)%8)//2)
                            yield [index,rep,order,shape,metric,pair,pos,PAIRS[pair][SIDES[pos]^order],(2*phase+1)*1000000//96]
                            index += 1


def reconstruct(path, claim_sha, recorded):
    groups, panel, count, total = {}, [], 0, 0
    with path.open() as stream:
        header = decode(next(stream))
        require(header['type'] == 'header' and header['protocol'] == PROTOCOL and
                header['claim'] == claim_sha and header['batch'] == 128,'raw header binding')
        require(len(header['fixtures']) == 6,'complete fixture roster')
        previous = 0
        for expected in coordinates():
            row = decode(next(stream))
            require(row['type'] == 'record' and row['coordinate'] == expected and
                    row['complete'] is True and row['checked'] is True,'exact complete raw chronology')
            clocks = row['observation']['clocks']
            require(len(clocks) == 6 and all(type(c) is int and c >= 0 for c in clocks) and
                    previous <= clocks[0] <= clocks[2] < clocks[3] <= clocks[5] and
                    0 <= clocks[4]-clocks[1] <= clocks[5]-clocks[0],'ordered raw clocks')
            previous = clocks[5]
            _,_,order,shape,metric,pair,pos,arm,_ = expected
            steps = header['fixtures'][shape]['arms'][arm]['steps']
            ledger = [128,2560,0,0,0,128] if metric == 0 else [0,0,128,128*steps[metric-1],128,128]
            require(row['counts'] == ledger and row['address_count'] == len(row['addresses']) == 128,
                    'complete own-endpoint API ledger')
            require(all(type(a) is int and 0 < a < 2**64 for a in row['addresses']),'valid recorded handles')
            count += 1
            duration = clocks[3]-clocks[2]
            total += duration
            panel.append(duration)
            if pos == 17:
                ratios = []
                for p in range(2,18,2):
                    side = SIDES[p]^order
                    ratio = math.log(panel[p])-math.log(panel[p+1])
                    ratios.append(ratio if side else -ratio)
                groups.setdefault((shape,metric,pair,order),[]).append(math.fsum(ratios)/8)
                panel = []
        footer = decode(next(stream))
        require(footer == dict(type='footer',complete=True,records=count,work_ns=total) and
                next(stream,None) is None and count == 77760 and total <= 180_000_000_000,'exact terminal footer and WORK cap')
    stats, controls, treatments = [], [], []
    for key,logs in sorted(groups.items()):
        shape,metric,pair,order = key
        ci = confidence(logs)
        item = dict(width=SHAPES[shape][0],tail=SHAPES[shape][1],metric=metric,comparison=pair,order=order,
                    estimate=ci,replicate_logs=logs,comparison_arms=[ARMS[a] for a in PAIRS[pair]])
        if pair < 6:
            passed = -math.log1p(.02) < ci['lower95_log'] and ci['upper95_log'] < math.log1p(.02)
            item['control_pass'] = passed
            if not passed: controls.append(list(key))
        else:
            passed = ci['upper95_log'] < 0
            item.update(treatment_pass=passed,upper_ratio_limit=1.0)
            if not passed: treatments.append(list(key))
        stats.append(item)
    require(len(stats) == 360,'complete separate cell statistics')
    compare(stats,recorded['statistics'])
    outcome = 'CONTROL_FAIL' if controls else 'FAIL' if treatments else 'PASS'
    require((outcome,controls,treatments) == (recorded['outcome'],recorded['failed_controls'],recorded['failed_treatments']),
            'unchanged original decision, no rescue')
    return dict(records=count,warmups=count//9,work_ns=total,cell_intervals=len(stats),
                failed_controls=controls,failed_treatments=treatments,
                fixtures_sha256=hashlib.sha256(canonical(header['fixtures'])).hexdigest())


def bundle(label, root, sealed_sha):
    require(digest(root/'COMPLETE.json') == sealed_sha,'fixed historical completion hash: '+label)
    seal = decode((root/'COMPLETE.json').read_bytes())
    require(set(seal) == {'files','outcome','protocol'} and seal['protocol'] == PROTOCOL,'sealed protocol')
    names = {'CLAIM.json','analysis.json','raw.jsonl','stderr.txt'}
    require({p.name for p in root.iterdir()} == names|{'COMPLETE.json'} and len(seal['files']) == 4,'exact bundle extent')
    paths = set()
    for record in seal['files']:
        path = Path(record['path'])
        require(set(record) == {'path','bytes','sha256'} and path.parent == root and path.name in names and
                path not in paths,'unique exact sealed member')
        paths.add(path)
        require(path.stat().st_size == record['bytes'] and digest(path) == record['sha256'],'sealed member hash/size')
    claim = decode((root/'CLAIM.json').read_bytes())
    analysis = decode((root/'analysis.json').read_bytes())
    require(claim['protocol'] == analysis['protocol'] == PROTOCOL and analysis['outcome'] == seal['outcome'],'consistent protocol/outcome')
    result = dict(namespace=label,bundle=str(root),head=claim['head'],complete_sha256=sealed_sha,
                  recorded_outcome=analysis['outcome'],speed_qualification_claimed_here=False)
    if analysis['outcome'] == 'INVALID':
        require(label == 'R17' and (root/'raw.jsonl').stat().st_size == 0 and
                (root/'stderr.txt').read_bytes() == b'INVALID: claim bytes\n','exact pre-codec claim failure')
        result.update(records=0,retained_failure=analysis['failure'])
    else:
        require((root/'stderr.txt').stat().st_size == 0,'empty successful-worker stderr')
        result.update(reconstruct(root/'raw.jsonl',digest(root/'CLAIM.json'),analysis))
    for record in seal['files']:
        path = Path(record['path'])
        require(path.stat().st_size == record['bytes'] and digest(path) == record['sha256'],'bundle unchanged during traversal')
    require(digest(root/'COMPLETE.json') == sealed_sha,'completion unchanged during traversal')
    return result


def path_only_change(before, after):
    changed = subprocess.check_output(['git','diff','--name-only',before,after],cwd=ROOT,text=True).splitlines()
    expected = ['bench/Wh2K4FreshCostBuildR1.py','bench/Wh2K4SerializedCostR1.cpp']
    require(changed == expected or changed == ['.beads/interactions.jsonl']+expected,
            'only claim-path code changes and optional issue metadata')
    for name in expected:
        old = subprocess.check_output(['git','show',before+':'+name],cwd=ROOT)
        new = subprocess.check_output(['git','show',after+':'+name],cwd=ROOT)
        old_id = b'R22' if before.startswith('a802b66') else b'R23'
        new_id = b'R23' if after.startswith('6b2851e') else b'R24'
        require(old.count(old_id) == new.count(new_id) == 1 and old.replace(old_id,new_id) == new,'claim-path-only source change')
    return dict(before=before,after=after,changed_files=changed,work_or_codec_change=False)


def audit():
    results = [bundle(*entry) for entry in BUNDLES]
    changes = [path_only_change(results[i]['head'],results[i+1]['head']) for i in (3,4)]
    require(len({r['fixtures_sha256'] for r in results if r['records']}) == 1,'identical retained workload fixture bytes')
    return dict(scope='retained hashes/coordinates/WORK clocks/API counts/statistics and source history; not phase/wait/rusage/CPU-identity/payload/producer qualification',
                namespaces=results,path_only_relaunches=changes,records=sum(r['records'] for r in results),
                new_codec_executions=0,production_promotion_claimed=False,
                conclusion='R23 CONTROL_FAIL followed by claim-path-only R24 PASS is a prohibited retry, not promotion evidence')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True,help='New external output file, never a historical bundle member')
    args = parser.parse_args()
    report = args.report.resolve()
    require(ROOT not in report.parents and all(root not in report.parents for _,root,_ in BUNDLES),'external non-evidence report')
    require(not report.exists(),'never overwrite a report')
    result = audit()
    with report.open('xb') as out: out.write(canonical(result))
    print(json.dumps(dict(records=result['records'],report=str(report),sha256=digest(report),conclusion=result['conclusion']),sort_keys=True))
