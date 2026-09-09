#!/usr/bin/env python3
"""Frozen dimension-two noncommuting GF256 screen, never native timing."""
import importlib.util
import itertools
from pathlib import Path
import resource
import sys


def sibling(name,filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


G=sibling('k2_generic_selection','Wh2K5ThueMorseR0.py')
R=sibling('k2_frozen_input_controller','Wh2K2ThueMorseRunR0.py')
K,F,require=G.K,G.F,G.require
PROTOCOL=R.C.PROTOCOL
WIDTHS,SCHEDULES,WORDS=G.WIDTHS,G.SCHEDULES,G.WORDS
MINORS=tuple(itertools.combinations(range(6),2))
LEGACY_PAIRS=((265,270),(2,966),(1056,1313))
STRIDES=(255,257,65537)


def fixed_feedback():
    # (x+1)(x+2) over the existing byte field, no extension arithmetic.
    F.init_field()
    return F.MUL[1][2],1^2


def checked_rank(rows):
    actual=F.matrix_rank(rows)
    require(actual==G.reference_rank(rows),'independent GF256 rank disagreement')
    return actual


def checked_pair(rows):
    require(len(rows)==2 and all(len(row)==2 for row in rows),'two-row determinant shape')
    rank=checked_rank(rows); a,b=rows
    determinant=F.MUL[a[0]][b[1]]^F.MUL[a[1]][b[0]]
    reference=K.multiply_polynomial(a[0],b[1])^K.multiply_polynomial(a[1],b[0])
    require(determinant==reference and (determinant!=0)==(rank==2),'independent pair determinant')
    return dict(rank=rank,determinant=determinant)


def trace(b,root,schedule):
    require(b in WIDTHS and K.H.root_value(root) and schedule in SCHEDULES,'frozen trace coordinates')
    state=(int(root,16)^2*0x9e3779b97f4a7c15^b*0xbf58476d1ce4e5b9)&K.MASK64
    if schedule!='iid': state^=0x10fade
    def uniform():
        nonlocal state
        state=(state+0x9e3779b97f4a7c15)&K.MASK64
        v=((state^(state>>30))*0xbf58476d1ce4e5b9)&K.MASK64
        v=((v^(v>>27))*0x94d049bb133111eb)&K.MASK64
        return ((v^(v>>31))>>11)*2.0**-53
    ids=[]; burst=0; loss=.1 if schedule=='iid' else .5
    for candidate in range(67072):
        if schedule=='burst' and burst: burst-=1; continue
        if uniform()<(loss/(8-7*loss) if schedule=='burst' else loss):
            if schedule=='burst': burst=7
            continue
        ids.append(K.MAX_ID-2*candidate if schedule=='adversarial' else
            2+candidate if schedule=='repair-only' else candidate)
        if len(ids)==6: return ids
    raise F.ScreenInvalid('trace candidate cap')


def trace_result(b,root,schedule,mapper):
    ids=trace(b,root,schedule); rows=[mapper.row(i) for i in ids]
    return dict(B=b,root=root,schedule=schedule,ids=ids,ranks=[checked_rank(rows[:n]) for n in range(2,7)])


def exclusions(inputs):
    roots=list(inputs['roots'])+list(K.H.MAIN_ROOTS)
    for protocol in ('wirehair.wh2.k3-thue-morse-r0','wirehair.wh2.k5-thue-morse-r0',
                     'wirehair.wh2.thue-morse-recovery-r0'):
        roots.extend('0x'+F.digest((protocol+':fresh/'+str(i)).encode())[:16] for i in range(512))
    return sorted(set(roots))


def fresh_roots(excluded):
    roots=['0x'+F.digest((PROTOCOL+':fresh/'+str(i)).encode())[:16] for i in range(512)]
    require(len(set(roots))==512 and not set(roots)&set(excluded),'fresh root collision')
    return roots


def stride_pairs(stride,count=512):
    require(stride in STRIDES and type(count) is int and 0<count<=512,'stride domain')
    pairs=[]
    for i in range(count):
        a=int(F.digest((PROTOCOL+':stride/'+str(stride)+'/'+str(i)).encode())[:16],16)%(2**32-stride)
        pairs.append((a,a+stride))
    require(len(set(pairs))==count,'stride anchor collision, no replacement')
    return pairs


def summarize_fresh(rows,per_cell=512):
    require(len(rows)==per_cell*12,'complete fresh denominator')
    cells=[]
    for b,schedule in itertools.product(WIDTHS,SCHEDULES):
        cell=[r for r in rows if r['B']==b and r['schedule']==schedule]
        require(len(cell)==per_cell and len({r['root'] for r in cell})==per_cell,'complete unique cell')
        for r in cell:
            require(len(r['ranks'])==5 and all(type(v) is int and 0<=v<=2 for v in r['ranks']) and
                all(0<=y-x<=1 for x,y in zip(r['ranks'],r['ranks'][1:])),'nested K2 ranks')
        failures=[sum(r['ranks'][h]<2 for r in cell) for h in range(5)]
        cells.append(dict(B=b,schedule=schedule,traces=per_cell,failures=failures))
    totals=[sum(c['failures'][h] for c in cells) for h in range(5)]
    return dict(cells=cells,failures=totals,fresh_pass=totals[0]*100<=len(rows) and
        all(c['failures'][0]*100<=per_cell for c in cells))


def run_screen(claim):
    budget=F.Budget()
    result=dict(protocol=PROTOCOL,claim_sha256=claim,outcome='INVALID',selection=[],pair=None,
        inputs=None,local=[],seams=[],hard=[],history=[],legacy=[],fresh=[],strides=[],evidence={},summary={})
    try:
        result['evidence']['field']=F.init_field()
        feedback=fixed_feedback(); result['feedback']=feedback
        pair=G.choose_pair(feedback,WORDS,MINORS,budget,result['selection'])
        if pair is None: result['outcome']='EXHAUSTED'; return result
        result['pair']=pair
        require(F.matrix_multiply(pair[0],pair[1])!=F.matrix_multiply(pair[1],pair[0]),'noncommuting pair')
        require(all(checked_rank(m)==2 for m in pair),'invertible pair')
        for word in WORDS:
            columns=K.local_columns(pair,word)
            require(all(checked_rank([columns[i] for i in s])==2 for s in MINORS),'selected local certificate')
            result['local'].append(dict(word=word,columns=columns,checked=15))
        mapper=K.Mapper(pair,budget)
        require(len(mapper.payload)==7168 and tuple(mapper.row(i) for i in range(2))==K.identity(2),'systematic packed K2 geometry')
        product=K.identity(2)
        for i in range(2049):
            budget.check(); require(mapper.row(i)==tuple(row[0] for row in product),'literal sequential oracle')
            product=F.matrix_multiply(product,pair[K.parity(i)])
        result['evidence'].update(lookup_bytes=7168,lookup_sha256=F.digest(mapper.payload),
            pair_sha256=F.digest(bytes(v for matrix in pair for row in matrix for v in row)))
        for start in [(1<<e)-4 for e in range(3,32)]+[K.MAX_ID-5]:
            budget.check(); ids=list(range(start,start+6)); rows=[mapper.row(i) for i in ids]
            result['seams'].append(dict(ids=ids,deficient=[list(s) for s in MINORS if checked_rank([rows[i] for i in s])<2]))
        result['inputs']=R.inventory_inputs(budget.deadline)
        for root,b,schedule in itertools.product(K.H.HARD_TRAINING_ROOTS+K.H.HARD_VALIDATION_ROOTS,WIDTHS,SCHEDULES):
            budget.check(); result['hard'].append(trace_result(b,root,schedule,mapper))
        for p in result['inputs']['prefixes']:
            budget.check(); result['history'].append(dict(p,rank=checked_rank([mapper.row(i) for i in p['ids']])))
        for ids in LEGACY_PAIRS:
            result['legacy'].append(dict(ids=ids,rank=checked_rank([mapper.row(i) for i in ids])))
        structural=(all(not r['deficient'] for r in result['seams']) and
            all(r['ranks'][0]==2 for r in result['hard']) and
            all(r['rank']==2 for r in result['history']+result['legacy']))
        result['summary']=dict(structural_pass=structural,fresh_entered=structural)
        if structural:
            roots=fresh_roots(exclusions(result['inputs']))
            for root,b,schedule in itertools.product(roots,WIDTHS,SCHEDULES):
                budget.check(); result['fresh'].append(trace_result(b,root,schedule,mapper))
            result['summary'].update(summarize_fresh(result['fresh']))
            for stride in STRIDES:
                pairs=[]
                for ids in stride_pairs(stride):
                    budget.check(); pairs.append(dict(ids=ids,**checked_pair([mapper.row(i) for i in ids])))
                failures=sum(r['rank']<2 for r in pairs)
                result['strides'].append(dict(stride=stride,pairs=pairs,failures=failures,passed=failures*100<=512))
            result['summary']['stride_pass']=all(r['passed'] for r in result['strides'])
        result['evidence']['unique_rows']=[dict(id=i,row=mapper.cache[i]) for i in sorted(mapper.cache)]
        result['counts']=dict(local_minors=150,seam_minors=450,hard_traces=len(result['hard']),
            history_prefixes=len(result['history']),legacy_pairs=len(result['legacy']),fresh_traces=len(result['fresh']),
            stride_pairs=sum(len(s['pairs']) for s in result['strides']),unique_rows=len(mapper.cache))
        require(len(result['local'])==10 and len(result['seams'])==30 and len(result['hard'])==72 and
            len(result['history'])==13 and len(result['legacy'])==3 and
            len(result['fresh'])==(6144 if structural else 0) and
            len(result['strides'])==(3 if structural else 0),'complete screen accounting')
        budget.check()
        result['outcome']='PASS' if structural and result['summary']['fresh_pass'] and result['summary']['stride_pass'] else 'FAIL'
    except Exception as error:
        result.update(outcome='INVALID',error=(type(error).__name__+': '+str(error))[:1024])
    return result


def main(argv):
    require(argv==['--worker'],'usage: Wh2K2ThueMorseR0.py --worker')
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    claim=R.claimed_inputs(); result=run_screen(claim); raw=F.canonical(result)+b'\n'
    require(len(raw)<=F.STDOUT_LIMIT,'worker output cap')
    sys.stdout.buffer.write(raw); sys.stdout.buffer.flush()
    return int(result['outcome']=='INVALID')


if __name__=='__main__': sys.exit(main(sys.argv[1:]))
