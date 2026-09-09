#!/usr/bin/env python3
"""Frozen finite K5 GF256 structural screen. No native codec or timing claim."""
import importlib.util
import itertools
from pathlib import Path
import resource
import sys


def sibling(name, filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


K = sibling('k5_generic_mapper','Wh2K3ThueMorseR0.py')
R = sibling('k5_inventory_reader','Wh2K5ThueMorseRunR0.py')
F, require = K.F, K.require
PROTOCOL = R.C.PROTOCOL
WIDTHS, SCHEDULES, WORDS = K.WIDTHS, K.SCHEDULES, K.WORDS
MINORS = tuple(itertools.combinations(range(9),5))
REFERENCE = None


def fixed_feedback():
    F.init_field()
    polynomial=[1]
    for root in F.EXP[:5]:
        out=[0]*(len(polynomial)+1)
        for i,value in enumerate(polynomial):
            out[i]^=F.MUL[root][value];out[i+1]^=value
        polynomial=out
    require(len(polynomial)==6 and polynomial[-1]==1 and polynomial[0]!=0,'monic invertible feedback')
    return tuple(polynomial[:-1])


def reference_rank(rows):
    """Opposite-column full elimination with independently built field tables."""
    global REFERENCE
    if REFERENCE is None:
        tables=tuple(bytes(K.multiply_polynomial(a,b) for b in range(256)) for a in range(256))
        inverses=(0,)+tuple(tables[a].index(1) for a in range(1,256))
        REFERENCE=tables,inverses
    tables,inverses=REFERENCE
    work=[list(row) for row in rows]
    if not work:return 0
    width=len(work[0])
    require(all(len(row)==width and all(type(v) is int and 0<=v<256 for v in row) for row in work),'rank input')
    pivot=0
    for column in range(width-1,-1,-1):
        chosen=next((i for i in range(len(work)-1,pivot-1,-1) if work[i][column]),None)
        if chosen is None:continue
        work[pivot],work[chosen]=work[chosen],work[pivot]
        work[pivot]=[tables[inverses[work[pivot][column]]][v] for v in work[pivot]]
        for i in range(len(work)):
            if i!=pivot:
                table=tables[work[i][column]]
                work[i]=[v^table[u] for v,u in zip(work[i],work[pivot])]
        pivot+=1
        if pivot==len(work):break
    return pivot


def checked_rank(rows):
    actual=F.matrix_rank(rows)
    require(actual==reference_rank(rows),'independent rank disagreement')
    return actual


def choose_pair(feedback,words,minors,budget,records):
    n=len(feedback)
    require(n>=2 and all(type(v) is int and 0<=v<256 for v in feedback) and feedback[0]!=0,'feedback')
    for value in range(1,256):
        if value==feedback[0]:continue
        budget.check()
        pair=(K.companion(feedback),K.companion((feedback[0]^value,)+tuple(feedback[1:])))
        entry=dict(parameter=value,checked=0,first_failure=None);records.append(entry)
        for word in words:
            columns=K.local_columns(pair,word)
            for selected in minors:
                entry['checked']+=1
                if F.matrix_rank([columns[i] for i in selected])!=n:
                    entry['first_failure']=dict(word=word,columns=list(selected))
                    break
            if entry['first_failure'] is not None:break
        if entry['first_failure'] is None:return pair
    return None


def trace(b,root,schedule):
    require(b in WIDTHS and K.H.root_value(root) and schedule in SCHEDULES,'trace coordinates')
    state=(int(root,16)^5*0x9e3779b97f4a7c15^b*0xbf58476d1ce4e5b9)&K.MASK64
    if schedule!='iid':state^=0x10fade
    def uniform():
        nonlocal state
        state=(state+0x9e3779b97f4a7c15)&K.MASK64
        v=((state^(state>>30))*0xbf58476d1ce4e5b9)&K.MASK64
        v=((v^(v>>27))*0x94d049bb133111eb)&K.MASK64
        return ((v^(v>>31))>>11)*2.0**-53
    ids=[];burst=0;loss=.1 if schedule=='iid' else .5
    for candidate in range(67840):
        if schedule=='burst' and burst:
            burst-=1;continue
        if uniform()<(loss/(8-7*loss) if schedule=='burst' else loss):
            if schedule=='burst':burst=7
            continue
        ids.append(K.MAX_ID-2*candidate if schedule=='adversarial' else 5+candidate if schedule=='repair-only' else candidate)
        if len(ids)==9:return ids
    raise F.ScreenInvalid('trace candidate bound')


def trace_result(b,root,schedule,mapper):
    ids=trace(b,root,schedule);rows=[mapper.row(i) for i in ids]
    return dict(B=b,root=root,schedule=schedule,ids=ids,ranks=[checked_rank(rows[:n]) for n in range(5,10)])


def fresh_roots(excluded):
    roots=['0x'+F.digest((PROTOCOL+':fresh/'+str(i)).encode())[:16] for i in range(512)]
    require(len(set(roots))==512 and not set(roots)&set(excluded),'fresh root collision')
    return roots


def exclusions(inputs):
    values=list(inputs['roots'])+list(K.H.MAIN_ROOTS)
    for protocol in ('wirehair.wh2.k3-thue-morse-r0','wirehair.wh2.thue-morse-recovery-r0'):
        values.extend('0x'+F.digest((protocol+':fresh/'+str(i)).encode())[:16] for i in range(512))
    return sorted(set(values))


def summarize_fresh(rows,per_cell=512):
    require(len(rows)==per_cell*12,'fresh complete denominator')
    cells=[]
    for b,schedule in itertools.product(WIDTHS,SCHEDULES):
        cell=[r for r in rows if r['B']==b and r['schedule']==schedule]
        require(len(cell)==per_cell and len({r['root'] for r in cell})==per_cell,'fresh cell denominator')
        for r in cell:
            require(len(r['ranks'])==5 and all(type(v) is int and 0<=v<=5 for v in r['ranks']) and
                    all(0<=y-x<=1 for x,y in zip(r['ranks'],r['ranks'][1:])),'nested ranks')
        failures=[sum(r['ranks'][oh]<5 for r in cell) for oh in range(5)]
        cells.append(dict(B=b,schedule=schedule,traces=per_cell,failures=failures,
            first_success=[per_cell-failures[0]]+[failures[i-1]-failures[i] for i in range(1,5)]+[failures[4]]))
    totals=[sum(c['failures'][i] for c in cells) for i in range(5)]
    return dict(cells=cells,failures=totals,fresh_pass=totals[0]*100<=len(rows) and
                all(c['failures'][0]*100<=per_cell for c in cells))


def run_screen(claim):
    budget=F.Budget()
    result=dict(protocol=PROTOCOL,claim_sha256=claim,outcome='INVALID',selection=[],pair=None,
        inputs=None,local=[],seams=[],hard=[],history=[],fresh=[],evidence={},summary={})
    try:
        result['evidence']['field']=F.init_field()
        feedback=fixed_feedback();result['feedback']=feedback
        pair=choose_pair(feedback,WORDS,MINORS,budget,result['selection'])
        if pair is None:
            result['outcome']='EXHAUSTED';return result
        result['pair']=pair
        # Selected once using local algebra. History/loss outcomes never enter selection.
        for word in WORDS:
            columns=K.local_columns(pair,word)
            require(all(checked_rank([columns[i] for i in selected])==5 for selected in MINORS),'selected local certificate')
            result['local'].append(dict(word=word,columns=columns,checked=126))
        mapper=K.Mapper(pair,budget)
        require(len(mapper.payload)==29440,'lookup byte count')
        require(tuple(mapper.row(i) for i in range(5))==K.identity(5),'systematic rows')
        product=K.identity(5)
        for packet_id in range(2049):
            budget.check()
            require(mapper.row(packet_id)==tuple(row[0] for row in product),'literal sequential product oracle')
            product=F.matrix_multiply(product,pair[K.parity(packet_id)])
        result['evidence'].update(lookup_bytes=len(mapper.payload),lookup_sha256=F.digest(mapper.payload),
            pair_sha256=F.digest(bytes(v for m in pair for row in m for v in row)))
        for start in [(1<<e)-4 for e in range(3,32)]+[K.MAX_ID-8]:
            budget.check();ids=list(range(start,start+9));rows=[mapper.row(i) for i in ids]
            deficient=[list(s) for s in MINORS if checked_rank([rows[i] for i in s])<5]
            result['seams'].append(dict(ids=ids,deficient=deficient))
        result['inputs']=R.inventory_inputs(budget.deadline)
        for root,b,schedule in itertools.product(K.H.HARD_TRAINING_ROOTS+K.H.HARD_VALIDATION_ROOTS,WIDTHS,SCHEDULES):
            budget.check();result['hard'].append(trace_result(b,root,schedule,mapper))
        for ids in result['inputs']['prefixes']:
            budget.check();result['history'].append(dict(ids=ids,rank=checked_rank([mapper.row(i) for i in ids])))
        structural=all(not row['deficient'] for row in result['seams']) and all(r['ranks'][0]==5 for r in result['hard']) and all(r['rank']==5 for r in result['history'])
        result['summary']=dict(structural_pass=structural,fresh_entered=structural)
        if structural:
            roots=fresh_roots(exclusions(result['inputs']))
            for root,b,schedule in itertools.product(roots,WIDTHS,SCHEDULES):
                budget.check();result['fresh'].append(trace_result(b,root,schedule,mapper))
            result['summary'].update(summarize_fresh(result['fresh']))
        result['evidence']['unique_rows']=[dict(id=i,row=mapper.cache[i]) for i in sorted(mapper.cache)]
        result['counts']=dict(local_minors=1260,seam_minors=3780,hard_traces=len(result['hard']),
            history_prefixes=len(result['history']),fresh_traces=len(result['fresh']),unique_rows=len(mapper.cache))
        require(len(result['local'])==10 and len(result['seams'])==30 and len(result['hard'])==72 and
                len(result['history'])==54 and len(result['fresh'])==(6144 if structural else 0),'complete screen accounting')
        budget.check()
        result['outcome']='PASS' if structural and result['summary']['fresh_pass'] else 'FAIL'
    except Exception as error:
        result['error']=(type(error).__name__+': '+str(error))[:1024]
        result['outcome']='INVALID'
    return result


def main(argv):
    require(argv==['--worker'],'usage: Wh2K5ThueMorseR0.py --worker')
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    claim=R.claimed_inputs()
    result=run_screen(claim);raw=F.canonical(result)+b'\n'
    require(len(raw)<=F.STDOUT_LIMIT,'worker output cap')
    sys.stdout.buffer.write(raw);sys.stdout.buffer.flush()
    return int(result['outcome']=='INVALID')


if __name__=='__main__':sys.exit(main(sys.argv[1:]))
