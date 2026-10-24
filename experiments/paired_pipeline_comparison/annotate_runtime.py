"""Identify concurrent storage work without subtracting it from stage wall time."""
import argparse
import csv
import datetime
import json
from pathlib import Path


def seconds(text):
    return datetime.datetime.fromisoformat(text).timestamp()


def overlap(intervals, start, end):
    total=0;last=start
    for lo,hi in sorted(intervals):
        lo=max(lo,start,last);hi=min(hi,end)
        if hi>lo:total+=hi-lo;last=hi
    return total


def main():
    p=argparse.ArgumentParser()
    for key in ('results','provenance','selection','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();selection=json.loads(a.selection.read_text())
    workers={selection['first_run']:'primary'}
    for shard in a.provenance.glob('*-shard.json'):
        x=json.loads(shard.read_text())
        for sample in x['runs']:
            assert sample not in workers;workers[sample]=x['name']
    intervals={}
    for pattern in ('*-reference-storage-check.json','*-reference-read-initialization.json','*-alignment-index-storage-check.json','*-validation-cache-storage-check.json'):
        for path in a.provenance.glob(pattern):
            x=json.loads(path.read_text());worker=x['worker'].replace('paired-study','primary').replace('paired-','')
            intervals.setdefault(worker,[]).append((seconds(x['started_utc']),seconds(x['finished_utc'])))
    rows=[]
    for path in sorted(a.results.glob('*/comparison-report/stage-wall-times.tsv')):
        label=path.parent.parent.name;sample=label.rsplit('_',1)[0];worker=workers.get(sample)
        if worker is None:raise ValueError('Missing worker provenance for '+label)
        with path.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                amount=overlap(intervals.get(worker,[]),seconds(row['started_utc']),seconds(row['finished_utc']))
                rows.append({'unit':label,'worker':worker,'worker_origin':'retained' if worker=='primary' else 'snapshot replica',
                             **row,'concurrent_storage_seconds':amount,'concurrent_storage_work':amount>0})
    a.out.mkdir(parents=True,exist_ok=True)
    if rows:
        with (a.out/'runtime-context.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (a.out/'runtime-context.json').write_text(json.dumps({'stages':len(rows),'stages_with_concurrent_storage_work':sum(r['concurrent_storage_work'] for r in rows),
        'interpretation':'Overlap is a union of recorded storage intervals. It does not establish slowdown and must not be subtracted from stage time. Stages without a recorded overlap are not proven warm-cache or steady-state measurements.'},indent=2)+'\n')


if __name__=='__main__':
    assert overlap([(0,10),(5,15)],2,12)==10
    assert overlap([(0,10),(2,5)],0,20)==10
    assert overlap([(0,2),(4,8)],1,6)==3
    assert overlap([(0,1)],2,4)==0
    main()
