"""Characterize original-read pooling using real fixture output and labeled constructed cap cases."""
import csv,importlib.util,json,tempfile
from pathlib import Path
spec=importlib.util.spec_from_file_location('runner',Path(__file__).with_name('run-query-only.py'));r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
f=r.m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=f/'query-only';records=r.m.fasta(f/'queries.fna');mapping=json.loads((d/'query-map.json').read_text())
rows,counts=r.normalize_pool(d/'tblastn-fast.raw.tsv',mapping,records)
assert counts['original_queries']==len(records) and counts['normalized_HSPs']==len(rows)
assert all(0<float(x['qcovs'])<=100 and x['qseqid'] in records for x in rows)
raw=list(csv.reader((d/'tblastn-fast.raw.tsv').open(),delimiter='\t'));source=next(x for x in raw if int(x[r.m.FIELDS.index('length')])==80)
# Repeated identical HSP coordinates must not inflate aggregate coverage; caps apply after frame pooling.
constructed=[]
for i in range(251):
    for j in range(6):
        x=dict(zip(r.m.FIELDS,source));x['sseqid']='constructed_subject_'+str(i);x['bitscore']=str(1000-i-j/100);constructed.append(x)
with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)/'constructed.tsv'
    with p.open('w') as out:csv.DictWriter(out,fieldnames=r.m.FIELDS,delimiter='\t').writerows(constructed)
    pooled,c=r.normalize_pool(p,mapping,records)
assert len(pooled)==250*5 and c['dropped_subjects_after_pool_cap']==1
assert len({x['sseqid'] for x in pooled})==250 and 'constructed_subject_250' not in {x['sseqid'] for x in pooled}
assert all(abs(float(x['qcovs'])-100*(abs(int(x['qend'])-int(x['qstart']))+1)/int(x['qlen']))<1e-9 for x in pooled)
result=dict(created_utc=r.m.now(),source_sha256=r.m.sha(__file__),runner_sha256=r.m.sha(Path(r.__file__)),raw_fixture_sha256=r.m.sha(d/'tblastn-fast.raw.tsv'),fixture_counts=counts,constructed_case_counts=c,invariants=['Every normalized HSP uses original read identifiers and valid nucleotide coordinates','All original queries remain counted, including no-hit reads','Per-read subject cap250 and per-read/subject HSP cap5 apply after pooling','Duplicate HSP spans do not inflate nucleotide aggregate coverage'],limits=['Constructed cap cases test normalization only; they are not biological search results','Search heuristic and score parity with tblastx are not asserted'])
r.m.save(d/'pooling-characterization.json',result);print(json.dumps(result))
