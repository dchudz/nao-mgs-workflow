"""Compare the completed pilot nucleotide rows with its already searched parent pool."""
import collections,csv,importlib.util,json
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r=m.ROOT/'cohort-storage';d=r/'translated-core-nt/pilot24';parent=r/'seed-diagnostic200';queries=m.fasta(d/'queries.fna');parents=m.fasta(parent/'queries.fna');assert all(parents[q]==s for q,s in queries.items())
fields=['qseqid','sseqid','staxids','bitscore','evalue','length','pident','qcovs'];result={}
for method in ['megablast','blastn']:
    assert (d/method/'complete.json').exists()
    old=collections.Counter(tuple(row) for row in csv.reader((parent/(method+'.tsv')).open(),delimiter='\t') if row[0] in queries)
    new=collections.Counter()
    for row in csv.reader((d/method/'alignments.tsv').open(),delimiter='\t'):
        x=dict(zip(m.FIELDS,row));new[tuple(x[k] for k in fields)]+=1
    removed=old-new;added=new-old
    result[method]=dict(parent_rows=sum(old.values()),pilot_rows=sum(new.values()),removed_rows=sum(removed.values()),added_rows=sum(added.values()),changed_queries=sorted({x[0] for x in set(removed)|set(added)}),row_multiset_equal=old==new)
    with (d/(method+'-parent-context-differences.tsv')).open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['direction','multiplicity']+fields)
        for direction,data in [('parent_only',removed),('pilot_only',added)]:
            for row,count in sorted(data.items()):w.writerow([direction,count,*row])
m.save(d/'parent-context-comparison.json',dict(created_utc=m.now(),methods=result,pilot_queries=len(queries),parent_queries=len(parents),note='Exact input sequence parity; common8outputfields compared as multisets. Prior parent search used8threads and pilot16; query pool and thread count both differ, so any difference is not attributable solely to batching.'))
print(json.dumps(result))
