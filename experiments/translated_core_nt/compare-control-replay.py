"""Compare unchanged-query nucleotide controls across the earlier shared8thread and current16thread runs."""
import collections,csv,hashlib,json,datetime
from pathlib import Path
r=Path('/home/ec2-user/experiment/cohort-storage');old=r/'protein-nt-panel';new=r/'translated-core-nt/protein232';assert hashlib.sha256((old/'queries.fna').read_bytes()).digest()==hashlib.sha256((new/'queries.fna').read_bytes()).digest()
result={}
for method in ['megablast','blastn']:
 assert (new/method/'complete.json').exists()
 left=collections.Counter(tuple(x) for x in csv.reader((old/(method+'.tsv')).open(),delimiter='\t'));right=collections.Counter(tuple(x[:8]) for x in csv.reader((new/method/'alignments.tsv').open(),delimiter='\t'))
 lo=left-right;ro=right-left
 result[method]={'old_rows':sum(left.values()),'new_rows':sum(right.values()),'alignment_multiset_equal_on_original8fields':left==right,'rows_only_old':sum(lo.values()),'rows_only_new':sum(ro.values()),'differing_query_mates':sorted({x[0] for x in lo.keys()|ro.keys()}),'old_sha256':hashlib.sha256((old/(method+'.tsv')).read_bytes()).hexdigest(),'new_sha256':hashlib.sha256((new/method/'alignments.tsv').read_bytes()).hexdigest()}
 if lo or ro:
  with (new/(method+'-control-replay-differences.tsv')).open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['side','count','qseqid','sseqid','staxids','bitscore','evalue','length','pident','qcovs'])
   for side,values in [('old',lo),('new',ro)]:
    for row,count in sorted(values.items()):w.writerow([side,count,*row])
out={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':result,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limits':['Same query bytes, frozen database and nucleotide task parameters. New runs use16threads and explicitE10; older8thread runs used pinned defaultE10.','Multiset comparison removes row-order differences and projects new19field output onto the exact original8fields.','Old wall times include shared-worker contention; time ratios are not algorithm speedups.']};(new/'control-replay-parity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
