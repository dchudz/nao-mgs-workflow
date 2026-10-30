"""Restore the experimental sorted-input boundary and measure query/alignment parity after the original workflow."""
import collections,csv,datetime,gzip,hashlib,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-cohort';session='98a3b25c-949c-4eac-adee-c0ef24cb3f78';archive=d/'failed-pooled-order-20261022T073824'
assert not (d/'complete.json').exists()
audit=json.loads((d/'sorting-recovery-audit.json').read_text());assert audit['production_hash_selection_order_invariance_fixture'] and len(audit['arms'])==5 and all(x['all_row_bytes_preserved'] and x['sorted_strictly_by_unique_seq_id'] for x in audit['arms'])
scope=json.loads((d/'scope.json').read_text());assert len(scope['completed_samples'])==1000 and not scope['missing_samples']
cmd=[str(r/'nextflow'),'-log',str(d/'nextflow.log'),'run','experiments/amino_acids/validation-pooled.nf','-resume',session,'-c','experiments/amino_acids/validation-cohort.config','-work-dir',str(d/'work'),'--manifest',str(d/'manifest-sorted.tsv'),'--ref_dir',str(r/'reference'),'--output',str(d/'output')]
(d/'sorting-recovery-protocol.json').write_text(json.dumps({'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'session':session,'production_modules_unchanged':True,'input_correction':'Use existing globally seq_id-sorted production outputs; all original read rows and duplicate labels verified byte-for-byte.','surviving_invariants':['Unique original IDs and complete row bytes','Order-independent production hash selection','Same per-ID FASTA queries, checked before committing the result','Unchanged search and annotation parameters','No validation extrapolation between original reads']},indent=2))
with (d/'run.log').open('w') as f:subprocess.run(cmd,cwd=r,stdout=f,stderr=f,check=True,timeout=28800)
assert len(list((d/'output/COPY_HITS').glob('*_validation_hits.tsv.gz')))==5
trace=list(csv.DictReader((d/'output/trace.tsv').open(),delimiter='\t'));parity=[]
def fasta(path):
 content=gzip.open(path,'rb').read();records={}
 for block in content.split(b'>')[1:]:
  lines=block.splitlines();identifier=lines[0].split()[0];assert identifier not in records;records[identifier]=(lines[0],b''.join(lines[1:]))
 return records,hashlib.sha256(content).hexdigest()
def hits(path):
 with gzip.open(path,'rt') as f:return collections.Counter(line.rstrip('\n') for line in f)
for arm in ['nt','nt22','nt20','aa12','union']:
 rows=[x for x in trace if x['name']==f'VALIDATE_VIRAL_ASSIGNMENTS:BLAST_FASTA:BLAST (id={arm})'];assert len(rows)==1 and rows[0]['status'] in ['COMPLETED','CACHED'] and rows[0]['exit']=='0'
 work=list((d/'work').glob(rows[0]['hash']+'*'));assert len(work)==1
 before,before_sha=fasta(archive/arm/f'{arm}_sampled_reads.fasta.gz');after,after_sha=fasta(work[0]/f'{arm}_sampled_reads.fasta.gz');assert before==after,(arm,'Selected per-ID query content changed')
 old=hits(archive/arm/f'{arm}_hits.blast.gz');new=hits(work[0]/f'{arm}_hits.blast.gz');removed=old-new;added=new-old
 parity.append({'arm':arm,'queries':len(after),'selected_query_IDs_headers_sequences_identical':True,'uncompressed_FASTA_bytes_identical':before_sha==after_sha,'old_query_sha256':before_sha,'new_query_sha256':after_sha,'alignment_multiset_identical':old==new,'old_alignment_rows':sum(old.values()),'new_alignment_rows':sum(new.values()),'removed_alignment_rows':sum(removed.values()),'added_alignment_rows':sum(added.values()),'differing_query_IDs':sorted({line.split('\t',1)[0] for line in removed.keys()|added.keys()})})
(d/'sorting-recovery-parity.json').write_text(json.dumps({'arms':parity,'all_selected_queries_preserved':True},indent=2))
with (d/'analysis.log').open('w') as f:subprocess.run(['python3',str(r/'summarize-native-validation.py'),str(d),'--nodes',str(r/'reference/results/taxonomy-nodes.dmp')],stdout=f,stderr=f,check=True)
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-cohort/'
aws=['aws','--profile','pipeline-codex','s3','cp'];subprocess.run(aws+[str(d)+'/',uri,'--recursive','--exclude','work/*','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_samples':1000,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resumed_session':session,'summary_complete':True,'sorted_input_contract_verified':True,'selected_queries_preserved':True},indent=2));subprocess.run(aws+[str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
print('COHORT_VALIDATION_COMPLETE',flush=True)
