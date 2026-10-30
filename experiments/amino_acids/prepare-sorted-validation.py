"""Recover the experimental input contract using existing production-sort outputs; preserve all prior evidence."""
import csv,gzip,hashlib,importlib.util,json,shutil
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-cohort';archive=d/'failed-pooled-order-20261022T073824';archive.mkdir(exist_ok=True)
assert not (d/'complete.json').exists() and not (d/'manifest-sorted.tsv').exists()
log=(d/'nextflow.log').read_text();assert 'File 1 is not sorted' in log and 'Session UUID: 98a3b25c-949c-4eac-adee-c0ef24cb3f78' in log
for name,path in [('nextflow.log',d/'nextflow.log'),('run.log',d/'run.log'),('trace.tsv',d/'output/trace.tsv'),('supervisor.log',r/'cohort-validation-resumed.log')]:shutil.copy2(path,archive/name)
trace=list(csv.DictReader((archive/'trace.tsv').open(),delimiter='\t'));arms=['nt','nt22','nt20','aa12','union'];audit=[];sorted_dir=d/'sorted-inputs';sorted_dir.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('sampler',r/'modules/local/downsampleTsvByHash/resources/usr/bin/downsample_tsv_by_hash.py');sampler=importlib.util.module_from_spec(spec);spec.loader.exec_module(sampler)
# The surviving production invariant: selected keys are independent of input order.
keys=['SRR38294902.84354065','SRR33038013.1363664172']+[f'read_{i}' for i in range(200)]
a,_=sampler.select_indices(iter(keys),20);rev=list(reversed(keys));b,_=sampler.select_indices(iter(rev),20);assert {keys[i] for i in a}=={rev[i] for i in b}
for arm in arms:
 tasks=[x for x in trace if x['name']==f'VALIDATE_VIRAL_ASSIGNMENTS:SPLIT_VIRAL_TSV_BY_SELECTED_TAXID:SORT_SEQ_ID (id={arm})'];assert len(tasks)==1 and tasks[0]['status'] in ['COMPLETED','CACHED']
 work=list((d/'work').glob(tasks[0]['hash']+'*'));assert len(work)==1
 source=work[0]/f'sorted_seq_id_{arm}-pooled.tsv.gz';assert source.exists()
 original=d/f'{arm}-pooled.tsv.gz';rows={};old_order_errors=0;last=None
 with gzip.open(original,'rb') as f:
  header=f.readline();column=header.rstrip(b'\n').split(b'\t').index(b'seq_id')
  for line in f:
   key=line.rstrip(b'\n').split(b'\t')[column];assert key not in rows;rows[key]=line
   if last is not None and key<last:old_order_errors+=1
   last=key
 count=len(rows);last=None;canonical=hashlib.sha256();canonical.update(header)
 with gzip.open(source,'rb') as f:
  assert f.readline()==header
  for line in f:
   key=line.rstrip(b'\n').split(b'\t')[column];assert last is None or last<key;last=key
   assert rows.pop(key)==line;canonical.update(line)
 assert not rows
 target=sorted_dir/f'{arm}-pooled.tsv.gz';shutil.copy2(source,target)
 prior=[x for x in trace if x['name']==f'VALIDATE_VIRAL_ASSIGNMENTS:BLAST_FASTA:BLAST (id={arm})'];assert len(prior)==1 and prior[0]['status']=='COMPLETED' and prior[0]['exit']=='0'
 bw=list((d/'work').glob(prior[0]['hash']+'*'));assert len(bw)==1
 old=archive/arm;old.mkdir()
 for name in [f'{arm}_sampled_reads.fasta.gz',f'{arm}_hits.blast.gz','.command.sh','.command.log','.command.trace']:shutil.copy2(bw[0]/name,old/name)
 audit.append({'arm':arm,'rows':count,'original_order_inversions':old_order_errors,'sorted_strictly_by_unique_seq_id':True,'all_row_bytes_preserved':True,'canonical_uncompressed_sha256':canonical.hexdigest(),'source_sort_task':tasks[0],'source_BLAST_task':prior[0],'sorted_file':str(target)})
 print(json.dumps({k:v for k,v in audit[-1].items() if k not in ['source_sort_task','source_BLAST_task']}),flush=True)
(d/'sorting-recovery-audit.json').write_text(json.dumps({'production_hash_selection_order_invariance_fixture':True,'arms':audit,'production_modules_unchanged':True,'contract':'Complete unique read rows and duplicate labels preserved byte-for-byte; global seq_id ordering restored. Subsequent selected-query and alignment parity is measured separately.'},indent=2))
with (d/'manifest-sorted.tsv').open('w') as f:
 f.write('id\tpath\n')
 for x in audit:f.write(x['arm']+'\t'+x['sorted_file']+'\n')
print('SORTED_INPUTS_VERIFIED',flush=True)
