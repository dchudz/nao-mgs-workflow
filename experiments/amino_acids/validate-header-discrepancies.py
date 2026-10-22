"""Search the eight observed filtering discrepancies against the unchanged comprehensive core_nt."""
import csv,datetime,hashlib,json,subprocess,sys,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
parent=r/'cohort-storage/preparation-discrepancies';d=parent/'nt-validation';d.mkdir()
selected={x['seq_id']:x for x in json.loads((parent/'selected.json').read_text())};assert len(selected)==8
with (d/'selection.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['query_pair','sample','stratum','viral_reference_taxid'],delimiter='\t');w.writeheader()
 for q,x in selected.items():w.writerow(dict(query_pair=q,sample=x['sample'],stratum='header_filter_discrepancy',viral_reference_taxid=28875))
with (d/'queries.fna').open('w') as f:
 for key,raw in pairs(parent/'input.fastq.gz'):
  q=key[0].decode().lstrip('@');assert q in selected
  for mate,seq in [(1,key[1]),(2,key[3])]:f.write(f'>{q}/{mate}\n{seq.decode()}\n')
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pairs':8,'queries':16,'query_sha256':hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),'selection':'Eight outcome-selected native RotavirusA assignments excluded by direct human filtering but retained after the native viral SAM-to-FASTQ header transformation.','purpose':'Compare the full nucleotide reference evidence for the same complete mate sequences; do not infer biological truth from the header-dependent human-filter decision.','reference':'Unchanged comprehensive core_nt; no reduced competitor set.','primary_results':'Unchanged. This is an additional mechanism diagnostic.'};(d/'protocol.json').write_text(json.dumps(protocol,indent=2))
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
for task in ['megablast','blastn']:
 name='aa-header-validation-'+task
 cmd=['docker','run','--rm','--name',name,'--network','none','--cpu-shares','2048','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{r}:{r}',image,'/usr/bin/time','-v','-o',str(d/f'{task}.time'),'blastn','-task',task,'-db',str(r/'scratch/blast_db_190b4ae1/blast_db'),'-query',str(d/'queries.fna'),'-out',str(d/f'{task}.tsv'),'-perc_identity','60','-max_hsps','5','-num_alignments','250','-qcov_hsp_perc','30','-num_threads','2','-outfmt','6 qseqid sseqid staxids bitscore evalue length pident qcovs']
 start=time.monotonic()
 try:
  with (d/f'{task}.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=7200)
 finally:subprocess.run(['docker','stop','-t','10',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 (d/f'{task}.resources.json').write_text(json.dumps({'wall_seconds':time.monotonic()-start,'command':cmd},indent=2));print(task+' complete',flush=True)
subprocess.run(['python3',str(r/'summarize-protein-nt-panel.py'),str(d),'--nodes',str(r/'reference/results/taxonomy-nodes.dmp'),'--selection-note','Eight outcome-selected header-filter discrepancies, using identical complete mate sequences in both nucleotide tasks. This is a mechanism diagnostic, not a prevalence or accuracy sample.'],check=True)
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/preparation-discrepancies/nt-validation/'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri,'--recursive','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
