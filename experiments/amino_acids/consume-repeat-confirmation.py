"""Prospectively fixed repeat-filter comparison on the untouched last500manifest samples."""
import csv,datetime,gzip,json,shutil,subprocess,sys,tarfile,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import union
protocol=json.loads((r/'repeat-confirmation-protocol.json').read_text());order=protocol['samples'];selected=set(order);d=r/'cohort-storage/repeat-confirmation';d.mkdir(exist_ok=True);(d/'batches').mkdir(exist_ok=True)
bucket='pipeline-lab-839993661825-us-east-1';prefix='runs/amino-acid-20261021/cohort-screens/'
config=r/'experiments/amino_acids/repeat-confirmation.config';config.write_text("includeConfig 'heldout.config'\nexecutor.cpus=4\nexecutor.memory='90 GB'\nprocess { withLabel: small { cpus=2 }; withLabel: max { cpus=4 } }\n")
deadline=datetime.datetime.fromisoformat('2026-10-23T12:00:00+00:00')
def done():return {s for p in (d/'batches').glob('*/complete.json') for s in json.loads(p.read_text())['samples']}
def available():
 keys=json.loads(subprocess.check_output(['aws','--profile','pipeline-codex','s3api','list-objects-v2','--bucket',bucket,'--prefix',prefix,'--query','Contents[?ends_with(Key, `/complete.json`)].Key','--output','json'],text=True)) or []
 ready={k.split('/')[-2] for k in keys};return [s for s in order if s in ready-done()]
def copy(src,dst):subprocess.run(['aws','--profile','pipeline-codex','s3','cp',src,str(dst),'--only-show-errors'],check=True)
while len(done())<len(order):
 if datetime.datetime.now(datetime.timezone.utc)>deadline:print('SCHEDULING_DEADLINE',flush=True);break
 ready=available()
 if len(ready)<25:time.sleep(30);continue
 samples=ready[:25];batch=d/'batches'/f'b{1+len(list((d/"batches").glob("*/complete.json"))):03d}';batch.mkdir();incoming=batch/'incoming';incoming.mkdir();screen=[];start=time.monotonic()
 with (batch/'manifest.tsv').open('w') as manifest:
  manifest.write('id\tpath\n')
  for sample in samples:
   s=incoming/sample;s.mkdir()
   for name in ['nt.fastq','aa12.fastq','complete.json']:copy(f's3://{bucket}/{prefix}{sample}/{name}',s/name)
   original=json.loads((s/'complete.json').read_text());assert original['observed_pairs']==5000000
   cmd=[str(r/'rust-tools/target-speed/release/screen_amino_acids'),'screen','--index',str(r/'reference/virus.aa12-repeat-filtered.bin'),'--reads',str(s/'aa12.fastq'),'--matched',str(s/'aa12-repeat.fastq'),'--threads','2']
   with (s/'screen.log').open('w') as log:subprocess.run(['/usr/bin/time','-f','{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"peak_rss_kib":%M}','-o',str(s/'resources.json')]+cmd,stdout=log,stderr=log,check=True,timeout=1800)
   stats=union(s/'nt.fastq',s/'aa12-repeat.fastq',s/'union.fastq',s/'rescued.fastq');assert stats['amino_pairs']<=original['amino_pairs'];assert stats['nucleotide_pairs']==original['nucleotide_pairs']
   screen.append(dict(sample=sample,original_aa12_pairs=original['amino_pairs'],original_union_pairs=original['union_pairs'],**stats,subset_screen_resources=json.loads((s/'resources.json').read_text())))
   manifest.write(f'{sample}_repeatunion\t{s/"union.fastq"}\n')
   for name in ['aa12.fastq','nt.fastq','aa12-repeat.fastq','rescued.fastq']:(s/name).unlink()
 (batch/'screen-metrics.json').write_text(json.dumps(screen,indent=2))
 cmd=[str(r/'nextflow'),'-log',str(batch/'nextflow.log'),'run','experiments/amino_acids/heldout-workflow.nf','-c',str(config),'-work-dir',str(batch/'work'),'--manifest',str(batch/'manifest.tsv'),'--ref_dir',str(r/'reference'),'--adapters',str(r/'ref/adapters.fasta'),'--output',str(batch/'output')]
 with (batch/'run.log').open('w') as log:subprocess.run(cmd,cwd=r,stdout=log,stderr=log,check=True,timeout=5400)
 for folder,pattern in [('hits','*_virus_hits.tsv.gz'),('marked','similarity-duplicate-marked*.tsv.gz'),('counts','*_clade_counts.tsv.gz')]:
  assert len(list((batch/'output'/folder).glob(pattern)))==len(samples),(folder,batch)
 with tarfile.open(batch/'task-diagnostics.tar.gz','w:gz') as tar:
  for p in (batch/'work').glob('*/*/.command*'):
   if p.is_file() and not p.is_symlink():tar.add(p,arcname=str(p.relative_to(batch/'work')))
 uri=f's3://{bucket}/runs/amino-acid-20261021/repeat-confirmation/{batch.name}/'
 subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(batch)+'/',uri,'--recursive','--exclude','work/*','--exclude','incoming/*','--only-show-errors'],check=True)
 marker={'samples':samples,'wall_seconds':time.monotonic()-start,'private_results':uri,'primary_protocol_unchanged':True};(batch/'complete.json').write_text(json.dumps(marker,indent=2));copy(str(batch/'complete.json'),uri+'complete.json')
 shutil.rmtree(batch/'work');shutil.rmtree(incoming);print(json.dumps({'completed_samples':len(done()),'planned':len(order),'batch':batch.name,'wall_seconds':marker['wall_seconds']}),flush=True)
(d/'scope.json').write_text(json.dumps({'completed':sorted(done()),'planned':order,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));copy(str(d/'scope.json'),f's3://{bucket}/runs/amino-acid-20261021/repeat-confirmation/scope.json')
print('REPEAT_CONFIRMATION_STOPPED',flush=True)
