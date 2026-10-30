"""Wait for first800 completion and verify private evidence before releasing this worker."""
import datetime,hashlib,json,subprocess,tarfile,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');held=r/'cohort/heldout';closure=r/'screen-closure';closure.mkdir(exist_ok=True)
rows=json.loads((r/'heldout-screen-first800.json').read_text());expected={x['run_accession'] for x in rows};assert len(expected)==800
stop=datetime.datetime.fromisoformat('2026-10-23T10:00:00+00:00')
while True:
 observed={p.parent.name for p in held.glob('*/complete.json')};assert observed<=expected
 failures=list(held.glob('*/failure.json'));assert not failures,failures
 result=r/'heldout-screen-remaining.json.results.json'
 if observed==expected and result.exists():break
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Incomplete producer before budget reserve'
 time.sleep(30)
results=json.loads(result.read_text());assert results and all('error' not in x for x in results)
for row in rows:
 evidence=json.loads((held/row['run_accession']/'complete.json').read_text())
 assert evidence['observed_pairs']==5000000
 assert evidence['spot_ranges']==row['spot_ranges']
 assert evidence['run_accession']==row['run_accession']
assert not list(held.glob('*/reads.fastq')) and not list(held.glob('*/*.partial'))
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
aws=['aws','--profile','pipeline-codex','--region','us-east-1']
def verify(root,prefix,name):
 subprocess.run(['python3',str(r/'verify-private-prefix.py'),str(root),'runs/amino-acid-20261021/'+prefix,'--output',str(closure/name)],check=True,stdout=subprocess.DEVNULL)
verify(held,'cohort-screens/','heldout-backup-verification.json')
# Existing development records and interrupted attempts are kept as operational evidence.
for directory,prefix in [('cohort/development','screen-development/'),('screen-resize-interruption','screen-resize-interruption/')]:
 subprocess.run(aws+['s3','cp',str(r/directory)+'/',uri+prefix,'--recursive','--only-show-errors'],check=True)
 verify(r/directory,prefix,directory.replace('/','-')+'-verification.json')
files=[p for p in r.iterdir() if p.is_file()]
for directory in ['experiments','screen-burden','cohort/interrupted-thread-run']:
 files.extend(p for p in (r/directory).rglob('*') if p.is_file() and not p.is_symlink())
manifest={str(p.relative_to(r)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
(closure/'source-log-manifest.json').write_text(json.dumps(manifest,indent=2))
archive=closure/'screen-source-logs.tar.gz'
with tarfile.open(archive,'w:gz') as tar:
 for p in files:tar.add(p,arcname=str(p.relative_to(r)))
record={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'instance':'i-0dd39a46c4efb6c4c','heldout_samples':800,'heldout_input_pairs':4000000000,'frozen_manifest_sha256':hashlib.sha256((r/'heldout-screen-first800.json').read_bytes()).hexdigest(),'source_log_archive_bytes':archive.stat().st_size,'source_log_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'raw_audits_relocated_to':'aa-experiment','safe_to_release_screen_producer':True}
(closure/'complete.json').write_text(json.dumps(record,indent=2))
subprocess.run(aws+['s3','cp',str(closure)+'/',uri+'screen-closure/','--recursive','--exclude','complete.json','--only-show-errors'],check=True)
subprocess.run(aws+['s3','cp',str(closure/'complete.json'),uri+'screen-closure/complete.json','--only-show-errors'],check=True)
verify(closure,'screen-closure/','closure-backup-verification.json')
print(json.dumps(record),flush=True)
