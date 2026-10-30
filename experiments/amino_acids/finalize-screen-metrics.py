"""Collect all frozen screen records and summarize their1000-sample candidate burden."""
import datetime,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');c=r/'cohort-storage';d=c/'screen-burden-inputs';out=c/'screen-burden';d.mkdir(exist_ok=True)
marker=json.loads((c/'producer-finished.json').read_text());assert marker['samples']==1000 and marker['private_screen_commits_verified']
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/';cmd=['aws','--profile','pipeline-codex','s3','cp',uri+'cohort-screens/',str(d)+'/', '--recursive','--exclude','*']
for name in ['complete.json','nt.log','nt22.log','nt20.log']:cmd+=['--include','*/'+name]
subprocess.run(cmd+['--only-show-errors'],check=True)
expected={x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())};observed={p.parent.name for p in d.glob('*/complete.json')};assert observed==expected and len(observed)==1000
for p in d.glob('*/complete.json'):assert json.loads(p.read_text())['observed_pairs']==5000000
subprocess.run(['python3',str(r/'analyze-screen-burden.py'),str(d),str(out)],check=True)
assert json.loads((out/'summary.json').read_text())['samples']==1000
(out/'collection-provenance.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':uri+'cohort-screens/','samples':1000,'raw_input_pairs':5000000000,'command':cmd,'manifest_sha256':marker['manifest_sha256']},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(out)+'/',uri+'screen-burden/','--recursive','--only-show-errors'],check=True)
print('FULL1000_SCREEN_BURDEN_COMPLETE',flush=True)
