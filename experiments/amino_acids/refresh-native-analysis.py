"""Analyze one immutable set of committed batches while later batches continue."""
import datetime,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');storage=r/'cohort-storage'
markers=sorted((storage/'batches').glob('*/complete.json'))
samples=[a for p in markers for a in json.loads(p.read_text())['samples']]
assert len(samples)==len(set(samples))
snapshot=storage/('analysis-inputs-'+str(len(samples)));snapshot.mkdir(exist_ok=True)
for p in markers:
 target=snapshot/p.parent.name
 if not target.exists():target.symlink_to(p.parent,target_is_directory=True)
assert len(list(snapshot.glob('*/complete.json')))==len(markers)
subprocess.run(['aws','--profile','pipeline-codex','s3','sync','s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/cohort-screens/',str(storage/'screen-metrics'),'--exclude','*','--include','*/complete.json','--include','*/nt.log','--include','*/nt20.log','--include','*/nt22.log','--only-show-errors'],check=True)
for cmd in [['analyze-cohort.py',str(snapshot),str(storage/'analysis')],['analyze-resources.py',str(snapshot)],['analyze-species.py',str(snapshot)],['summarize-species-support.py',str(storage/'analysis')],['analyze-clustered.py'],['analyze-screen-burden.py',str(storage/'screen-metrics'),str(storage/'screen-burden')]]:
 with (storage/'analysis'/(cmd[0]+'.log')).open('w') as log:subprocess.run(['python3',str(r/cmd[0]),*cmd[1:]],stdout=log,stderr=log,check=True)
d=storage/'analysis';summary=json.loads((d/'summary.json').read_text());assert summary['completed_samples']==len(samples)
for name in ['resource-summary.json','species-summary.json']:assert json.loads((d/name).read_text())['samples']==len(samples),name
assert all(x['samples']==len(samples) for x in json.loads((d/'clustered-summary.json').read_text())['arms'].values())
(d/'analysis-snapshot.json').write_text(json.dumps({'samples':samples,'batches':[p.parent.name for p in markers],'analyzed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print(json.dumps({'analyzed_samples':len(samples),'arms':{k:v['totals'] for k,v in summary['arms'].items()}}),flush=True)
