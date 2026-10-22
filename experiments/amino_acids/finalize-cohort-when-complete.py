"""Finalize existing cohort analyses only after the frozen complete sample sets exist."""
import datetime,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');c=r/'cohort-storage';out=c/'analysis';out.mkdir(exist_ok=True)
expected={x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())};assert len(expected)==1000
repeat=set(json.loads((r/'repeat-confirmation-protocol.json').read_text())['samples']);assert len(repeat)==500 and repeat<=expected
stop=datetime.datetime.fromisoformat('2026-10-23T12:00:00+00:00')
def completed(root):
 values=[s for p in root.glob('*/complete.json') for s in json.loads(p.read_text())['samples']]
 assert len(values)==len(set(values)),'Duplicate committed sample'
 return set(values)
def run(name,*args):
 with (out/(name+'.log')).open('w') as log:subprocess.run(['python3',str(r/name),*map(str,args)],stdout=log,stderr=log,check=True)
while True:
 seen=completed(c/'batches');assert seen<=expected
 if seen==expected:break
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Cohort incomplete before budget reserve'
 time.sleep(30)
for sample in expected:assert (c/'screen-metrics'/sample/'complete.json').exists(),sample
run('analyze-cohort.py',c/'batches',out)
run('analyze-resources.py')
run('analyze-species.py')
run('analyze-clustered.py')
run('summarize-species-support.py',out)
assert json.loads((out/'summary.json').read_text())['completed_samples']==1000
assert json.loads((out/'resource-summary.json').read_text())['samples']==1000
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(out)+'/',uri+'heldout-analysis/','--recursive','--only-show-errors'],check=True)
print('FULL1000_NATIVE_ANALYSES_COMPLETE',flush=True)
while True:
 seen=completed(c/'repeat-confirmation/batches');assert seen<=repeat
 if seen==repeat:break
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Repeat cohort incomplete before budget reserve'
 time.sleep(30)
run('analyze-repeat-confirmation.py')
rd=c/'repeat-confirmation/analysis';assert json.loads((rd/'summary.json').read_text())['completed_paired_samples']==500
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(rd)+'/',uri+'repeat-confirmation/analysis/','--recursive','--only-show-errors'],check=True)
marker={'native_samples':1000,'repeat_samples':500,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(out/'final-cohort-analyses.json').write_text(json.dumps(marker,indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(out/'final-cohort-analyses.json'),uri+'heldout-analysis/final-cohort-analyses.json','--only-show-errors'],check=True)
print(json.dumps(marker),flush=True)
