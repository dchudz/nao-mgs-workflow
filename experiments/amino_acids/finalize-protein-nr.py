"""Complete NR-side descriptive analyses and verified backups before releasing the frozen standalone search."""
import datetime,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/combined-validation';out=r/'protein/final-analysis';out.mkdir(exist_ok=True)
base='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/';aws=['aws','--profile','pipeline-codex','--region','us-east-1']
stop=datetime.datetime.fromisoformat('2026-10-23T09:00:00+00:00')
assert subprocess.check_output(aws+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
def run(script,*args):
 command=['python3',str(r/script),*map(str,args)]
 with (out/'commands.jsonl').open('a') as f:f.write(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':command})+'\n')
 print('START '+script,flush=True)
 with (out/(Path(script).name+'.log')).open('w') as f:subprocess.run(command,cwd=r,stdout=f,stderr=f,check=True)
 print('COMPLETE '+script,flush=True)
def backup(directory,prefix):
 subprocess.run(aws+['s3','cp',str(directory)+'/',base+prefix+'/', '--recursive','--only-show-errors'],check=True)
 run('verify-private-prefix.py',directory,'runs/amino-acid-20261021/'+prefix,'--output',directory/'backup-verification.json')
while True:
 if (d/'complete.json').exists() and (d/'search-raw-complete.json').exists():
  result=subprocess.run(aws+['s3api','head-object','--bucket','pipeline-lab-839993661825-us-east-1','--key','runs/amino-acid-20261021/protein-finalization/complete.json'],capture_output=True,text=True)
  if result.returncode==0:break
  if '(404)' not in result.stderr and 'Not Found' not in result.stderr:raise RuntimeError(result.stderr)
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'No complete primary analysis before the financial safeguard'
 time.sleep(60)
for source,dest in [('protein-finalization/complete.json','primary-complete.json'),('protein200-native-comparison-nt20/pair-evidence.tsv','nt20-pair-evidence.tsv')]:subprocess.run(aws+['s3','cp',base+source,str(out/dest),'--only-show-errors'],check=True)
assert json.loads((out/'primary-complete.json').read_text())['private_backups_verified']
run('summarize-diamond-stages.py',d/'nr-search.log',d/'stage-timings.json')
run('annotate-nr-competitors.py')
run('describe-nr-resolution.py',out/'nt20-pair-evidence.tsv',r/'protein/current-nr',r/'protein/resolution')
run('audit-protein-cases.py','--evidence',out/'nt20-pair-evidence.tsv','--protein',r/'cohort-storage/protein200','--nr',d/'nr.tsv','--taxonomy',r/'protein/current-nr','--output',r/'protein/cases')
backup(r/'protein/resolution','protein-resolution')
backup(r/'protein/cases','protein-cases')
backup(d,'protein-nr-combined')
# Back up analysis provenance before the gate; its receipt is copied separately.
subprocess.run(aws+['s3','cp',str(out)+'/',base+'protein-nr-final-analysis/','--recursive','--only-show-errors'],check=True)
marker={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_joins_complete':True,'nr_descriptive_analyses_complete':True,'private_backups_verified':True,'standalone_scope_unchanged':True}
staged=r/'nr-primary-analysis-ready.json';staged.write_text(json.dumps(marker,indent=2))
subprocess.run(aws+['s3','cp',str(staged),base+'protein-nr-final-analysis/complete.json','--only-show-errors'],check=True)
staged.replace(r/'nr-primary-analysis-complete.json')
print('NR_PRIMARY_ANALYSIS_COMPLETE_STANDALONE_RELEASED',flush=True)
