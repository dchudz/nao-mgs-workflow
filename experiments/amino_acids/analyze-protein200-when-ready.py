"""Collect completed candidate metadata and compare it with native complete pairs."""
import datetime,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/protein200';d.mkdir(exist_ok=True)
base='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
files=['all.tsv.gz','all.members.tsv.gz','viral-pairs.tsv','viral-pairs.tsv.json']
for name in files:
 while True:
  result=subprocess.run(['aws','--profile','pipeline-codex','s3','cp',base+'protein200/'+name,str(d/name),'--only-show-errors'],stderr=subprocess.PIPE,text=True)
  if result.returncode==0:break
  assert '404' in result.stderr or 'does not exist' in result.stderr,result.stderr
  assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime.fromisoformat('2026-10-22T02:00:00+00:00')
  time.sleep(30)
out=d/'native-comparison';out.mkdir(exist_ok=True)
subprocess.run(['python3',str(r/'analyze-protein200.py'),'--protein',str(d),'--output',str(out)],check=True)
subprocess.run(['python3',str(r/'analyze-protein-overlap.py'),str(out/'pair-evidence.tsv'),str(out/'overlap-summary.json')],check=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(out)+'/',base+'protein200-native-comparison/','--recursive','--only-show-errors'],check=True)
