"""Preserve finished raw search evidence before potentially lengthy aggregation."""
import datetime,json,re,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/combined-validation'
deadline=datetime.datetime.fromisoformat('2026-10-23T09:00:00+00:00')
while True:
 p=d/'nr-search.log';tail=''
 if p.exists():
  with p.open('rb') as f:f.seek(max(0,p.stat().st_size-8192));tail=f.read().decode(errors='replace')
 if re.search(r'^\d+ queries aligned\.',tail,re.M) and 'Exit status: 0' in tail and 'Command terminated by signal' not in tail and 'Command exited with non-zero status' not in tail:break
 if datetime.datetime.now(datetime.timezone.utc)>deadline:raise RuntimeError('No completed NR search before raw-backup deadline')
 time.sleep(30)
base='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/protein-nr-combined/'
files={name:(d/name).stat().st_size for name in ['nr.tsv','nr-search.log','search-protocol.json','combined-queries.fna']}
for name in files:subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/name),base+name,'--only-show-errors'],check=True)
record={'search_completed_and_raw_backed_up_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_bytes':files,'interpretation':'Raw full-reference search finished successfully. This marker does not assert completed taxonomy aggregation or final scientific analysis.'}
(d/'search-raw-complete.json').write_text(json.dumps(record,indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'search-raw-complete.json'),base+'search-raw-complete.json','--only-show-errors'],check=True)
print(json.dumps(record),flush=True)
