"""Adopt extra threads only after exact single-volume alignment parity."""
import collections,datetime,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/nr-thread-parity';d.mkdir(parents=True,exist_ok=True)
deadline=datetime.datetime.fromisoformat('2026-10-22T02:00:00+00:00')
while not (r/'protein/current-nr/download-complete.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<deadline,'Reference not ready before acceleration deadline'
 time.sleep(30)
for name in ['original.fna','original.tsv','summary.json']:
 subprocess.run(['aws','--profile','pipeline-codex','s3','cp','s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/nr-dedup-benchmark/'+name,str(d/name),'--only-show-errors'],check=True)
old=json.loads((d/'summary.json').read_text());cmd=old['runs']['original']['command']
assert cmd[cmd.index('--threads')+1]=='4';cmd[cmd.index('--threads')+1]='96';cmd[cmd.index('--query')+1]=str(d/'original.fna');cmd[cmd.index('--out')+1]=str(d/'threads96.tsv')
start=time.monotonic()
with (d/'threads96.log').open('w') as f:subprocess.run(['/usr/bin/time','-v']+cmd,stdout=f,stderr=f,check=True,timeout=3600)
def rows(path):
 with open(path) as f:return collections.Counter(line.rstrip('\n') for line in f)
expected=rows(d/'original.tsv');actual=rows(d/'threads96.tsv');removed=expected-actual;added=actual-expected
result={'command':cmd,'wall_seconds':time.monotonic()-start,'expected_rows':sum(expected.values()),'actual_rows':sum(actual.values()),'removed_rows':sum(removed.values()),'added_rows':sum(added.values()),'exact_alignment_multiset_parity':expected==actual,'hardware':'c6a.24xlarge,96vCPU AMD; original r7i.4xlarge,4threads Intel','query_scope':'54468 development and full-depth candidate reads; single nr.000 volume, same frozen search parameters.'}
(d/'parity.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/nr-thread-parity/','--recursive','--only-show-errors'],check=True)
assert expected==actual,'Thread/hardware change failed raw parity; do not adopt accelerated full search'
del expected,actual,removed,added
subprocess.run(['python3',str(r/'run-development-nr-fast.py')],check=True)
