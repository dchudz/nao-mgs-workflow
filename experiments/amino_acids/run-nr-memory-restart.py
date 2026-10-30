"""Verify exact hardware/thread parity before restarting the unchanged full-reference search."""
import collections,datetime,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');prior=r/'protein/nr-thread-parity-on-demand';d=r/'protein/nr-thread-parity-memory-restart';d.mkdir(exist_ok=True)
assert os.cpu_count()==32
assert not (d/'parity.json').exists(),'Do not repeat a completed parity trial'
old=json.loads((prior/'parity.json').read_text());cmd=old['command'].copy();assert cmd[cmd.index('--threads')+1]=='24'
shutil.copy2(prior/'original.fna',d/'queries.fna');cmd[cmd.index('--threads')+1]='32';cmd[cmd.index('--query')+1]=str(d/'queries.fna');cmd[cmd.index('--out')+1]=str(d/'threads32.tsv')
start=time.monotonic()
with (d/'run.log').open('w') as log:subprocess.run(['/usr/bin/time','-v']+cmd,stdout=log,stderr=log,check=True,timeout=3600)
def rows(path):
 with path.open() as f:return collections.Counter(line.rstrip('\n') for line in f)
expected=rows(prior/'original.tsv');actual=rows(d/'threads32.tsv')
result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'wall_seconds':time.monotonic()-start,'expected_rows':sum(expected.values()),'actual_rows':sum(actual.values()),'removed_rows':sum((expected-actual).values()),'added_rows':sum((actual-expected).values()),'exact_alignment_multiset_parity':expected==actual,'query_sha256':hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),'hardware':'r7a.8xlarge,32physicalAMDcores,256GiB memory; comparisonoriginal r7i.4xlarge4threads and earlier m8azn.6xlarge24threads','scope':'Same54468development/full-anchorquerymates; complete nr.000 volume; frozen sensitivity/block/chunk/scoring/output parameters.'}
(d/'parity.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/','s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/nr-thread-parity-memory-restart/','--recursive','--only-show-errors'],check=True)
assert expected==actual,'Do not adopt resource change without exact fixture alignment parity'
del expected,actual
subprocess.run(['python3','-u',str(r/'run-resized-nr.py')],check=True)
