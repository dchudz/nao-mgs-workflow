"""Release the final partial native batch only after every planned screen is privately committed."""
import datetime,hashlib,importlib.util,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');spec=importlib.util.spec_from_file_location('consumer',r/'consume-cohort-fast.py');consumer=importlib.util.module_from_spec(spec);spec.loader.exec_module(consumer)
expected=set(consumer.order);assert len(expected)==1000
deadline=datetime.datetime.fromisoformat('2026-10-23T04:00:00+00:00')
while consumer.ready()!=expected:
 if datetime.datetime.now(datetime.timezone.utc)>deadline:raise RuntimeError('Not every planned screen was committed before the handoff deadline')
 time.sleep(60)
marker=r/'cohort-storage/producer-finished.json';data={'samples':1000,'private_screen_commits_verified':True,'manifest_sha256':hashlib.sha256((r/'heldout-ranges.json').read_bytes()).hexdigest(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
temporary=marker.with_suffix('.partial');temporary.write_text(json.dumps(data,indent=2));temporary.replace(marker)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(marker),'s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/producer-finished.json','--only-show-errors'],check=True)
print(json.dumps(data),flush=True)
