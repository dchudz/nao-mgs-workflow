"""Verify and stage a frozen private benchmark; does not start a search."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

name, instance, launched = sys.argv[1:]
root = Path('/home/ec2-user/paired-study-20261024')
source = root/'sources/c77cdde0'
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
prefix = 's3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/'
account = json.loads(subprocess.check_output(aws+['sts','get-caller-identity']))['Account']
assert account == '839993661825'
root.mkdir(exist_ok=True)
def fetch(key, path, expected):
    subprocess.run(aws+['s3','cp',prefix+key,str(path),'--only-show-errors'],check=True)
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual == expected, (str(path), actual)
    return actual
archive=root/'source-c77cdde0.tar.gz'
archive_hash=fetch('source/c77cdde0.tar.gz', archive, '119e9c3697bb5f20ee5227e4447f7eb58b976c39c914200a6ae35ab825a832c9')
source.mkdir(parents=True,exist_ok=True)
with tarfile.open(archive) as tar:
    tar.extractall(source,filter='data')
link=source/'analysis-tools'
if not link.exists(): link.symlink_to('experiments/amino_acids')
(root/'tools').mkdir(exist_ok=True)
sampler=root/'tools/select_pairs'
sampler_hash=fetch('tools/select_pairs-linux',sampler,'e12efbcbbf3cb2ab87bf390b3a1a21ed7efbcad30aae64d202264d3edba29e0d')
sampler.chmod(0o755)
(root/'source-commit.txt').write_text('c77cdde08b03834a9c1f54f3d80b0c9668bf9a83\n')
bound={'stop_utc':'2026-10-25T05:07:00+00:00','instance_id':instance,'started_utc':launched,
       'worker_rate_bound_USD_per_hour':2,'shared_new_compute_allowance_USD':200,
       'all_worker_compute_allocation_bound_USD':192,'shutdown_behavior':'terminate, retain volume',
       'deadline_type':'Renewable spending safeguard, not a user deadline'}
(root/'execution-bound.json').write_text(json.dumps(bound,indent=2)+'\n')
report_image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/pandas:91e85c3fd8e9026c'
subprocess.run(['docker','pull',report_image],check=True)
image_id=subprocess.check_output(['docker','image','inspect',report_image,'--format','{{.Id}}'],text=True).strip()
assert image_id=='sha256:7ab4d47babe91a05ed38d629d7081a3c9a5278c9a1b24127d2d717efa71cfdb4'
sys.path.insert(0,str(source/'experiments/paired_pipeline_comparison'))
from synthetic_control_gate import verify_control
control=verify_control(source,'/home/ec2-user/experiment/protein/tools/diamond','pipeline-e2e-rust:local')
timer=subprocess.check_output(['systemctl','list-timers','paired-study-budget-stop.timer','--no-pager'],text=True)
assert '2026-10-25 05:07:00 UTC' in timer
proof={'worker':name,'instance_id':instance,'account':account,'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'source_sha256':archive_hash,'sampler_sha256':sampler_hash,'synthetic_control':control,
       'report_container_id':image_id,'timer':timer,'execution_bound':bound}
path=root/'bootstrap-verification.json'
path.write_text(json.dumps(proof,indent=2)+'\n')
subprocess.run(aws+['s3','cp',str(path),prefix+'workers/'+name+'/'+path.name,'--only-show-errors'],check=True)
print(json.dumps(proof),flush=True)
