"""Install a cost-bounded guard and record placement; launch no analysis."""
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

payload=json.loads(sys.argv[1]);worker=payload['worker'];plan=payload['plan']
assert worker in {'primary','spot2'} and plan['status']=='accepted-for-execution'
root=Path('/home/ec2-user/paired-study-20261024');source=root/'sources/c77cdde0'
assert (root/'source-commit.txt').read_text().strip()=='c77cdde08b03834a9c1f54f3d80b0c9668bf9a83'
assert hashlib.sha256((source/'experiments/paired_pipeline_comparison/selection.json').read_bytes()).hexdigest()==plan['selection_sha256']
selection=json.loads((source/'experiments/paired_pipeline_comparison/selection.json').read_text())
rows=plan['workers'][worker];runs=rows['runs']
assert len(runs)==20 and len(set(runs))==20 and not set(runs)&set(selection['pilot_runs'])
assert set(runs)<={r['run_accession'] for r in selection['selected']}
for sample in runs:
    input_path=root/'inputs'/sample
    if worker=='spot2':
        scratch=Path('/mnt/paired-study-inputs')
        assert scratch.is_mount()
        target=scratch/sample;target.mkdir(exist_ok=True)
        if input_path.is_symlink():assert input_path.resolve()==target
        else:
            assert not input_path.exists();input_path.symlink_to(target,target_is_directory=True)
    for depth in ([5,25] if sample in selection['depth_runs'] else [5]):
        case=root/'comparisons'/f'{sample}_{depth}m';case.mkdir(exist_ok=True)
        assert not (case/'complete.json').exists(), 'Never reassign completed units'
        placement={'worker':worker,'instance_id':payload['instance_id'],'shard':worker+'-expansion','reason':'Frozen size-balanced allocation after the feasibility gate'}
        (case/'execution-worker.json').write_text(json.dumps(placement,indent=2)+'\n')
shard={'name':worker+'-expansion','runs':runs,'stage':'frozen-main','source_commit':'c77cdde08b03834a9c1f54f3d80b0c9668bf9a83','selection_sha256':plan['selection_sha256']}
(root/(worker+'-expansion-shard.json')).write_text(json.dumps(shard,indent=2)+'\n')
old=(root/'execution-bound.json').read_text();(root/'execution-bound-before-expansion.json').write_text(old)
stop=dt.datetime.fromisoformat(plan['stop_utc'])
assert 24*3600 < (stop-dt.datetime.now(dt.timezone.utc)).total_seconds() <=25*3600
unit='paired-study-budget-stop-expansion'
subprocess.run(['sudo','systemd-run','--unit='+unit,'--on-calendar='+stop.strftime('%Y-%m-%d %H:%M:%S UTC'),'/sbin/shutdown','-h','now'],check=True,capture_output=True,text=True)
for oldunit in ['paired-study-budget-stop.timer','paired-study-budget-stop-recovery.timer','paired-study-budget-stop-24h.timer']:
    subprocess.run(['sudo','systemctl','stop',oldunit],capture_output=True)
bound={'stop_utc':plan['stop_utc'],'instance_id':payload['instance_id'],'deadline_type':'Renewable shared cost guard, not a user deadline','new_compute_allowance_USD':200,'guard_compute_ceiling_USD':plan['guard_compute_ceiling_USD'],'spend_at_allocation_ceiling_USD':plan['elapsed_compute_ceiling_USD'],'remaining_rate_ceiling_USD_per_hour':6,'shutdown_behavior':'stop, retain original volume' if worker=='primary' else 'terminate Spot instance, retain volumes'}
(root/'execution-bound.json').write_text(json.dumps(bound,indent=2)+'\n')
(root/'expansion-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
timer=subprocess.check_output(['systemctl','list-timers',unit+'.timer','--no-pager'],text=True)
assert unit+'.timer' in timer
aws=['aws','--profile','pipeline-codex','--region','us-east-1'];prefix='s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/workers/'+worker+'-expansion/'
for name in [worker+'-expansion-shard.json','execution-bound.json','expansion-plan.json']:
    subprocess.run(aws+['s3','cp',str(root/name),prefix+name,'--only-show-errors'],check=True)
print(json.dumps({'worker':worker,'instance_id':payload['instance_id'],'stop_utc':plan['stop_utc'],'runs':runs,'runtime_config_sha256':hashlib.sha256((root/'runtime.config').read_bytes()).hexdigest(),'timer':timer,'analysis_launched':False}))
