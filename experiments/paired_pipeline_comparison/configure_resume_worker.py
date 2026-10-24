"""Restore execution metadata and cost guards; change no scientific configuration."""
import datetime as dt,hashlib,json,subprocess,sys
from pathlib import Path
payload=json.loads(sys.argv[1]);worker=payload['worker'];plan=payload['plan'];record=plan['workers'][worker]
root=Path('/home/ec2-user/paired-study-20261024');source=root/'sources/c77cdde0';retained=Path('/home/ec2-user/experiment')
assert (root/'source-commit.txt').read_text().strip()==plan['source_commit']
assert hashlib.sha256((source/'experiments/paired_pipeline_comparison/selection.json').read_bytes()).hexdigest()==plan['selection_sha256']
assert hashlib.sha256((root/'runtime.config').read_bytes()).hexdigest()=='d061955a7dc733fb7fd1eea58ea0b1e30fb7bcdb726ac73059e0b82e1161f29a'
sys.path.insert(0,str(source/'experiments/paired_pipeline_comparison'))
from synthetic_control_gate import verify_control
control=verify_control(source,retained/'protein/tools/diamond','pipeline-e2e-rust:local')
assert control['passed']
selection=json.loads((source/'experiments/paired_pipeline_comparison/selection.json').read_text())
assert len(record['runs'])==len(set(record['runs']))==19
assert set(record['levels'])==set(record['runs'])
archived=[]
for sample in record['runs']:
    input_path=root/'inputs'/sample
    if worker=='spot2':
        scratch=Path('/mnt/paired-study-inputs');assert scratch.is_mount()
        target=scratch/sample;assert target.is_dir() and input_path.is_symlink() and input_path.resolve()==target
    for level in record['levels'][sample]:
        assert level=='shallow' or level=='deep' and sample in selection['depth_runs']
        label=sample+('_5m' if level=='shallow' else '_25m');case=root/'comparisons'/label
        assert not (case/'completed-evidence-backup.json').exists() and not (case/'complete.json').exists()
        old=json.loads((case/'execution-worker.json').read_text()) if (case/'execution-worker.json').exists() else None
        if case.exists() and {p.name for p in case.iterdir()}-{'execution-worker.json'}:
            dest=root/'interrupted-attempts'/('halted-'+label);dest.parent.mkdir(exist_ok=True);assert not dest.exists();case.rename(dest);archived.append(str(dest))
        case.mkdir(exist_ok=True)
        (case/'execution-worker.json').write_text(json.dumps({'worker':worker,'instance_id':record['instance_id'],'shard':record['name'],'reason':'Resume unchanged frozen comparison after an erroneous halt','prior_placement':old,'retained_reference_volume':record['boot_volume_id'] if worker=='primary' else 'vol-05ca1445d265f5325'},indent=2)+'\n')
shard={'name':record['name'],'runs':record['runs'],'levels':record['levels'],'source_commit':plan['source_commit'],'selection_sha256':plan['selection_sha256']}
(root/(worker+'-resume-shard.json')).write_text(json.dumps(shard,indent=2)+'\n')
stop=dt.datetime.fromisoformat(plan['stop_utc']);assert 3600<(stop-dt.datetime.now(dt.timezone.utc)).total_seconds()<24*3600
unit='paired-study-budget-stop-resume'
subprocess.run(['sudo','systemd-run','--unit='+unit,'--on-calendar='+stop.strftime('%Y-%m-%d %H:%M:%S UTC'),'/sbin/shutdown','-h','now'],check=True,capture_output=True,text=True)
for oldunit in ['paired-study-budget-stop','paired-study-budget-stop-recovery','paired-study-budget-stop-24h','paired-study-budget-stop-expansion']:
    subprocess.run(['sudo','systemctl','stop',oldunit+'.timer'],capture_output=True)
old=root/'execution-bound-before-resume.json';assert not old.exists();old.write_bytes((root/'execution-bound.json').read_bytes())
(root/'execution-bound.json').write_text(json.dumps({'stop_utc':plan['stop_utc'],'instance_id':record['instance_id'],'deadline_type':'Renewable shared cost guard, not a user deadline','new_compute_allowance_USD':200,'guard_compute_ceiling_USD':plan['guard_compute_ceiling_USD'],'shutdown_behavior':'stop' if worker=='primary' else 'terminate; retain disks'},indent=2)+'\n')
(root/'resume-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
timer=subprocess.check_output(['systemctl','list-timers',unit+'.timer','--no-pager'],text=True);assert unit+'.timer' in timer
proof={'worker':worker,'instance_id':record['instance_id'],'verified_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'control':control,'source_commit':plan['source_commit'],'runtime_config_sha256':hashlib.sha256((root/'runtime.config').read_bytes()).hexdigest(),'driver_sha256':hashlib.sha256((root/'resume-tools/run_remaining.py').read_bytes()).hexdigest(),'archived_interrupted_attempts':archived,'units':sum(len(x) for x in record['levels'].values()),'timer':timer,'analysis_launched':False}
(root/'resume-configuration.json').write_text(json.dumps(proof,indent=2)+'\n')
aws=['aws','--profile','pipeline-codex','--region','us-east-1'];assert json.loads(subprocess.check_output(aws+['sts','get-caller-identity'],text=True))['Account']=='839993661825'
prefix='s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/workers/'+record['name']+'/'
for name in [worker+'-resume-shard.json','execution-bound.json','resume-plan.json','resume-configuration.json']:
    subprocess.run(aws+['s3','cp',str(root/name),prefix+name,'--only-show-errors'],check=True)
print(json.dumps(proof))
