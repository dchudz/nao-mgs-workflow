"""Observe already-running workers after client exit; never launch scientific work."""
import datetime as dt,json,shlex,subprocess,time
from pathlib import Path
L=Path(__file__).resolve().parent
W=Path('/Users/david/.codex/worktrees/amino-acid-screen/pipeline')
S=W/'experiments/paired_pipeline_comparison'
R='/home/ec2-user/paired-study-20261024'
AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
workers={'primary':('paired-study','i-02d118993a6f4989a',10343),'spot2':('paired-spot2','i-0638e9b2262f8f52c',10881)}
def now():return dt.datetime.now(dt.timezone.utc).isoformat()
def emit(**x):print(json.dumps({'utc':now(),**x}),flush=True)
def read(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def run(args):return subprocess.check_output(args,text=True,timeout=240)
def observe(worker):
    host,instance,pid=workers[worker]
    code="""import json,os
from pathlib import Path
r=Path(%r); p=r/'batch-logs'/%r
def read(n):
 f=p/n
 return json.loads(f.read_text()) if f.exists() else None
print(json.dumps(dict(alive=Path('/proc/%s').exists(),failure=read('failure.json'),finished=read('finished.json'),completed=sorted(x.parent.name for x in (r/'comparisons').glob('*/completed-evidence-backup.json')))))
"""%(R,worker+'-resume',pid)
    return json.loads(run(['ssh','-F',str(L/'ssh_config'),host,'python3.12 -c '+shlex.quote(code)]))
def collect(worker):
    result=run(['python3',str(S/'collect_worker_reports.py'),'--worker',worker,'--repo',str(W),'--local',str(L),'--ssh-config',str(L/'ssh_config'),'--commit'])
    emit(collection=worker,result=result.strip())
def stop(worker,finished):
    host,instance,pid=workers[worker]
    save(L/(worker+'-resume-finished.json'),finished)
    inv=json.loads(run(AWS+['ec2','describe-instances','--instance-ids',instance]))
    assert all(not d['Ebs']['DeleteOnTermination'] for d in inv['Reservations'][0]['Instances'][0]['BlockDeviceMappings'])
    action='stop-instances' if worker=='primary' else 'terminate-instances'
    result=json.loads(run(AWS+['ec2',action,'--instance-ids',instance]))
    save(L/(worker+'-resume-shutdown-request.json'),dict(requested_utc=now(),inventory=inv,result=result))
    emit(worker=worker,shutdown_requested=action,volumes_retained=True)
    target='stopped' if worker=='primary' else 'terminated'
    while True:
        inv=json.loads(run(AWS+['ec2','describe-instances','--instance-ids',instance]))
        if inv['Reservations'][0]['Instances'][0]['State']['Name']==target:break
        time.sleep(15)
    stamp=now();save(L/(worker+'-resume-final-state.json'),dict(observed_utc=stamp,inventory=inv))
    ledger=read(L/'worker-intervals.json');periods=[r for r in ledger['intervals'] if r['instance_id']==instance and r['ended_utc_upper'] is None]
    assert len(periods)==1
    periods[0].update(ended_utc_upper=stamp,evidence=worker+'-resume-final-state.json');save(L/'worker-intervals.json',ledger)
    emit(worker=worker,shutdown_verified=True,compute_interval_closed=stamp)
def spend():
    path=L/('spend-supervised-'+dt.datetime.now(dt.timezone.utc).strftime('%H%M')+'.json')
    emit(spend=run(['python3',str(S/'record_spend.py'),'--intervals',str(L/'worker-intervals.json'),'--out',str(path)]).strip())
    assert read(path)['elapsed_compute_ceiling_USD']<200
done={};last={};next_spend=0
save(L/'client-restart-observation.json',dict(observed_utc=now(),meaning='Reattached to surviving remote output pipes; monitor only. No jobs launched.',workers=workers))
while len(done)<len(workers):
    for worker in workers:
        if worker in done:continue
        x=observe(worker)
        assert not x['failure'],x
        if any(not (S/'results'/label/'completed-evidence-backup.json').exists() for label in x['completed']):collect(worker)
        if x['finished']:
            assert x['finished']['all_assigned_complete'],x
            collect(worker);stop(worker,x['finished']);done[worker]=dict(exit_code=0);continue
        assert x['alive'],'Worker exited without final receipt: '+worker
        if last.get(worker)!=len(x['completed']):emit(worker=worker,preserved_comparisons=len(x['completed']),surviving_process=True);last[worker]=len(x['completed'])
    if time.monotonic()>=next_spend:spend();next_spend=time.monotonic()+1800
    if len(done)<len(workers):time.sleep(45)
spend()
save(L/'resume-supervisor-finished.json',dict(finished_utc=now(),workers=done,monitor='watch_surviving_workers.py'))
emit(all_resume_batches_finished=True,shutdown_verified_for_both=True)
