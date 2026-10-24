"""Foreground resume supervisor: preserve each completed unit and stop idle workers."""
import argparse
import datetime as dt
import json
import queue
import shlex
import subprocess
import threading
import time
from pathlib import Path

AWS = ['aws','--profile','pipeline-codex','--region','us-east-1']
REMOTE = '/home/ec2-user/paired-study-20261024'
SOURCE = REMOTE+'/sources/c77cdde0'
COMMIT = 'c77cdde08b03834a9c1f54f3d80b0c9668bf9a83'


def now(): return dt.datetime.now(dt.timezone.utc).isoformat()
def read(p): return json.loads(p.read_text())
def save(p,x): p.write_text(json.dumps(x,indent=2)+'\n')
def emit(x): print(json.dumps({'utc':now(),**x}),flush=True)


def main():
    parser=argparse.ArgumentParser()
    for name in ('repo','local','ssh-config','plan'):parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args();study=a.repo/'experiments/paired_pipeline_comparison'
    plan=read(a.plan);assert plan['status']=='accepted-for-execution'
    hosts={worker:(record['ssh_alias'],record['instance_id']) for worker,record in plan['workers'].items()}
    assert plan['guard_compute_ceiling_USD']<200
    assert subprocess.check_output(['git','status','--porcelain'],cwd=a.repo,text=True).strip()==''
    # Re-check the existing acceptance gate immediately before launching work.
    import sys
    sys.path.insert(0,str(study));from aggregate_reports import collect
    selection=read(study/'selection.json');rows=collect(study/'results',selection,read(study/'validation-order.json'))
    pilot=[r for r in rows if r['feasibility_run']]
    assert len(pilot)==12 and all(r['status']=='complete' for r in pilot)
    planned=[sample for w in plan['workers'].values() for sample in w['runs']]
    assert len(planned)==len(set(planned))
    units={sample+('_5m' if level=='shallow' else '_25m') for record in plan['workers'].values() for sample,levels in record['levels'].items() for level in levels}
    assert len(units)==plan['remaining_comparisons']==43
    assert units=={r['label'] for r in rows if r['status']!='complete'}
    q=queue.Queue();jobs={};ended={}
    def ssh(worker,command):return ['ssh','-F',str(a.ssh_config),hosts[worker][0],command]
    def pump(worker,process,log):
        with log.open('w') as f:
            for line in process.stdout:
                f.write(line);f.flush();q.put((worker,line.rstrip()))
        q.put((worker,{'exit_code':process.wait()}))
    def collect_worker(worker):
        result=subprocess.run(['python3',str(study/'collect_worker_reports.py'),'--worker',worker,'--repo',str(a.repo),'--local',str(a.local),'--ssh-config',str(a.ssh_config),'--commit'],text=True,capture_output=True,timeout=180)
        if result.returncode:raise RuntimeError('Report collection failed: '+result.stderr[-2000:])
        if result.stdout.strip():emit({'collection':worker,'result':result.stdout.strip()})
    def remote_completed(worker):
        code="from pathlib import Path; import json; print(json.dumps(sorted(p.parent.name for p in Path('"+REMOTE+"/comparisons').glob('*/completed-evidence-backup.json'))))"
        return json.loads(subprocess.check_output(ssh(worker,'python3.12 -c '+shlex.quote(code)),text=True,timeout=45))
    def stop_worker(worker):
        instance=hosts[worker][1]
        inventory=json.loads(subprocess.check_output(AWS+['ec2','describe-instances','--instance-ids',instance],text=True))
        inst=inventory['Reservations'][0]['Instances'][0]
        assert all(not item['Ebs']['DeleteOnTermination'] for item in inst['BlockDeviceMappings'])
        action='stop-instances' if worker=='primary' else 'terminate-instances'
        result=json.loads(subprocess.check_output(AWS+['ec2',action,'--instance-ids',instance],text=True))
        save(a.local/(worker+'-resume-shutdown-request.json'),dict(requested_utc=now(),inventory=inventory,result=result))
        emit({'worker':worker,'shutdown_requested':action,'volumes_retained':True})
        condition='instance-stopped' if worker=='primary' else 'instance-terminated'
        subprocess.run(AWS+['ec2','wait',condition,'--instance-ids',instance],check=True,timeout=600)
        final=json.loads(subprocess.check_output(AWS+['ec2','describe-instances','--instance-ids',instance],text=True))
        stamp=now();save(a.local/(worker+'-resume-final-state.json'),dict(observed_utc=stamp,inventory=final))
        ledger_path=a.local/'worker-intervals.json';ledger=read(ledger_path)
        periods=[x for x in ledger['intervals'] if x['instance_id']==instance and x['ended_utc_upper'] is None]
        assert len(periods)==1
        periods[0].update(ended_utc_upper=stamp,evidence=worker+'-resume-final-state.json')
        save(ledger_path,ledger)
        emit({'worker':worker,'shutdown_verified':True,'compute_interval_closed':stamp})
    for worker in hosts:
        command=['python3.12',REMOTE+'/resume-tools/run_remaining.py','--root',REMOTE,'--source',SOURCE,'--retained','/home/ec2-user/experiment','--shard',REMOTE+'/'+worker+'-resume-shard.json','--source-commit',COMMIT]
        process=subprocess.Popen(ssh(worker,shlex.join(command)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        jobs[worker]=process
        threading.Thread(target=pump,args=(worker,process,a.local/(worker+'-resume-foreground.log')),daemon=False).start()
        emit({'launched':worker,'instance':hosts[worker][1],'runs':len(plan['workers'][worker]['runs'])})
    next_collect=time.monotonic()+90;next_spend=time.monotonic()+1800
    while len(ended)<len(jobs):
        try:worker,event=q.get(timeout=5)
        except queue.Empty:event=None
        if event is not None:
            emit({'worker':worker,'event':event})
            if isinstance(event,dict):
                ended[worker]=event
                if event['exit_code']!=0:raise RuntimeError(worker+' failed; preserve its current state and inspect before retry')
                collect_worker(worker)
                finished=json.loads(subprocess.check_output(ssh(worker,'cat '+REMOTE+'/batch-logs/'+worker+'-resume/finished.json'),text=True,timeout=45))
                save(a.local/(worker+'-resume-finished.json'),finished)
                if not finished['all_assigned_complete']:raise RuntimeError(worker+' reached the guard before finishing; inspect its preserved progress')
                stop_worker(worker)
        if time.monotonic()>=next_collect:
            for worker in jobs:
                if worker in ended:continue
                labels=remote_completed(worker)
                if any(not (study/'results'/label/'completed-evidence-backup.json').exists() for label in labels):collect_worker(worker)
            next_collect=time.monotonic()+120
        if time.monotonic()>=next_spend:
            dest=a.local/('spend-supervised-'+dt.datetime.now(dt.timezone.utc).strftime('%H%M')+'.json')
            result=subprocess.run(['python3',str(study/'record_spend.py'),'--intervals',str(a.local/'worker-intervals.json'),'--out',str(dest)],capture_output=True,text=True,timeout=120)
            if result.returncode:raise RuntimeError('Cost recording failed: '+result.stderr[-2000:])
            emit({'spend':result.stdout.strip()});assert read(dest)['elapsed_compute_ceiling_USD']<plan['guard_compute_ceiling_USD']
            next_spend=time.monotonic()+1800
    save(a.local/'resume-supervisor-finished.json',dict(finished_utc=now(),workers=ended))
    emit({'all_resume_batches_finished':True,'shutdown_requested_for_both':True})


if __name__=='__main__':main()
