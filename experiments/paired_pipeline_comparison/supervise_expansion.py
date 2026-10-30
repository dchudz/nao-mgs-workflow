"""Foreground supervisor: preserve each completed unit and stop idle workers."""
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
HOSTS = {'primary':('paired-study','i-02d118993a6f4989a'), 'spot2':('paired-spot2','i-0ee9da3f09d197be6')}


def now(): return dt.datetime.now(dt.timezone.utc).isoformat()
def read(p): return json.loads(p.read_text())
def save(p,x): p.write_text(json.dumps(x,indent=2)+'\n')
def emit(x): print(json.dumps({'utc':now(),**x}),flush=True)


def main():
    parser=argparse.ArgumentParser()
    for name in ('repo','local','ssh-config','plan'):parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args();study=a.repo/'experiments/paired_pipeline_comparison'
    plan=read(a.plan);assert plan['status']=='accepted-for-execution'
    assert plan['guard_compute_ceiling_USD']<200
    assert subprocess.check_output(['git','status','--porcelain'],cwd=a.repo,text=True).strip()==''
    # Re-check the existing acceptance gate immediately before launching work.
    import sys
    sys.path.insert(0,str(study));from aggregate_reports import collect
    selection=read(study/'selection.json');rows=collect(study/'results',selection,read(study/'validation-order.json'))
    pilot=[r for r in rows if r['feasibility_run']]
    assert len(pilot)==12 and all(r['status']=='complete' for r in pilot)
    planned=[sample for w in plan['workers'].values() for sample in w['runs']]
    assert len(planned)==len(set(planned))==40
    assert set(planned)=={r['run_accession'] for r in selection['selected']}-set(selection['pilot_runs'])
    q=queue.Queue();jobs={};ended={}
    def ssh(worker,command):return ['ssh','-F',str(a.ssh_config),HOSTS[worker][0],command]
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
        instance=HOSTS[worker][1]
        inventory=json.loads(subprocess.check_output(AWS+['ec2','describe-instances','--instance-ids',instance],text=True))
        inst=inventory['Reservations'][0]['Instances'][0]
        assert all(not item['Ebs']['DeleteOnTermination'] for item in inst['BlockDeviceMappings'])
        action='stop-instances' if worker=='primary' else 'terminate-instances'
        result=json.loads(subprocess.check_output(AWS+['ec2',action,'--instance-ids',instance],text=True))
        save(a.local/(worker+'-expansion-shutdown-request.json'),dict(requested_utc=now(),inventory=inventory,result=result))
        emit({'worker':worker,'shutdown_requested':action,'volumes_retained':True})
    for worker in HOSTS:
        command=['python3.12',SOURCE+'/experiments/paired_pipeline_comparison/run_batch.py','--root',REMOTE,'--source',SOURCE,'--retained','/home/ec2-user/experiment','--shard',REMOTE+'/'+worker+'-expansion-shard.json','--source-commit',COMMIT]
        process=subprocess.Popen(ssh(worker,shlex.join(command)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        jobs[worker]=process
        threading.Thread(target=pump,args=(worker,process,a.local/(worker+'-expansion-foreground.log')),daemon=False).start()
        emit({'launched':worker,'instance':HOSTS[worker][1],'runs':len(plan['workers'][worker]['runs'])})
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
                finished=json.loads(subprocess.check_output(ssh(worker,'cat '+REMOTE+'/batch-logs/'+worker+'-expansion/finished.json'),text=True,timeout=45))
                save(a.local/(worker+'-expansion-finished.json'),finished)
                if not finished['all_assigned_complete']:raise RuntimeError(worker+' reached the guard before finishing; inspect its preserved progress')
                stop_worker(worker)
        if time.monotonic()>=next_collect:
            for worker in jobs:
                if worker in ended:continue
                labels=remote_completed(worker)
                if any(not (study/'results'/label/'completed-evidence-backup.json').exists() for label in labels):collect_worker(worker)
            next_collect=time.monotonic()+120
        if time.monotonic()>=next_spend and not ended:
            dest=a.local/('spend-supervised-'+dt.datetime.now(dt.timezone.utc).strftime('%H%M')+'.json')
            result=subprocess.run(['python3',str(study/'record_spend.py'),'--intervals',str(a.local/'worker-intervals.json'),'--out',str(dest)],capture_output=True,text=True,timeout=120)
            if result.returncode:raise RuntimeError('Cost recording failed: '+result.stderr[-2000:])
            emit({'spend':result.stdout.strip()});assert read(dest)['elapsed_compute_ceiling_USD']<plan['guard_compute_ceiling_USD']
            next_spend=time.monotonic()+1800
    save(a.local/'expansion-supervisor-finished.json',dict(finished_utc=now(),workers=ended))
    emit({'all_expansion_batches_finished':True,'shutdown_requested_for_both':True})


if __name__=='__main__':main()
