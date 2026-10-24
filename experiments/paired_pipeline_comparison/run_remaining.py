"""Resume explicitly unfinished units with the frozen scientific driver.

No selection depends on biological results. Each run uses the frozen sampler and
pipeline driver. Complete input subsets are backed up and read back before the
re-downloadable full-source cache is removed.
"""
import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import tarfile
import time
from pathlib import Path

p=argparse.ArgumentParser()
for key in ('root','source','retained','shard'):p.add_argument('--'+key,type=Path,required=True)
p.add_argument('--source-commit',required=True);a=p.parse_args()
root=a.root;code=a.source/'experiments/paired_pipeline_comparison'
aws=['aws','--profile','pipeline-codex','--region','us-east-1']
prefix='s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/'
assert json.loads(subprocess.check_output(aws+['sts','get-caller-identity'],text=True))['Account']=='839993661825'
selection=json.loads((code/'selection.json').read_text())
shard=json.loads(a.shard.read_text());assert len(set(shard['runs']))==len(shard['runs'])
assert set(shard['levels'])==set(shard['runs'])
for sample, levels in shard['levels'].items():
    expected=['shallow','deep'] if sample in selection['depth_runs'] else ['shallow']
    assert levels and len(levels)==len(set(levels)) and set(levels)<=set(expected)
    for level in levels:
        label=sample+('_5m' if level=='shallow' else '_25m')
        assert not (root/'comparisons'/label/'completed-evidence-backup.json').exists(), 'Never rerun a preserved completed unit'
assert set(shard['runs'])<={r['run_accession'] for r in selection['selected']}
logdir=root/'batch-logs'/shard['name'];logdir.mkdir(parents=True,exist_ok=True)

def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc)
def remaining():
    bound=json.loads((root/'execution-bound.json').read_text())
    return (datetime.datetime.fromisoformat(bound['stop_utc'])-now()).total_seconds()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def upload(p,key):
    subprocess.run(aws+['s3','cp',str(p),prefix+key,'--only-show-errors'],check=True,timeout=1800)
def run(label,args):
    seconds=remaining()-1800
    if seconds<300:raise RuntimeError('Preservation reserve reached')
    record={'command':list(map(str,args)),'started_utc':now().isoformat(),'source_commit':a.source_commit,'execution_driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    save(logdir/(label+'.started.json'),record);start=time.monotonic()
    try:
        with (logdir/(label+'.log')).open('w') as f:
            result=subprocess.run(record['command'],stdout=f,stderr=subprocess.STDOUT,timeout=seconds)
        record['exit_code']=result.returncode
    except Exception as e:
        record.update(exit_code=-1,exception=str(e));raise
    finally:
        record.update(seconds=time.monotonic()-start,finished_utc=now().isoformat())
        save(logdir/(label+'.json'),record)
        upload(logdir/(label+'.json'),'workers/'+shard['name']+'/'+label+'.json')
    if record['exit_code']:raise RuntimeError('Failed: '+label)
    print(json.dumps({'finished':label,'seconds':record['seconds']}),flush=True)
def readback(p,key):
    upload(p,key);h=hashlib.sha256();size=0
    proc=subprocess.Popen(aws+['s3','cp',prefix+key,'-','--only-show-errors'],stdout=subprocess.PIPE)
    for block in iter(lambda:proc.stdout.read(8*1024**2),b''):h.update(block);size+=len(block)
    assert proc.wait()==0 and size==p.stat().st_size and h.hexdigest()==sha(p)
    return {'s3':prefix+key,'bytes':size,'sha256':h.hexdigest(),'readback_verified':True}

done=[]
try:
    for sample in shard['runs']:
        deep='deep' in shard['levels'][sample]
        # Reserve time for a full paired unit, not just downloading its input.
        if remaining()<(5*3600 if deep else 3*3600):
            save(logdir/'bounded-stop.json',{'next_run':sample,'reason':'Insufficient remaining guarded time for another complete unit','completed_runs':done})
            break
        inputs=root/'inputs'/sample;inputs.mkdir(parents=True,exist_ok=True)
        if not (inputs/'verified.json').exists():
            archive=inputs/sample;partial=inputs/(sample+'.partial')
            if not archive.exists():
                run(sample+'-download',aws+['s3','cp',f's3://sra-pub-run-odp/sra/{sample}/{sample}',partial,'--no-sign-request','--only-show-errors'])
                partial.rename(archive)
            run(sample+'-sample',['python3.12',code/'acquire_subset.py','--root',root,'--sample',sample,
                '--selection',code/'selection.json','--tools',a.retained/'sratoolkit.3.4.1-ubuntu64/bin','--selector',root/'tools/select_pairs'])
        acq=json.loads((inputs/'verified.json').read_text());objects=[]
        for item in acq['files']:
            f=Path(item['path']);assert sha(f)==item['sha256'] and f.stat().st_size==item['bytes']
            objects.append(readback(f,f'inputs/{sample}/{item["level"]}/{f.name}'))
        provenance=inputs/'sampling-provenance.tar.gz'
        if not provenance.exists():
            with tarfile.open(provenance,'w:gz') as tar:
                for path in sorted(inputs.rglob('*')):
                    if path.is_file() and (path.name in ('ids.tsv','complete.json','verified.json','sra-stat.xml','vdb-validate.log') or path.suffix=='.log'):
                        tar.add(path,arcname=str(path.relative_to(inputs)))
        objects.append(readback(provenance,f'inputs/{sample}/{provenance.name}'))
        upload(inputs/'verified.json',f'inputs/{sample}/verified.json')
        save(inputs/'backup-readback.json',{'objects':objects,'verified_utc':now().isoformat()})
        upload(inputs/'backup-readback.json',f'inputs/{sample}/backup-readback.json')
        # Only remove this study's reproducible source/decoding cache. Selected
        # input subsets and every pipeline result remain both local and private-S3.
        removed=[]
        for path in (inputs/sample,inputs/'chunks'):
            if path.exists():
                if path.is_dir():shutil.rmtree(path)
                else:path.unlink()
                removed.append(str(path))
        save(inputs/'source-cache-evicted.json',{'paths':removed,'source_sha256':acq['source_sha256'],'selected_inputs_retained':True,'backup_readback_verified':True})
        upload(inputs/'source-cache-evicted.json',f'inputs/{sample}/source-cache-evicted.json')
        for level in shard['levels'][sample]:
            label=sample+('_5m' if level=='shallow' else '_25m');d=root/'comparisons'/label
            if not (d/'complete.json').exists():
                run(label+'-pair',['python3.12',code/'run_pair.py','--root',root,'--source',a.source,
                    '--retained',a.retained,'--sample',sample,'--level',level,'--source-commit',a.source_commit])
            run(label+'-integrity',['python3.12',code/'check_bowtie_outputs.py',d])
            if not (d/'comparison-report/accounting.json').exists():
                run(label+'-report',['python3.12',code/'report_pair.py','--root',root,'--source',a.source,'--retained',a.retained,'--sample',label])
            summary=d/'completed-evidence.tar.gz'
            with tarfile.open(summary,'w:gz') as tar:
                for path in sorted(d.rglob('*')):
                    if not path.is_file() or path.is_symlink() or path==summary:continue
                    relative=path.relative_to(d)
                    if 'work' not in relative.parts and '.nextflow' not in relative.parts:
                        if any(part.startswith('.') for part in relative.parts):continue
                        tar.add(path,arcname=str(relative))
                    elif path.name in ('.command.sh','.command.log','.command.err','.command.out','.command.trace','.exitcode'):
                        tar.add(path,arcname=str(relative))
            proof=readback(summary,f'comparisons/{label}/{summary.name}')
            save(d/'completed-evidence-backup.json',proof)
            upload(d/'completed-evidence-backup.json',f'comparisons/{label}/completed-evidence-backup.json')
        done.append(sample)
        save(logdir/'progress.json',{'completed_runs':done,'remaining_runs':[s for s in shard['runs'] if s not in done]})
        upload(logdir/'progress.json','workers/'+shard['name']+'/progress.json')
    save(logdir/'finished.json',{'completed_runs':done,'assigned_runs':shard['runs'],'all_assigned_complete':done==shard['runs'],'finished_utc':now().isoformat()})
except Exception as e:
    save(logdir/'failure.json',{'error':str(e),'completed_runs':done,'failed_utc':now().isoformat()})
    raise
finally:
    subprocess.run(aws+['s3','sync',str(logdir),prefix+'workers/'+shard['name']+'/','--only-show-errors'],check=True,timeout=1800)
