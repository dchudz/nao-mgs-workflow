"""Collect verified reports, preserve each new unit in Git, and refresh summaries."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

HOSTS = {'primary':'paired-study', 'spot1':'paired-spot1', 'spot2':'paired-spot2'}

def read(path):
    return json.loads(path.read_text())

def run(args, cwd=None, **kwargs):
    return subprocess.run(list(map(str,args)), cwd=cwd, check=True, **kwargs)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--worker', choices=HOSTS, required=True)
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--local', type=Path, required=True)
    p.add_argument('--ssh-config', type=Path, required=True)
    p.add_argument('--commit', action='store_true', required=True)
    a = p.parse_args(); repo=a.repo.resolve(); local=a.local.resolve()
    study=repo/'experiments/paired_pipeline_comparison'; results=study/'results'
    # This collector is invoked sequentially; concurrent Git mutations are unsafe.
    pending=subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)
    assert not pending.strip(), 'Commit or isolate existing work before collecting reports'
    host=HOSTS[a.worker]
    response=run(['ssh','-F',a.ssh_config,host,'python3.12','-',a.worker],
                 input=(study/'export_completed_reports.py').read_text(),text=True,capture_output=True)
    metadata=json.loads(response.stdout)
    remote=f'/home/ec2-user/paired-study-20261024/report-export-{a.worker}.tar.gz'
    assert metadata['archive']==remote
    assert re.fullmatch(r'[0-9a-f]{64}',metadata['sha256'])
    assert all(re.fullmatch(r'SRR[0-9]+_(5|25)m',label) for label in metadata['labels'])
    selection=read(study/'selection.json'); orders=read(study/'validation-order.json')
    planned={f"{r['run_accession']}_{n}m" for r in selection['selected']
             for n in ([5,25] if r['run_accession'] in selection['depth_runs'] else [5])}
    assert len(metadata['labels'])==len(set(metadata['labels']))
    assert set(metadata['labels'])<=planned
    local.mkdir(parents=True,exist_ok=True)
    archive=local/f"report-export-{a.worker}-{metadata['sha256'][:12]}.tar.gz"
    run(['scp','-F',a.ssh_config,f'{host}:{remote}',archive])
    assert archive.stat().st_size==metadata['bytes']
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==metadata['sha256']
    (archive.with_suffix('.receipt.json')).write_text(json.dumps(metadata,indent=2)+'\n')
    new=[]; updated=[]
    from aggregate_reports import collect
    with tempfile.TemporaryDirectory(prefix='accepted-reports-',dir=local) as directory:
        extracted=Path(directory)
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                path=Path(member.name)
                assert not path.is_absolute() and '..' not in path.parts and path.parts[0] in metadata['labels']
                assert not member.issym()
                if member.islnk():
                    target=Path(member.linkname)
                    assert not target.is_absolute() and '..' not in target.parts and target.parts[0] in metadata['labels']
            tar.extractall(extracted,filter='data')
        accepted={r['label']:r for r in collect(extracted,selection,orders)}
        for label in metadata['labels']:
            assert accepted[label]['status']=='complete', (label,accepted[label])
            source=extracted/label; dest=results/label
            if not (dest/'complete.json').exists():
                # A launch may record placement before completed reports arrive.
                assert not dest.exists() or {p.name for p in dest.iterdir()}=={'execution-worker.json'}
                shutil.copytree(source,dest,dirs_exist_ok=True)
                run(['git','add','--',dest],cwd=repo)
                run(['git','commit','--quiet','-m',f'Record {label} paired results and verified preservation'],cwd=repo)
                new.append(label)
            else:
                for name in ('complete.json','completed-evidence-backup.json','comparison-report/accounting.json'):
                    assert read(source/name)==read(dest/name), f'Previously preserved result changed: {label}/{name}'
                for name in ('input-acquisition.json','execution-overhead.json'):
                    if not (dest/name).exists() or read(source/name)!=read(dest/name):
                        shutil.copyfile(source/name,dest/name);updated.append(f'{label}/{name}')
    if new or updated:
        commands=[('aggregate_reports.py',['--results',results,'--selection',study/'selection.json','--orders',study/'validation-order.json','--out',results/'aggregate']),
                  ('annotate_runtime.py',['--results',results,'--selection',study/'selection.json','--provenance',study/'runtime-provenance/parallel-feasibility','--out',results/'aggregate']),
                  ('report_depth.py',['--aggregate',results/'aggregate/runs.tsv','--selection',study/'selection.json','--out',results/'aggregate']),
                  ('write_progress.py',['--study',study]),
                  ('report_compute_split.py',['--study',study])]
        for script,args in commands:run(['python3',study/script,*args],cwd=repo)
        run(['git','add','--',results],cwd=repo)
        run(['git','commit','--quiet','-m','Refresh accepted comparison inventory and complete execution receipts'],cwd=repo)
    print(json.dumps({'worker':a.worker,'new_units':new,'updated_receipts':updated,
                      'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()}))

if __name__=='__main__':
    main()
