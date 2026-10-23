"""Read back every saved sample artifact before releasing reconstructable scratch."""
import argparse,datetime,hashlib,json,shutil,subprocess,tarfile
from pathlib import Path
R=Path('/home/ec2-user/experiment');PREFIX='s3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/'
AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
p=argparse.ArgumentParser();p.add_argument('sample');p.add_argument('--release-scratch',action='store_true');a=p.parse_args()
allowed={r['run_accession'] for r in json.loads((R/'experiments/amino_acids/end-to-end-full-depth/sample-selection.json').read_text())['initial_samples']}
assert a.sample in allowed
D=R/'e2e'/a.sample;complete=json.loads((D/'complete.json').read_text());assert complete['all_four_outputs_completed']
identity=json.loads(subprocess.check_output(AWS+['sts','get-caller-identity'],text=True));assert identity['Account']=='839993661825'
def included(path):
    parts=path.relative_to(D).parts
    return 'work' not in parts and '.nextflow' not in parts and not path.name.endswith('.partial') and path.name not in {'file-manifest.json','backup-verification.json'}
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
# Preserve generated count/parameter sidecars that are not workflow outputs.
for work in D.glob('*/work'):
    archive=work.parent/'generated-task-metrics.tar.gz'
    metrics=list(work.glob('*/*/*.json'))+list(work.glob('*/*/*.stats.txt'))
    with tarfile.open(archive,'w:gz') as tar:
        for item in metrics:
            if item.is_file() and not item.is_symlink():tar.add(item,arcname=str(item.relative_to(work)))
    subprocess.run(AWS+['s3','cp',str(archive),PREFIX+'samples/'+a.sample+'/'+str(archive.relative_to(D)),'--only-show-errors'],check=True)
files=[p for p in D.rglob('*') if p.is_file() and included(p)]
size=sum(p.stat().st_size for p in files);download=R/'verify-download'/a.sample
assert shutil.disk_usage(R).free>size+50*1024**3, 'Insufficient room for independent backup readback'
download.mkdir(parents=True,exist_ok=True)
subprocess.run(AWS+['s3','sync',PREFIX+'samples/'+a.sample+'/',str(download)+'/', '--only-show-errors'],check=True)
manifest=[];local_hashes={}
for path in sorted(files):
    rel=path.relative_to(D);copy=download/rel
    assert copy.is_file() and copy.stat().st_size==path.stat().st_size,rel
    stat=path.stat();key=(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns)
    if key not in local_hashes:local_hashes[key]=digest(path)
    sha=local_hashes[key];assert digest(copy)==sha,rel
    manifest.append({'path':str(rel),'bytes':stat.st_size,'sha256':sha})
(D/'file-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'verified':True,'sample':a.sample,'files':len(manifest),'bytes':size,'source_commit':complete['source_commit'],
         'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'method':'Every retained artifact independently downloaded from private S3 and compared by SHA256.',
         'manifest_sha256':digest(D/'file-manifest.json'),'released_scratch':False}
(D/'backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name in ['file-manifest.json','backup-verification.json']:
    subprocess.run(AWS+['s3','cp',str(D/name),PREFIX+'samples/'+a.sample+'/'+name,'--only-show-errors'],check=True)
shutil.rmtree(download)
if a.release_scratch:
    # Published results and per-task diagnostic archives have all passed readback.
    # Shared reference caches, source code, data receipts and scientific outputs stay.
    for work in D.glob('*/work'):
        assert (work.parent/'task-diagnostics.tar.gz').is_file()
        shutil.rmtree(work)
    raw=R/'data'/a.sample
    for name in [a.sample,'R1.fastq.gz','R2.fastq.gz','R1.partial','R2.partial']:
        path=raw/name
        if path.exists():path.unlink()
    receipt['released_scratch']=True
    (D/'backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    subprocess.run(AWS+['s3','cp',str(D/'backup-verification.json'),PREFIX+'samples/'+a.sample+'/backup-verification.json','--only-show-errors'],check=True)
print(json.dumps(receipt))
