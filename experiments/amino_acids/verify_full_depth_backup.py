"""Read back every saved sample artifact before releasing reconstructable scratch."""
import argparse,datetime,hashlib,json,os,shutil,subprocess,tarfile,tempfile
from pathlib import Path
R=Path('/home/ec2-user/experiment');PREFIX='s3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/'
AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
p=argparse.ArgumentParser();p.add_argument('sample');p.add_argument('--release-scratch',action='store_true');p.add_argument('--release-published-hsps',action='store_true');a=p.parse_args()
assert not a.release_published_hsps or a.release_scratch,'Release HSP copies only with verified scratch cleanup'
allowed={r['run_accession'] for r in json.loads((R/'experiments/amino_acids/end-to-end-full-depth/sample-selection.json').read_text())['initial_samples']}
assert a.sample in allowed
D=R/'e2e'/a.sample
previous=D/'backup-verification.json'
if previous.exists():
    assert not any(json.loads(previous.read_text()).get(key) for key in ['locally_evicted_hsps','local_hsp_eviction_plan']),'Restore the recorded S3 HSPs before repeating full verification'
complete=json.loads((D/'complete.json').read_text());assert complete['all_four_outputs_completed']
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
# Prune scratch before traversal, and stream one independent S3 readback at a
# time. Verification must not require a second complete copy of large HSP files.
files=[]
for directory,dirs,names in os.walk(D):
    dirs[:]=[name for name in dirs if name not in {'work','.nextflow'}]
    files.extend(p for name in names if (p:=Path(directory)/name).is_file() and included(p))
size=sum(p.stat().st_size for p in files)
manifest=[];local_hashes={}
for path in sorted(files):
    rel=path.relative_to(D)
    stat=path.stat();key=(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns)
    if key not in local_hashes:local_hashes[key]=digest(path)
    sha=local_hashes[key];remote_hash=hashlib.sha256();remote_bytes=0
    with tempfile.TemporaryFile() as errors:
        proc=subprocess.Popen(AWS+['s3','cp',PREFIX+'samples/'+a.sample+'/'+str(rel),'-','--only-show-errors'],
                              stdout=subprocess.PIPE,stderr=errors)
        try:
            for block in iter(lambda:proc.stdout.read(16*1024**2),b''):
                remote_hash.update(block);remote_bytes+=len(block)
            proc.stdout.close();code=proc.wait(timeout=60)
        finally:
            if proc.poll() is None:proc.kill();proc.wait()
        if code:
            errors.seek(0);raise RuntimeError(f'S3 readback failed for {rel}: '+errors.read(4096).decode(errors='replace'))
    assert remote_bytes==stat.st_size and remote_hash.hexdigest()==sha,('Independent readback differs',rel)
    assert path.stat().st_mtime_ns==stat.st_mtime_ns,('Artifact changed during verification',rel)
    manifest.append({'path':str(rel),'bytes':stat.st_size,'sha256':sha})
(D/'file-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'verified':True,'sample':a.sample,'files':len(manifest),'bytes':size,'source_commit':complete['source_commit'],
         'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'method':'Every retained artifact independently streamed from private S3 and compared by byte count and SHA256, using bounded memory and no duplicate disk copy.',
         'manifest_sha256':digest(D/'file-manifest.json'),'released_scratch':False}
(D/'backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name in ['file-manifest.json','backup-verification.json']:
    subprocess.run(AWS+['s3','cp',str(D/name),PREFIX+'samples/'+a.sample+'/'+name,'--only-show-errors'],check=True)
if a.release_scratch:
    # Published results and per-task diagnostic archives have all passed readback.
    # Shared reference caches, source code, data receipts and scientific outputs stay.
    # A retained symlink into scratch would become invalid after cleanup.
    removals=list(D.glob('*/work'))+[R/'data'/a.sample/name for name in [a.sample,'R1.fastq.gz','R2.fastq.gz','R1.partial','R2.partial']]
    for path in files:
        if path.is_symlink():assert not any(path.resolve().is_relative_to(p) for p in removals),path
    for work in D.glob('*/work'):
        assert (work.parent/'task-diagnostics.tar.gz').is_file()
        shutil.rmtree(work)
    raw=R/'data'/a.sample
    for name in [a.sample,'R1.fastq.gz','R2.fastq.gz','R1.partial','R2.partial']:
        path=raw/name
        if path.exists():path.unlink()
    receipt['released_scratch']=True
    if a.release_published_hsps:
        # These exact immutable outputs have just passed independent S3 readback.
        # Keep their original manifest entries and object keys after local eviction.
        names={f'protein-run/published/{stage}/protein-hits.tsv' for stage in ['PROTEIN_ALIGN','PROTEIN_CANDIDATES']}
        selected=[record for record in manifest if record['path'] in names]
        assert {record['path'] for record in selected}==names
        assert len({(record['bytes'],record['sha256']) for record in selected})==1
        for record in selected:
            path=D/record['path']
            assert path.is_file() and not path.is_symlink() and path.stat().st_size==record['bytes']
        receipt['local_hsp_eviction_plan']=[{**record,'s3_uri':PREFIX+'samples/'+a.sample+'/'+record['path']} for record in selected]
        # Persist the recovery plan before unlinking, including across interruption.
        (D/'backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
        subprocess.run(AWS+['s3','cp',str(D/'backup-verification.json'),PREFIX+'samples/'+a.sample+'/backup-verification.json','--only-show-errors'],check=True)
        receipt['locally_evicted_hsps']=[]
        for record in receipt['local_hsp_eviction_plan']:
            (D/record['path']).unlink()
            receipt['locally_evicted_hsps'].append(record)
    (D/'backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    subprocess.run(AWS+['s3','cp',str(D/'backup-verification.json'),PREFIX+'samples/'+a.sample+'/backup-verification.json','--only-show-errors'],check=True)
print(json.dumps(receipt))
