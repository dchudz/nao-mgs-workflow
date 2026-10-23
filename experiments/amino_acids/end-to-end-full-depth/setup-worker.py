"""Stage frozen references and pinned tools on the dedicated private worker."""
import concurrent.futures,datetime,gzip,hashlib,json,shutil,subprocess,tarfile,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'setup';D.mkdir(exist_ok=True)
AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
PREFIX='s3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/'

def run(name,cmd):
    previous=[json.loads(s) for s in (D/'events.jsonl').read_text().splitlines()] if (D/'events.jsonl').exists() else []
    if any(r['stage']==name for r in previous):
        print(json.dumps({'stage':name,'reusing_completed_command':True}),flush=True)
        return
    t=time.monotonic()
    with (D/(name+'.log')).open('w') as log:
        subprocess.run(list(map(str,cmd)),cwd=R,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=7200)
    receipt={'stage':name,'seconds':time.monotonic()-t,'command':list(map(str,cmd))}
    with (D/'events.jsonl').open('a') as f:f.write(json.dumps(receipt)+'\n')
    print(json.dumps(receipt),flush=True)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()

def references():
    if (D/'reference-verification.json').exists() and (R/'reference/genome-lengths.json').exists() and (R/'reference/virus-genomes-human-masked.fasta.gz').exists():
        assert json.loads((D/'reference-verification.json').read_text())['verified']
        return
    run('references',AWS+['s3','sync','s3://pipeline-lab-839993661825-us-west-2/runs/full-reference-20261020a/index/output/',R/'reference','--only-show-errors'])
    frozen=json.loads((D/'core-reference-archive-complete.json').read_text())
    for item in frozen['files']:
        path=R/'reference/results/blast_db'/item['name']
        assert path.stat().st_size==item['bytes']
        assert sha(path)==item['source_sha256'],path
    # Also pin the viral reference, metadata, taxonomy and native indexes.
    checked=[]
    for item in json.loads((D/'reference-results.json').read_text())['files']:
        path=R/item['path']
        if 'bytes' in item:
            assert path.stat().st_size==item['bytes'] and sha(path)==item['sha256'],path
            checked.append(item['path'])
    (D/'reference-verification.json').write_text(json.dumps({'core_files':len(frozen['files']),'other_files':checked,'verified':True},indent=2))
    lengths={};name=None;length=0
    with gzip.open(R/'reference/results/virus-genomes-masked.fasta.gz','rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:lengths[name]=length
                name=line[1:].split()[0];length=0
            else:length+=len(line.strip())
        if name is not None:lengths[name]=length
    (R/'reference/genome-lengths.json').write_text(json.dumps(lengths))
    run('human-masked',AWS+['s3','cp','s3://pipeline-lab-839993661825-us-west-2/runs/full-reference-20261020a/index/work/ee/68c35db6c192c580af569f6d9237d8/virus-genomes-human-masked.fasta.gz',R/'reference/virus-genomes-human-masked.fasta.gz','--only-show-errors'])

def protein():
    info=json.loads((D/'protein.json').read_text());archive=D/'protein.tar.gz'
    wanted={'protein/tools/diamond','protein/viral-alignment.dmnd','protein/viral-alignment.faa'}
    checks=[f for f in info['files'] if f['path'] in wanted]
    if all((R/f['path']).is_file() and sha(R/f['path'])==f['sha256'] for f in checks):return
    run('protein-restore',AWS+['s3','cp',info['archive']['uri'],archive,'--only-show-errors'])
    assert sha(archive)==info['archive']['sha256']
    with tarfile.open(archive) as tar:
        for name in wanted:tar.extract(tar.getmember(name),R,filter='data')
    (R/'protein/tools/diamond').chmod(0o755)
    archive.unlink()

def build():
    image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/rust-tools@sha256:e69b93d9a1a0ae41c59ec4c0a6a0841f3a7f4871fc2f9b333eaa094563612270'
    run('rust-build',['docker','run','--rm','-v',str(R/'rust-tools')+':/build','-w','/build','rust:1.88-bookworm','cargo','build','--release','--locked','-p','screen_amino_acids','-p','mark_duplicates','-p','mark_duplicates_similarity'])
    context=D/'rust-image';context.mkdir(exist_ok=True)
    for binary in ['mark_duplicates','mark_duplicates_similarity']:
        shutil.copy2(R/'rust-tools/target/release'/binary,context/binary)
    dockerfile=context/'Dockerfile'
    dockerfile.write_text(f'FROM {image}\nCOPY mark_duplicates mark_duplicates_similarity /usr/local/bin/\n')
    run('rust-image',['docker','build','-f',dockerfile,'-t','pipeline-e2e-rust:local',context])
    run('nextflow',['curl','-fL','--retry','3','https://github.com/nextflow-io/nextflow/releases/download/v26.04.6/nextflow-26.04.6-dist','-o',R/'nextflow'])
    (R/'nextflow').chmod(0o755)

identity=json.loads(subprocess.check_output(AWS+['sts','get-caller-identity'],text=True));assert identity['Account']=='839993661825'
with concurrent.futures.ThreadPoolExecutor(3) as pool:
    tasks=[pool.submit(f) for f in [references,protein,build]]
    for task in tasks:task.result()
run('aa-index',['docker','run','--rm','-v',str(R)+':'+str(R),'rust:1.88-bookworm',R/'rust-tools/target/release/screen_amino_acids','index','--reference',R/'reference/virus-genomes-human-masked.fasta.gz','--output',R/'reference/virus.aa12.bin','--k','12'])
(D/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified':True}))
run('backup-setup',AWS+['s3','cp',str(D)+'/',PREFIX+'setup/','--recursive','--exclude','*.tar.gz','--only-show-errors'])
