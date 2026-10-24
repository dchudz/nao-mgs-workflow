"""Run the two frozen pipelines on one verified input subset.

Reuse production modules and previously verified experimental modules unchanged.
Never treat a failed, missing, or unfinished output as an empty biological result.
"""
import argparse,csv,datetime,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
from synthetic_control_gate import verify_control
p=argparse.ArgumentParser()
for flag in ('root','source','retained'):p.add_argument('--'+flag,type=Path,required=True)
p.add_argument('--sample',required=True);p.add_argument('--level',choices=['shallow','deep'],default='shallow')
p.add_argument('--source-commit',required=True);a=p.parse_args()
R=a.source;root=a.root;retained=a.retained;ref=retained/'reference'
control_proof=verify_control(R,retained/'protein/tools/diamond','pipeline-e2e-rust:local')
sample=a.sample+('_5m' if a.level=='shallow' else '_25m');pairs=5_000_000 if a.level=='shallow' else 25_000_000
D=root/'comparisons'/sample;D.mkdir(parents=True,exist_ok=True)
(D/'synthetic-control-gate.json').write_text(json.dumps(control_proof,indent=2)+'\n')
data=root/'inputs'/a.sample/a.level;AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
prefix='s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/'
def now():return datetime.datetime.now(datetime.timezone.utc)
def save(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
deadline=datetime.datetime.fromisoformat(json.loads((root/'execution-bound.json').read_text())['stop_utc'])
def backup():
    # Hardlink only completed regular published artifacts; Nextflow FIFOs never enter staging.
    staging=root/'checkpoint'/sample
    if staging.exists():
        import shutil;shutil.rmtree(staging)
    staging.mkdir(parents=True)
    for directory,dirs,files in os.walk(D):
        dirs[:]=[s for s in dirs if s not in {'work','.nextflow'}]
        for name in files:
            source=Path(directory)/name
            if not source.is_file() or source.is_symlink() or name.endswith('.partial'):continue
            target=staging/source.relative_to(D);target.parent.mkdir(parents=True,exist_ok=True)
            os.link(source,target)
    subprocess.run(AWS+['s3','sync',str(staging)+'/',prefix+'comparisons/'+sample+'/','--only-show-errors'],check=True,timeout=1800)
def run(stage,cmd):
    receipt=D/(stage+'.json');command=list(map(str,cmd))
    if receipt.exists():
        old=json.loads(receipt.read_text());assert old['command']==command and old['source_commit']==a.source_commit and old['exit_code']==0
        return
    remaining=(deadline-now()).total_seconds()-1800
    assert remaining>1800,'Reached shutdown/preservation reserve'
    started=now();t=time.monotonic()
    record={'stage':stage,'command':command,'source_commit':a.source_commit,'started_utc':started.isoformat()}
    save(D/(stage+'.started.json'),record)
    failure=None
    try:
        with (D/(stage+'.log')).open('w') as f:
            result=subprocess.run(command,cwd=R,stdout=f,stderr=subprocess.STDOUT,timeout=remaining)
        code=result.returncode
    except Exception as error:
        code=-1;failure=error;record['exception']=str(error)
    record.update(finished_utc=now().isoformat(),wall_seconds=time.monotonic()-t,exit_code=code)
    save(receipt if code==0 else D/(stage+'.failure.json'),record)
    backup();print(json.dumps(record),flush=True)
    if failure:raise failure
    if code:raise RuntimeError('Stage failed: '+stage)
def one(pattern):
    values=list(D.glob(pattern));assert len(values)==1,(pattern,values);return values[0]
cfg=root/'runtime.config';publish=D/'publish.config'
publish.write_text('process.publishDir = [path: { "${params.base_dir}/published/${task.process.tokenize(\':\')[-1]}" }, mode: \'copy\']\n')
def nf(stage,script,configuration,params,published=False):
    target=D/stage;target.mkdir(exist_ok=True)
    values={**params,'base_dir':str(target),'ref_dir':str(ref)};save(target/'params.json',values)
    cmd=[retained/'nextflow','-log',target/'nextflow.log','run',R,'-main-script',script,'-lib',R/'lib']
    if configuration:cmd+=['-c',R/configuration]
    cmd+=['-c',cfg]
    if published:cmd+=['-c',publish]
    cmd+=['-params-file',target/'params.json','-work-dir',target/'work','-ansi-log','false']
    run(stage,cmd)
assert json.loads(subprocess.check_output(AWS+['sts','get-caller-identity'],text=True))['Account']=='839993661825'
assert (root/'source-commit.txt').read_text().strip()==a.source_commit
acq=json.loads((root/'inputs'/a.sample/'verified.json').read_text())
for item in [r for r in acq['files'] if r['level']==a.level]:
    path=Path(item['path']);assert item['pairs']==pairs and path.stat().st_size==item['bytes'] and sha(path)==item['sha256']
assert len([r for r in acq['files'] if r['level']==a.level])==2
assert json.loads((retained/'integration/parallel-adapter/verification.json').read_text())['verified']
proof=json.loads((retained/'integration/parallel-adapter/verification.json').read_text())
assert proof['implementation_sha256']==sha(R/'experiments/amino_acids/parallel_diamond_adapter.py')
assert proof['frozen_adapter_sha256']==sha(R/'experiments/amino_acids/diamond_adapter.py')
save(D/'protocol.json',{'source_commit':a.source_commit,'sample':sample,'source_run':a.sample,'level':a.level,'pairs':pairs,'acquisition_sha256':sha(root/'inputs'/a.sample/'verified.json'),'arms':['production','protein_variant'],'validation_policy':'Original production settings, equal in both arms','group':'One run and depth per group','source_reference':str(ref),'runtime_config_sha256':sha(cfg),'started_utc':now().isoformat()})
sheet=D/'samplesheet.csv'
with sheet.open('w') as f:
    w=csv.writer(f);w.writerow(['sample','fastq_1','fastq_2']);w.writerow([sample,data/'R1.fastq.gz',data/'R2.fastq.gz'])
nf('baseline-run','main.nf','configs/run.config',{'platform':'illumina','sample_sheet':str(sheet),'adapters':str(R/'ref/adapters.fasta'),'n_reads_profile':1000000,'random_seed':'17310','bt2_score_threshold':20,'bracken_threshold':1,'taxid_artificial':'81077'})
baseline=D/'baseline-run/output/results'
with (baseline/f'{sample}_read_counts.tsv').open() as f:measured=next(csv.DictReader(f,delimiter='\t'))
assert int(float(measured['n_read_pairs']))==pairs
run('aa-screen',[retained/'rust-tools/target/release/screen_amino_acids','screen','--index',ref/'virus.aa12.bin','--reads',data/'R1.fastq.gz','--reads2',data/'R2.fastq.gz','--threads','32','--matched',D/'aa.fastq'])
native=one('baseline-run/output/intermediates/reads/raw_viral/*_nucleaze_match.fastq.gz')
run('union',['python3.12',R/'experiments/amino_acids/union_pairs.py',native,D/'aa.fastq',D/'union.fastq',D/'rescued.fastq',D/'union-counts.json'])
nf('prepare','experiments/paired_pipeline_comparison/prepare.nf',None,{'sample':sample,'reads':str(D/'union.fastq'),'adapters':str(R/'ref/adapters.fasta')},True)
trimmed=one('prepare/published/FASTP/*_fastp.fastq.gz')
nf('protein-run','experiments/amino_acids/protein-parallel.nf',None,{'sample':sample,'reads':str(trimmed),'adapters':str(R/'ref/adapters.fasta'),'validate':False},True)
hits={'production':baseline/f'{sample}_virus_hits.tsv.gz','protein_variant':one('protein-run/published/PROTEIN_ADAPTER/*_virus_hits.tsv.gz')}
for arm,hit in hits.items():
    view=D/arm/'run-results';view.mkdir(parents=True,exist_ok=True)
    for source in baseline.iterdir():
        if source.name==f'{sample}_virus_hits.tsv.gz':continue
        target=view/source.name
        if not target.exists():os.link(source,target)
    target=view/f'{sample}_virus_hits.tsv.gz'
    if not target.exists():os.link(hit,target)
    groups=D/arm/'groups.tsv';groups.write_text(f'sample\tgroup\n{sample}\t{sample}\n')
    manifest=D/arm/'input.csv'
    with manifest.open('w') as f:
        w=csv.writer(f);w.writerow(['label','run_results_dir','groups_tsv']);w.writerow(['cohort',view,groups])
    save(D/arm/'shared-inputs.json',{'profile_and_qc_from':str(baseline),'reason':'Same original input and production profiling configuration; independent of viral alignment','viral_hits_sha256':sha(hit)})
    nf(arm+'-downstream','main.nf','configs/downstream.config',{'platform':'illumina','input_file':str(manifest),'aln_dup_deviation':1,'validation_n_sample':20,'blast_task':None})
left=json.loads((D/'production-downstream/params.json').read_text());right=json.loads((D/'protein_variant-downstream/params.json').read_text())
for k in ('base_dir','input_file'):left.pop(k);right.pop(k)
assert left==right
save(D/'complete.json',{'sample':sample,'source_commit':a.source_commit,'pairs':pairs,'both_arms_complete':True,'downstream_parameters_equal':True,'completed_utc':now().isoformat()})
backup()
