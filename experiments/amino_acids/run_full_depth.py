"""Execute the frozen ten-sample comparison with per-stage private checkpoints.

Run from the repository root, with one selected sample first. All four outputs
share full-depth raw input. Production RUN QC/profile outputs are reused because
those branches do not depend on viral screening or alignment.
"""
import argparse,csv,datetime,gzip,hashlib,json,os,shutil,subprocess,tarfile,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');E=R/'experiments/amino_acids';PLAN=E/'end-to-end-full-depth'
AWS=['aws','--profile','pipeline-codex','--region','us-east-1'];PREFIX='s3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/'
p=argparse.ArgumentParser();p.add_argument('--sample',required=True);p.add_argument('--source-commit',required=True);a=p.parse_args()
assert len(a.source_commit)==40 and all(c in '0123456789abcdef' for c in a.source_commit)
assert (R/'source-commit.txt').read_text().strip()==a.source_commit
selection=json.loads((PLAN/'sample-selection.json').read_text())['initial_samples'];row=next(x for x in selection if x['run_accession']==a.sample)
D=R/'e2e'/a.sample;D.mkdir(parents=True,exist_ok=True);data=R/'data'/a.sample
bound=json.loads((R/'execution-bound.json').read_text());deadline=datetime.datetime.fromisoformat(bound['absolute_stop_utc'])-datetime.timedelta(hours=2)
def now():return datetime.datetime.now(datetime.timezone.utc)
def save(path,obj):Path(path).write_text(json.dumps(obj,indent=2)+'\n')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
def backup():
    # Only generated work directories and public raw reads are excluded. Published
    # candidates, real HSPs, read-level outputs, commands and resource traces remain.
    for work in D.glob('*/work'):
        diagnostic=work.parent/'task-diagnostics.tar.gz'
        # Rebuild only when this stage's work tree changes.
        tasks=list(work.glob('*/*/.command.*'))+list(work.glob('*/*/.exitcode'))
        newest=max((f.stat().st_mtime for f in tasks if f.is_file()),default=0)
        if not diagnostic.exists() or diagnostic.stat().st_mtime<newest:
            with tarfile.open(diagnostic,'w:gz') as tar:
                for f in tasks:
                    if f.is_file() and not f.is_symlink():tar.add(f,arcname=str(f.relative_to(work)))
    cmd=AWS+['s3','sync',str(D)+'/',PREFIX+f'samples/{a.sample}/','--exclude','*/work/*','--exclude','*/.nextflow/*','--exclude','*.partial','--only-show-errors']
    subprocess.run(cmd,check=True,timeout=1800)
def run(stage,cmd,cwd=R):
    receipt=D/(stage+'.json')
    if receipt.exists():
        old=json.loads(receipt.read_text())
        assert old['command']==list(map(str,cmd)),('Completed stage command changed',stage)
        assert old['source_commit']==a.source_commit,('Source changed after stage completion',stage)
        return
    remaining=(deadline-now()).total_seconds();assert remaining>1800,'Preserve/stop budget reserve reached'
    started=now();begin=time.monotonic()
    with (D/(stage+'.log')).open('w') as log:
        proc=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=log,stderr=subprocess.STDOUT,timeout=remaining-900)
    record={'stage':stage,'source_commit':a.source_commit,'command':list(map(str,cmd)),
            'started_utc':started.isoformat(),'finished_utc':now().isoformat(),'wall_seconds':time.monotonic()-begin,'exit_code':proc.returncode}
    if proc.returncode:
        save(D/(stage+'.failure.json'),record);backup();raise RuntimeError(f'{stage} failed; inspect its preserved log')
    save(receipt,record);backup();print(json.dumps(record),flush=True)
def nf(stage,script,config,params):
    target=D/stage;target.mkdir(exist_ok=True)
    parameters={**params,'base_dir':str(target),'ref_dir':str(R/'reference')}
    save(target/'params.json',parameters)
    cmd=[R/'nextflow','-log',target/'nextflow.log','run',R,'-main-script',script,'-lib',str(R/'lib'),
         '-c',R/config,'-c',E/'e2e-worker.config','-params-file',target/'params.json','-work-dir',target/'work','-ansi-log','false']
    run(stage,cmd)
def single(pattern):
    paths=list(D.glob(pattern));assert len(paths)==1,(pattern,paths);return paths[0]
def counts(path):
    with open(path) as f:return next(csv.DictReader(f,delimiter='\t'))
def results_view(arm,hits):
    dest=D/arm/'run-results';dest.mkdir(parents=True,exist_ok=True)
    baseline=D/'baseline-run/output/results'
    for source in baseline.iterdir():
        if source.name==f'{a.sample}_virus_hits.tsv.gz':continue
        target=dest/source.name
        if not target.exists():os.link(source,target)
    target=dest/f'{a.sample}_virus_hits.tsv.gz'
    if not target.exists():os.link(hits,target)
    save(D/arm/'shared-run-output-provenance.json',{'shared_from':str(baseline),
         'reason':'Production QC and profiling consume the same full-depth raw input and fixed profile subsample, independently of viral alignment.',
         'viral_hits':str(hits),'viral_hits_sha256':sha(hits),'grouping':'Each original sample is its own downstream group in every arm.'})
    groups=D/arm/'groups.tsv';groups.write_text(f'sample\tgroup\n{a.sample}\t{a.sample}\n')
    manifest=D/arm/'input.csv'
    with manifest.open('w') as f:
        w=csv.writer(f);w.writerow(['label','run_results_dir','groups_tsv']);w.writerow(['cohort',dest,groups])
    return manifest

assert json.loads(subprocess.check_output(AWS+['sts','get-caller-identity'],text=True))['Account']=='839993661825'
assert json.loads((R/'setup/complete.json').read_text())['verified']
assert json.loads((R/'integration/real-diamond/verification.json').read_text())['verified']
assert json.loads((R/'integration/baseline-replay.json').read_text())['verified']
assert json.loads((R/'integration/protein-workflow/verification.json').read_text())['verified']
assert json.loads((R/'integration/adapter-bounded/verification.json').read_text())['verified']
acquisition=json.loads((data/'verified.json').read_text());assert acquisition['spot_count']==int(row['read_count'])
assert not (D/'complete.json').exists(),'Sample already completed; do not silently rerun it'
save(D/'protocol.json',{'source_commit':a.source_commit,'sample':a.sample,'selection':row,
    'started_utc':now().isoformat(),'acquisition':acquisition,'logical_outputs':json.loads((PLAN/'protocol.json').read_text())['logical_outputs'],
    'grouping':'One original sample per group; 20 hash-selected joint duplicate exemplars per species per group.',
    'validation':'Original production merged/joined nucleotide query preparation; individual read results, never extrapolated.'})
sheet=D/'samplesheet.csv'
with sheet.open('w') as f:
    w=csv.writer(f);w.writerow(['sample','fastq_1','fastq_2']);w.writerow([a.sample,data/'R1.fastq.gz',data/'R2.fastq.gz'])
nf('baseline-run','main.nf','configs/run.config',{'platform':'illumina','sample_sheet':str(sheet),'adapters':str(R/'ref/adapters.fasta'),'n_reads_profile':1000000,'random_seed':'17310','bt2_score_threshold':20,'bracken_threshold':1,'taxid_artificial':'81077'})
measured=counts(D/'baseline-run/output/results'/f'{a.sample}_read_counts.tsv')
assert int(float(measured['n_read_pairs']))==acquisition['spot_count'],measured
run('aa-screen',[R/'rust-tools/target/release/screen_amino_acids','screen','--index',R/'reference/virus.aa12.bin',
                 '--reads',data/'R1.fastq.gz','--reads2',data/'R2.fastq.gz','--threads','32','--matched',D/'aa12.fastq'])
native=single('baseline-run/output/intermediates/reads/raw_viral/*_nucleaze_match.fastq.gz')
run('union',['python3.12',E/'union_pairs.py',native,D/'aa12.fastq',D/'union.fastq',D/'rescued.fastq',D/'union-counts.json'])
# Publish experimental stage outputs using their original module/process names.
extra_config=D/'extra.config';extra_config.write_text("process.publishDir = [path: { \"${params.base_dir}/published/${task.process.tokenize(':')[-1]}\" }, mode: 'copy']\n")
def extra(stage,script,reads):
    target=D/stage;target.mkdir(exist_ok=True)
    save(target/'params.json',{'base_dir':str(target),'ref_dir':str(R/'reference'),'sample':a.sample,'reads':str(reads),'adapters':str(R/'ref/adapters.fasta'),'validate':False})
    run(stage,[R/'nextflow','-log',target/'nextflow.log','run',R,'-main-script',script,'-lib',str(R/'lib'),'-c',E/'e2e-worker.config','-c',extra_config,
               '-params-file',target/'params.json','-work-dir',target/'work','-ansi-log','false'])
extra('screen-run','experiments/amino_acids/postscreen.nf',D/'union.fastq')
trimmed=single('screen-run/published/FASTP/*_fastp.fastq.gz')
extra('protein-run','experiments/amino_acids/protein-production.nf',trimmed)
hit_paths={'production':D/'baseline-run/output/results'/f'{a.sample}_virus_hits.tsv.gz',
           'validation_control':D/'baseline-run/output/results'/f'{a.sample}_virus_hits.tsv.gz',
           'screen_control':single('screen-run/published/RENAME_VIRUS_HITS/*_virus_hits.tsv.gz'),
           'protein_variant':single('protein-run/published/PROTEIN_ADAPTER/*_virus_hits.tsv.gz')}
for arm,hits in hit_paths.items():
    manifest=results_view(arm,hits)
    nf(arm+'-downstream','main.nf','configs/downstream.config',{'platform':'illumina','input_file':str(manifest),
       'aln_dup_deviation':1,'validation_n_sample':20,'blast_task':None if arm=='production' else 'blastn'})
# The validator contrast must have identical candidate evidence and sampling.
base=D/'production-downstream/output/results_downstream'/f'{a.sample}_validation_hits.tsv.gz'
control=D/'validation_control-downstream/output/results_downstream'/f'{a.sample}_validation_hits.tsv.gz'
with gzip.open(base,'rt') as left,gzip.open(control,'rt') as right:
    l=csv.DictReader(left,delimiter='\t');r=csv.DictReader(right,delimiter='\t')
    keys=[k for k in l.fieldnames if not k.startswith('validation_')]
    import itertools
    checked=0
    for x,y in itertools.zip_longest(l,r):
        assert x is not None and y is not None
        assert all(x[k]==y[k] for k in keys),('Validator contrast changed upstream evidence',x['seq_id'])
        checked+=1
save(D/'validation-control-parity.json',{'verified':True,'rows':checked,'compared_columns':keys})
save(D/'complete.json',{'sample':a.sample,'source_commit':a.source_commit,'full_depth_pairs':acquisition['spot_count'],
     'completed_utc':now().isoformat(),'all_four_outputs_completed':True,'budget_stop_utc':deadline.isoformat()})
backup()
