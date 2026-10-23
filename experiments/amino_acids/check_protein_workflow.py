"""Production-module integration on a frozen, enriched software fixture only."""
import csv,gzip,hashlib,json,subprocess,sys
from pathlib import Path
from union_pairs import pairs
R=Path('/home/ec2-user/experiment');D=R/'integration/protein-workflow';D.mkdir(parents=True,exist_ok=True)
source='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/cohort-screens/SRR37076045/aa12.fastq'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',source,str(D/'saved-aa12.fastq'),'--only-show-errors'],check=True)
with gzip.open(D/'fixture.fastq.gz','wb') as out:
    count=0
    for key,raw in pairs(D/'saved-aa12.fastq'):
        if count==1000:break
        out.write(raw);count+=1
config=D/'publish.config';config.write_text("process.publishDir = [path: { \"${params.base_dir}/published/${task.process.tokenize(':')[-1]}\" }, mode:'copy']\n")
def run(name,script,extra):
    stage=D/name;stage.mkdir(exist_ok=True)
    command=[str(R/'nextflow'),'-log',str(stage/'nextflow.log'),'run',str(R),'-main-script',script,
        '-c',str(R/'experiments/amino_acids/e2e-worker.config'),'-c',str(config),'-work-dir',str(stage/'work'),
        '--base_dir',str(stage),'--ref_dir',str(R/'reference'),'--sample','adapter_fixture','--adapters',str(R/'ref/adapters.fasta'),'-ansi-log','false',*extra]
    with (stage/'run.log').open('w') as f:subprocess.run(command,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
    return command
cmd1=run('native','experiments/amino_acids/postscreen.nf',['--reads',str(D/'fixture.fastq.gz'),'--validate','false'])
trimmed=next((D/'native/published/FASTP').glob('*_fastp.fastq.gz'))
cmd2=run('protein','experiments/amino_acids/protein-production.nf',['--reads',str(trimmed)])
hits=D/'protein/published/PROTEIN_ADAPTER/adapter_fixture_virus_hits.tsv.gz'
with gzip.open(hits,'rt') as f:before=list(csv.DictReader(f,delimiter='\t'))
assert before,'The integration fixture must exercise actual candidate rows'
cmd3=run('downstream','experiments/amino_acids/check-protein-downstream.nf',['--hits',str(hits)])
final=next((D/'downstream/published/COPY_HITS').glob('*_validation_hits.tsv.gz'))
with gzip.open(final,'rt') as f:after=list(csv.DictReader(f,delimiter='\t'))
assert {r['seq_id'] for r in before}=={r['seq_id'] for r in after}
by_id={r['seq_id']:r for r in before}
for row in after:
    for key in ['query_seq','query_seq_rev','query_qual','query_qual_rev']:
        assert row[key]==by_id[row['seq_id']][key]
    assert row['validation_status'] in {'aligned','no_alignment','not_sampled'}
receipt={'verified':True,'source':source,'source_sha256':hashlib.sha256((D/'saved-aa12.fastq').read_bytes()).hexdigest(),
    'fixture_rule':'First1000 original pairs of a preserved AA-screen-selected sample, fixed before its new alignments. Enriched software fixture, not an unbiased biological comparison.',
    'input_pairs':count,'protein_rows_preserved_through_downstream':len(after),'commands':[cmd1,cmd2,cmd3],
    'checks':['Real HSPs accepted by adapter','Both original sequences and qualities preserved','Production alignment and similarity duplicate marking completed',
              'Production hash sampling, merged/joined query preparation, sensitive core_nt validation and per-read statuses completed'],
    'biological_interpretation':'None. This fixture tests software integration and is excluded from full-depth cohort results.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
