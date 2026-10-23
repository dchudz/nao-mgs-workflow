"""Replay baseline BLAST module before/after optional task plumbing."""
import gzip,hashlib,json,random,subprocess
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/blast-task';D.mkdir(parents=True,exist_ok=True)
db=D/'db';db.mkdir(exist_ok=True);rng=random.Random(920)
seq=''.join(rng.choices('ACGT',k=1000))
(D/'reference.fna').write_text('>reference\n'+seq+'\n')
q=seq[200:400];mut=''.join(base if i%15 else {'A':'C','C':'G','G':'T','T':'A'}[base] for i,base in enumerate(q))
(D/'query.fna').write_text('>exact\n'+q+'\n>divergent\n'+mut+'\n')
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast:b7be9dc397b73870'
with (D/'makeblastdb.log').open('w') as f:
    subprocess.run(['docker','run','--rm','-v',str(R)+':'+str(R),image,'makeblastdb','-in',str(D/'reference.fna'),'-dbtype','nucl','-parse_seqids','-taxid','10239','-out',str(db/'blast_db')],stdout=f,stderr=subprocess.STDOUT,check=True)
(D/'publish.config').write_text("process.publishDir = [path: { \"${params.base_dir}/published/${task.process.tokenize(':')[-1]}\" }, mode:'copy']\nprocess { withLabel: blast_resources { cpus=4; memory='4 GB' } }\n")
cmd=[str(R/'nextflow'),'-log',str(D/'nextflow.log'),'run',str(R),'-main-script','experiments/amino_acids/check-blast-task.nf',
     '-c',str(R/'experiments/amino_acids/e2e-worker.config'),'-c',str(D/'publish.config'),'-work-dir',str(D/'work'),
     '--base_dir',str(D),'--query',str(D/'query.fna'),'--database',str(db),'-ansi-log','false']
with (D/'run.log').open('w') as f:subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
outputs={k:gzip.open(D/'published'/k/'fixture_hits.blast.gz','rb').read() for k in ['ORIGINAL','DEFAULT','SENSITIVE']}
assert outputs['ORIGINAL'] and outputs['ORIGINAL']==outputs['DEFAULT']
def names(data):return {line.split(b'\t')[0].decode() for line in data.splitlines()}
assert names(outputs['DEFAULT'])=={'exact'},names(outputs['DEFAULT'])
assert names(outputs['SENSITIVE'])=={'exact','divergent'},names(outputs['SENSITIVE'])
receipt={'verified':True,'baseline_module_commit':'fa7cd464','command':cmd,
    'decompressed_output_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in outputs.items()},
    'query_matches':{k:sorted(names(v)) for k,v in outputs.items()},
    'scope':'Actual BLAST execution: absent task option reproduces original module bytes; explicit blastn recovers synthetic dispersed substitutions. Baseline RUN wiring remains unchanged. This is software verification, not biological evidence.'}
(R/'integration/baseline-replay.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
