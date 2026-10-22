"""Test whether original FASTP processing creates an NT20 match absent in raw pairs."""
import argparse,collections,csv,gzip,json,subprocess,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(r/'cohort-storage/aa-over-nt20'));args=parser.parse_args()
out=Path(args.output);selected=json.loads((out/'pairs.json').read_text());by=collections.defaultdict(set)
for x in selected:by[x['sample']].add(x['seq_id'])
batches={a:p.parent for p in (r/'cohort-storage/batches').glob('*/complete.json') for a in json.loads(p.read_text())['samples']};raws={};trimmed={}
for a,wanted in by.items():
 p=out/(a+'_aa12.fastq');subprocess.run(['aws','--profile','pipeline-codex','s3','cp',f's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/cohort-screens/{a}/aa12.fastq',str(p),'--only-show-errors'],check=True)
 for source,dest in [(p,raws),(batches[a]/'output/fastp'/f'{a}_union_fastp.fastq.gz',trimmed)]:
  for key,data in pairs(source):
   q=key[0].decode().lstrip('@')
   if q in wanted:dest[q]=(key,data)
 p.unlink()
assert set(raws)==set(trimmed)=={x['seq_id'] for x in selected}
results={}
for name,data in [('raw',raws),('trimmed',trimmed)]:
 p=out/(name+'-diagnostic.fastq');p.write_bytes(b''.join(data[q][1] for q in sorted(data)));hit=out/(name+'-nt20.fastq')
 cmd=[str(r/'tools/nucleaze'),'--in',str(p),'--interinput','--outm',str(hit),'--binref',str(r/'reference/virus.nt20.bin'),'--k','20','--canonical','--minhits','1','--threads','2']
 with (out/(name+'-nt20.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
 results[name]={key[0].decode().lstrip('@') for key,data in pairs(hit)}
rows=[]
for q in sorted(raws):
 distances=[]
 for i in [1,3]:
  a=raws[q][0][i];b=trimmed[q][0][i];assert len(b)<=len(a)
  distances.append(min(sum(x!=y for x,y in zip(a[start:start+len(b)],b)) for start in range(len(a)-len(b)+1)))
 rows.append({'seq_id':q,'raw_nt20_match':q in results['raw'],'after_FASTP_nt20_match':q in results['trimmed'],'minimum_substitutions_after_end_trimming_per_mate':distances})
summary={'pairs':rows,'interpretation':'A new exact20-base match after FASTP, when absent in the same raw pair, cannot be created by end trimming alone. Original FASTP includes paired-read correction; this is an isolated diagnostic, not a changed primary pipeline.'};(out/'correction-diagnostic.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
