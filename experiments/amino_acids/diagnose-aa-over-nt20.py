"""Attribute AA-union-only final calls to the initial screen or later processing."""
import argparse,collections,csv,gzip,hashlib,json,subprocess,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(r/'cohort-storage/aa-over-nt20'));args=parser.parse_args()
out=Path(args.output);out.mkdir(exist_ok=True)
selected=[x for x in csv.DictReader(open(r/'cohort-storage/analysis/added-pairs.tsv'),delimiter='\t') if x['arm']=='union' and x['also_nt20']=='False'];by_sample=collections.defaultdict(list)
for x in selected:by_sample[x['sample']].append(x)
batches={a:p.parent for p in (r/'cohort-storage/batches').glob('*/complete.json') for a in json.loads(p.read_text())['samples']}
with gzip.open(r/'reference/results/total-virus-db-annotated.tsv.gz','rt') as f:names={x['taxid']:x.get('name','') for x in csv.DictReader(f,delimiter='\t')}
def selected_pairs(p,wanted):
 result={}
 for key,raw in pairs(p):
  q=key[0].decode().lstrip('@')
  if q in wanted:result[q]=key
 return result
result=[]
for a,rows in by_sample.items():
 wanted={x['seq_id'] for x in rows};d=batches[a];raw=out/(a+'_nt20.fastq')
 subprocess.run(['aws','--profile','pipeline-codex','s3','cp',f's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/cohort-screens/{a}/nt20.fastq',str(raw),'--only-show-errors'],check=True)
 initial=selected_pairs(raw,wanted);base=selected_pairs(d/'output/fastp'/f'{a}_nt20_fastp.fastq.gz',wanted);aa=selected_pairs(d/'output/fastp'/f'{a}_union_fastp.fastq.gz',wanted)
 for x in rows:
  q=x['seq_id'];assert q in aa
  row=dict(x);row.update(name=names.get(x['taxid'],''),stage='screen' if q not in initial else 'fastp' if q not in base else 'later_processing',nt20_trimmed_lengths=None if q not in base else [len(base[q][1]),len(base[q][3])],union_trimmed_lengths=[len(aa[q][1]),len(aa[q][3])],trimmed_sequences_equal=None if q not in base else (base[q][1],base[q][3])==(aa[q][1],aa[q][3]),trimmed_qualities_equal=None if q not in base else (base[q][2],base[q][4])==(aa[q][2],aa[q][4]));result.append(row)
 raw.unlink()
summary={'completed_cohort_samples':json.loads((r/'cohort-storage/analysis/summary.json').read_text())['completed_samples'],'selection_table_sha256':hashlib.sha256((r/'cohort-storage/analysis/added-pairs.tsv').read_bytes()).hexdigest(),'added_beyond_nt20':len(result),'samples':len(by_sample),'stage':dict(collections.Counter(x['stage'] for x in result)),'interpretation':'Compare exact original pair IDs at the raw screen and FASTP boundaries. Final nucleotide assignments are a comparator, not known truth.'};(out/'pairs.json').write_text(json.dumps(result,indent=2));(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
