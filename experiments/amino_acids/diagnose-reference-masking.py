"""Inspect exact nucleotide seeds in the assigned reference before/after human masking."""
import argparse,csv,gzip,json,re,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
parser=argparse.ArgumentParser();parser.add_argument('--output',default=str(r/'cohort-storage/aa-over-nt20'));args=parser.parse_args()
out=Path(args.output);wanted={x['seq_id'] for x in json.loads((out/'pairs.json').read_text())};hits={}
for p in (r/'cohort-storage/batches').glob('*/output/hits/*_union_virus_hits.tsv.gz'):
 with gzip.open(p,'rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   if x['seq_id'] in wanted:hits[x['seq_id']]=x
assert set(hits)==wanted
ids={q:set(re.findall(r'[A-Za-z]+_?\d+\.\d+',x['prim_align_genome_id_all'])) for q,x in hits.items()};gids=set().union(*ids.values());assert all(ids.values())
def read_selected(path):
 result={};name=None;parts=[]
 with gzip.open(path,'rt') as f:
  for line in f:
   if line.startswith('>'):
    if name in gids:result[name]=''.join(parts).upper();parts=[]
    name=line[1:].split()[0]
   elif name in gids:parts.append(line.strip())
  if name in gids:result[name]=''.join(parts).upper()
 assert set(result)==gids,(path,gids-set(result));return result
refs={label:read_selected(r/path) for label,path in [('screen_human_masked','reference/virus-genomes-human-masked.fasta.gz'),('alignment_DUST_masked','reference/results/virus-genomes-masked.fasta.gz')]}
queries={key[0].decode().lstrip('@'):key for key,data in pairs(out/'trimmed-diagnostic.fastq')}
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def longest(q,ref):
 lo=0;hi=len(q)
 while lo<hi:
  k=(lo+hi+1)//2
  found=any(q[i:i+k] in ref or rc(q[i:i+k]) in ref for i in range(len(q)-k+1))
  if found:lo=k
  else:hi=k-1
 return lo
result=[]
for q in sorted(wanted):
 for gid in sorted(ids[q]):
  row={'seq_id':q,'assigned_genome':gid,'edit_distance':hits[q]['prim_align_edit_distance'],'edit_distance_rev':hits[q]['prim_align_edit_distance_rev']}
  for label,records in refs.items():row[label+'_longest_exact_nt_per_mate']=[longest(queries[q][i].decode(),records[gid]) for i in [1,3]]
  result.append(row)
summary={'pairs':result,'interpretation':'Longest exact nucleotide substring in either strand of the reported primary reference genome. The screen uses human-masked references; Bowtie2 uses DUST-masked references. This is a local reference diagnostic, not a change to the frozen comparison.'};(out/'reference-masking-diagnostic.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
