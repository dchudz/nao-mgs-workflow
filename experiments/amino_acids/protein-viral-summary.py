"""Native protein evidence: near-best taxonomic LCA; no nucleotide score substitution."""
import argparse,csv,gzip,json,sys,collections
import gc
# Full-depth replay verified identical TSV, JSON, and FASTA; avoid repeated scans of retained objects.
gc.disable()
from pathlib import Path
r=Path('/home/ec2-user/experiment')
sys.path.insert(0,str(r/'modules/local/lcaTsv/resources/usr/bin'))
from lca_tsv import parse_nodes_db,find_lca_set,path_to_root
p=argparse.ArgumentParser();p.add_argument('hits');p.add_argument('metadata');p.add_argument('output');a=p.parse_args()
with gzip.open(r/'reference/results/virus-genome-metadata-gid.tsv.gz','rt') as f:genomes={row['genome_id']:int(row['taxid']) for row in csv.DictReader(f,delimiter='\t')}
parents,_=parse_nodes_db(str(r/'reference/results/taxonomy-nodes.dmp'),81077);cache={}
def lca(taxa):
 global cache
 result,cache=find_lca_set(set(taxa),parents,cache);return result
def viral(taxid):
 global cache
 ancestors,cache=path_to_root(taxid,parents,cache);return 10239 in ancestors
hits=collections.defaultdict(list)
with open(a.hits) as f:
 for row in csv.reader(f,delimiter='\t'):
  q,subject,score,evalue,length,identity,qlen,qstart,qend,sstart,send,frame=row
  if int(length)<25:continue
  gid=subject.split('~')[0];taxid=genomes[gid]
  hits[q].append((float(score),taxid,float(identity),int(length),float(evalue)))
query_calls={}
for q,alignments in hits.items():
 best=max(x[0] for x in alignments);near=[x for x in alignments if x[0]>=0.9*best]
 query_calls[q]={'taxid':lca(x[1] for x in near),'bitscore':best,'alignments':len(near),'best_identity':max(x[2] for x in near),'best_length':max(x[3] for x in near),'evalue':min(x[4] for x in near)}
with gzip.open(a.metadata,'rt') as f:meta={row['query_pair']:row for row in csv.DictReader(f,delimiter='\t')}
rows=[]
for name in sorted({q.rsplit('/',1)[0] for q in query_calls}):
 mates=[query_calls[q] for q in [name+'/1',name+'/2'] if q in query_calls];taxid=lca(x['taxid'] for x in mates);m=meta[name]
 rows.append({'query_pair':name,'example_seq_id':m['example_seq_id'],'occurrences':int(m['occurrences']),'protein_taxid_lca':taxid,'viral_lca':viral(taxid),'aligned_mates':len(mates),'best_bitscore':max(x['bitscore'] for x in mates),'minimum_evalue':min(x['evalue'] for x in mates),'maximum_alignment_aa':max(x['best_length'] for x in mates),'maximum_identity_pct':max(x['best_identity'] for x in mates),'broad_validation':'pending'})
with open(a.output,'w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');writer.writeheader();writer.writerows(rows)
summary={'unique_query_pairs':len(meta),'aligned_query_reads':len(query_calls),'aligned_pairs':len(rows),'viral_lca_pairs':sum(x['viral_lca'] for x in rows),'two_mate_pairs':sum(x['aligned_mates']==2 for x in rows),'rules':{'minimum_alignment_amino_acids':25,'near_best_bitscore_fraction':0.9,'pair_taxonomy':'LCA of eligible near-best taxa across both aligned mates','independent_NR_validation':'pending'}}
Path(a.output+'.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
# Include both mates of every candidate pair in the independent broad-reference query.
with open(a.output+'.fna','w') as f:
 for row in rows:
  m=meta[row['query_pair']]
  for mate,key in [(1,'query_seq'),(2,'query_seq_rev')]:f.write(f'>{row["query_pair"]}/{mate}\n{m[key]}\n')
