"""Classify broad protein competition, retaining missing-taxonomy uncertainty."""
import argparse,collections,csv,json
from pathlib import Path

def summarize(hits,parents,viral_calls):
 cache={}
 def lineage(t):
  if t not in parents:return None
  original=t;result=[];seen=set()
  while t not in seen:
   if t not in parents:return None
   result.append(t);seen.add(t)
   if parents[t]==t:break
   t=parents[t]
  cache[original]=tuple(result)
  return cache[original]
 def ancestors(t):
  if t in cache:return cache[t]
  return lineage(t)
 def lca(taxa):
  paths=[ancestors(t) for t in taxa]
  if not paths or any(x is None for x in paths):return 0
  common=set(paths[0]).intersection(*(set(x) for x in paths[1:]))
  return next((t for t in paths[0] if t in common),0)
 queries=collections.defaultdict(list)
 for row in hits:
  q,subject,score,evalue,length,identity,taxids=row
  if int(length)<25:continue
  taxa={int(x) for x in taxids.split(';') if x.isdigit()} or {0}
  queries[q].append((float(score),taxa))
 calls={}
 for q,alignments in queries.items():
  best=max(score for score,taxa in alignments)
  taxa=set().union(*(taxa for score,taxa in alignments if score>=0.9*best))
  calls[q]=(taxa,best)
 result=[]
 for pair,original_taxid in viral_calls.items():
  mates=[calls[q] for q in (pair+'/1',pair+'/2') if q in calls]
  taxa=set().union(*(x[0] for x in mates)) if mates else set()
  unknown={t for t in taxa if ancestors(t) is None}
  known=taxa-unknown;viral={t for t in known if 10239 in ancestors(t)};nonviral=known-viral
  if not mates:status='no_alignment'
  elif nonviral and viral:status='mixed_viral_nonviral'
  elif unknown:status='unknown_taxonomy'
  elif nonviral:status='nonviral_only'
  else:status='viral_only'
  taxid=lca(known) if not unknown else 0
  original=int(original_taxid);a=ancestors(original);b=ancestors(taxid)
  agreement='unknown' if a is None or b is None else ('exact' if original==taxid else 'ancestor_related' if original in b or taxid in a else 'divergent')
  result.append(dict(query_pair=pair,viral_reference_taxid=original,nr_status=status,nr_taxid_lca=taxid,taxonomic_relation=agreement,aligned_mates=len(mates),viral_taxa=len(viral),nonviral_taxa=len(nonviral),unknown_taxa=len(unknown),best_bitscore=max((x[1] for x in mates),default=0)))
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--nodes',required=True);p.add_argument('--hits',required=True);p.add_argument('--viral-calls',nargs='+',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 parents={}
 with open(a.nodes) as f:
  for line in f:
   fields=line.split('|');parents[int(fields[0])]=int(fields[1])
 viral_calls={}
 for path in a.viral_calls:
  with open(path) as f:
   for row in csv.DictReader(f,delimiter='\t'):
    assert row['query_pair'] not in viral_calls,'Query collision'
    viral_calls[row['query_pair']]=int(row['protein_taxid_lca'])
 with open(a.hits) as f:result=summarize(csv.reader(f,delimiter='\t'),parents,viral_calls)
 with open(a.output,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
 summary={'pairs':len(result),'status':dict(collections.Counter(x['nr_status'] for x in result)),'relation':dict(collections.Counter(x['taxonomic_relation'] for x in result)),'interpretation':'Broad-reference evidence, not known-truth validation; no alignment is unvalidated. Missing or unrecognized taxonomy remains unknown.'}
 Path(a.output+'.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
