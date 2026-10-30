"""Summarize each read's own original validation evidence without extrapolation."""
import argparse,collections,csv,gzip,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('validation');p.add_argument('--nodes',required=True);p.add_argument('--arms',nargs='+',default=['nt','nt22','nt20','aa12','union']);a=p.parse_args()
d=Path(a.validation);parents={}
for line in open(a.nodes):
 f=line.split('|');parents[int(f[0])]=int(f[1])
cache={}
def path(t):
 if t in cache:return cache[t]
 original=t;seen=set();result=[]
 while t in parents and t not in seen:
  seen.add(t);result.append(t)
  if parents[t]==t:cache[original]=result;return result
  t=parents[t]
 cache[original]=None;return None
def rows(arm):
 with gzip.open(d/'output/COPY_HITS'/f'{arm}_validation_hits.tsv.gz','rt') as f:yield from csv.DictReader(f,delimiter='\t')
baseline={x['seq_id'] for x in rows(a.arms[0])};summaries={};evidence=[]
for arm in a.arms:
 groups={'all':collections.Counter(),'added_to_baseline':collections.Counter()};seen=set()
 for x in rows(arm):
  q=x['seq_id'];assert q not in seen;seen.add(q);status=x['validation_status'];assert status in ['aligned','not_sampled','no_alignment']
  natural=x['validation_staxid_lca'];artificial=x['validation_staxid_lca_artificial'];category=status;relation='unvalidated'
  if status=='aligned':
   lineage=path(int(natural)) if natural not in ['', 'NA'] else None
   if natural in ['', 'NA']:category='aligned_no_natural_LCA'
   elif lineage is None:category='aligned_unknown_natural_taxonomy'
   else:category='aligned_natural_LCA_viral' if 10239 in lineage else 'aligned_natural_LCA_outside_viral_clade'
   original=path(int(x['aligner_taxid_lca']))
   if lineage is not None and original is not None:
    relation='exact' if natural==x['aligner_taxid_lca'] else 'ancestor_related' if int(natural) in original or int(x['aligner_taxid_lca']) in lineage else 'divergent'
   else:relation='unknown'
  for group in ['all']+(['added_to_baseline'] if q not in baseline else []):
   groups[group]['pairs']+=1;groups[group]['status/'+status]+=1;groups[group]['evidence/'+category]+=1
   if status=='aligned':groups[group]['relation/'+relation]+=1
   if artificial not in ['', 'NA']:groups[group]['also_has_artificial_LCA']+=1
  if status!='not_sampled':evidence.append(dict(arm=arm,seq_id=q,added_to_baseline=q not in baseline,status=status,natural_LCA=natural,artificial_LCA=artificial,natural_evidence=category,relation=relation))
 summaries[arm]={k:dict(v) for k,v in groups.items()}
with (d/'sampled-evidence.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(evidence[0]) if evidence else ['arm'],delimiter='\t');w.writeheader();w.writerows(evidence)
summary={'arms':summaries,'interpretation':['Original production joined/merged-mate validation, with 20 selected exemplars per species per pooled arm.','Only each read’s own result is summarized; unsampled reads remain unvalidated.','An LCA outside the viral clade can represent conflicting viral and nonviral evidence; it is not necessarily nonviral-only.','Artificial-sequence evidence is reported separately and can overlap natural evidence.','Species-balanced selection is not an unbiased estimate of read-level accuracy.']}
(d/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
