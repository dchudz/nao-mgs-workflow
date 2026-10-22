"""Describe two predeclared mate-support restrictions on frozen similarity groups."""
import argparse,collections,csv,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=False)
def rows(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:yield from csv.DictReader(f,delimiter='\t')
r=a.results;e=r/'protein200-compacted-nt20/pair-evidence.tsv';rp=r/'protein-resolution-compacted-v2/resolution.tsv';mp=r/'protein-similarity-support-compacted-v2-nt20/cluster-membership.tsv.gz';cp=r/'protein-similarity-support-compacted-v2-nt20/sample-species-calls.tsv'
evidence={x['query_pair']:x for x in rows(e)};resolution={x['query_pair']:x for x in rows(rp)};prior=list(rows(cp));b=collections.defaultdict(set);groups=collections.defaultdict(set);seen=set();counts=collections.Counter()
for x in rows(mp):
 sample,cluster=x['sample'],x['similarity_exemplar']
 for sp in x['native_species'].split(';'):
  if sp:b[sample,sp].add(cluster)
 q=x['protein_query']
 if not q:continue
 assert q not in seen;seen.add(q);v=resolution[q];z=evidence[q];assert v['sample']==z['sample']==sample and z['nr_status']=='viral_only'
 policies=['original']
 if z['nr_aligned_mates']=='2':policies.append('NR_both_mates')
 if z['nr_aligned_mates']=='2' and z['viral_reference_aligned_mates']=='2':policies.append('NR_and_viral_reference_both_mates')
 for policy in policies:counts[policy]+=1
 if not v['nr_species']:continue
 tiers=['NR_resolved_species']+(['NR_and_viral_reference_same_species'] if v['species_relation']=='same_species' else [])
 for policy in policies:
  for tier in tiers:groups[policy,tier,sample,v['nr_species']].add(cluster)
assert seen==resolution.keys()
# Characterization: reconstruct every previously reported support count first.
for x in prior:
 key=x['sample'],x['species'];base=b[key];prot=groups['original',x['tier'],*key]
 assert len(base)==int(x['baseline_similarity_clusters']) and len(prot)==int(x['protein_similarity_clusters']) and len(base|prot)==int(x['union_similarity_clusters'])
out=[];summary={}
for policy in ['NR_both_mates','NR_and_viral_reference_both_mates']:
 for tier in ['NR_resolved_species','NR_and_viral_reference_same_species']:
  selected=[]
  for x in prior:
   if x['tier']!=tier:continue
   key=x['sample'],x['species'];base=b[key];prot=groups[policy,tier,*key];gained=len(base)<3<=len(prot);union_gained=len(base)<3<=len(base|prot)
   row=dict(policy=policy,tier=tier,sample=x['sample'],species=x['species'],human_annotation=x['human_annotation'],baseline_support=len(base),protein_support=len(prot),union_support=len(base|prot),gained_protein_call=gained,gained_union_call=union_gained,previous_protein_support=x['protein_similarity_clusters'],previous_gained_protein_call=x['gained_protein_call']);out.append(row);selected.append(row)
  summary[policy+'/'+tier]={'protein_calls':sum(x['protein_support']>=3 for x in selected),'gained_protein_calls':sum(x['gained_protein_call'] for x in selected),'human_annotated_gains':sum(x['gained_protein_call'] and x['human_annotation']=='1' for x in selected),'gains_from_zero':sum(x['gained_protein_call'] and x['baseline_support']==0 for x in selected),'gained_union_calls':sum(x['gained_union_call'] for x in selected)}
with (a.output/'sample-species-calls.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
result={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_hashes':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [e,rp,mp,cp]},'prior_counts_reconstructed_exactly':len(prior),'viral_only_pair_counts':dict(counts),'counts':summary,'original_four_human_gains':[x for x in out if x['tier']=='NR_and_viral_reference_same_species' and x['human_annotation']=='1' and x['previous_gained_protein_call']=='True'],'limits':['Post-hoc descriptive robustness check; thresholds were not changed to maximize gains.','Uses original common similarity partition unchanged, without rerunning clustering after filtering.','Each retained pair has its own viral-only evidence; no verdict propagated through a cluster.','No production alignment-coordinate duplicate stage or biological truth labels.','Native baseline was not symmetrically required to pass full NR validation on both mates.']};(a.output/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps({'counts':summary,'original_four_human_gains':result['original_four_human_gains']}))
