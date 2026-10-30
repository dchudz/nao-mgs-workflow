"""Aggregate frozen classifications and memberships; no sequence access or reclustering."""
import argparse, collections, csv, gzip, hashlib, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--annotation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def rows(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f: yield from csv.DictReader(f,delimiter='\t')
def sha(p): return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
r=a.results
files={'evidence':r/'protein200-compacted-nt20/pair-evidence.tsv','resolution':r/'protein-resolution-compacted-v2/resolution.tsv','membership':r/'protein-similarity-support-compacted-v2-nt20/cluster-membership.tsv.gz','prior':r/'protein-similarity-support-compacted-v2-nt20/sample-species-calls.tsv','annotation':a.annotation}
e={x['query_pair']:x for x in rows(files['evidence'])};v={x['query_pair']:x for x in rows(files['resolution'])};prior=list(rows(files['prior']))
ann={x['taxid']:x for x in rows(a.annotation)}
human={x['species']:x['human_annotation'] for x in prior}
for sp,x in ann.items():
 if x['rank']=='species':
  if sp in human: assert human[sp]==x['infection_status_human']
  human[sp]=x['infection_status_human']
b=collections.defaultdict(set);g=collections.defaultdict(set);unresolved=collections.defaultdict(set);nr_rank_pairs=collections.Counter();seen=set();paircounts=collections.Counter()
for x in rows(files['membership']):
 sample,cl=x['sample'],x['similarity_exemplar']
 for sp in x['native_species'].split(';'):
  if sp: b[sample,sp].add(cl)
 q=x['protein_query']
 if not q:continue
 assert q not in seen;seen.add(q);z=e[q];t=v[q]
 assert z['nr_status']=='viral_only' and z['sample']==t['sample']==sample
 policies=['pair_viral_only']+(['both_mates_viral_only'] if z['nr_aligned_mates']=='2' else [])
 for policy in policies:
  paircounts[policy]+=1;nr_rank_pairs[policy,t['nr_rank']]+=1
  labels={'NR_species':t['nr_species'],'reference_species':t['reference_species'],'species_agreement':t['nr_species'] if t['species_relation']=='same_species' else ''}
  for tier,sp in labels.items():
   if sp:g[policy,tier,sample,sp].add(cl)
  if not t['reference_species']:
   tax=t['viral_reference_taxid'];rank=ann.get(tax,{}).get('rank','unknown')
   unresolved[policy,rank,sample,tax].add(cl)
assert seen==v.keys()
for x in prior:
 tier={'NR_resolved_species':'NR_species','NR_and_viral_reference_same_species':'species_agreement'}[x['tier']]
 key=x['sample'],x['species'];prot=g['pair_viral_only',tier,*key]
 assert len(b[key])==int(x['baseline_similarity_clusters'])
 assert len(prot)==int(x['protein_similarity_clusters'])
 assert len(b[key]|prot)==int(x['union_similarity_clusters'])
summary=[]
for policy in ['pair_viral_only','both_mates_viral_only']:
 for tier in ['NR_species','species_agreement','reference_species']:
  keys=set(b)|{(s,sp) for pol,ti,s,sp in g if pol==policy and ti==tier}
  calls=gains=zero=uniongains=losses=0;annotations=collections.Counter()
  for key in keys:
   base=b[key];prot=g[policy,tier,*key]
   calls+=len(prot)>=3;losses+=len(base)>=3>len(prot)
   if len(base)<3<=len(prot):
    gains+=1;zero+=not base;annotations[human.get(key[1],'NA')]+=1
   uniongains+=len(base)<3<=len(base|prot)
  summary.append(dict(evidence_policy=policy,taxonomy_policy=tier,protein_calls=calls,gained_calls=gains,gains_from_zero=zero,human_annotation_of_gains=dict(annotations),union_gains=uniongains,descriptive_losses=losses))
assert next(x for x in summary if x['evidence_policy']=='pair_viral_only' and x['taxonomy_policy']=='NR_species')['gained_calls']==701
strict=next(x for x in summary if x['evidence_policy']=='both_mates_viral_only' and x['taxonomy_policy']=='species_agreement')
assert (strict['protein_calls'],strict['gained_calls'],strict['gains_from_zero'],strict['human_annotation_of_gains']['1'])==(622,53,21,2)
# Broad assignments remain broad. Native membership stores species only, so no broad-rank gain is inferred.
broad=[]
for policy in ['pair_viral_only','both_mates_viral_only']:
 for rank in sorted({k[1] for k in unresolved if k[0]==policy}):
  selected=[(s,t,cs) for (pol,rk,s,t),cs in unresolved.items() if pol==policy and rk==rank]
  broad.append(dict(evidence_policy=policy,reference_rank=rank,sample_taxon_groups=len(selected),groups_with_at_least_three_frozen_clusters=sum(len(cs)>=3 for _,_,cs in selected),distinct_taxa=len({t for _,t,_ in selected}),gains_vs_native=None))
result={'scope':'Aggregate bookkeeping of frozen classifications and common similarity memberships only. No sequences, searches, new alignments, candidate lists, or duplicate reruns.','source_sha256':sha(Path(__file__)),'inputs':{k:{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for k,p in files.items()},'prior_rows_reconstructed':len(prior),'threshold_frozen_similarity_clusters':3,'baseline_species_calls':sum(len(x)>=3 for x in b.values()),'counts':summary,'pair_counts':dict(paircounts),'NR_rank_counts':[{'evidence_policy':p,'rank':rank,'pairs':n} for (p,rank),n in sorted(nr_rank_pairs.items())],'reference_above_species':broad,'human_annotation_definition':{'1':'reference annotation compatible with human infection','0':'reference annotation inconsistent with human infection','2':'unclear','NA':'unavailable'},'limits':['Reference-only species labels do not upgrade NR evidence to species resolution. Both labels remain distinct.','All rows use the saved common similarity partition; neither arm is recounted with the production alignment-coordinate stage.','The native arm is NT20 retained output, not the production NT24 screen and not symmetrically NR validated.','Viral-only status depends on reported hits and filtering; NR and core_nt validation are not interchangeable.','Losses are descriptive output differences under asymmetric validation, not estimates of false negatives.','Broad-rank assignments are retained in rank summaries; native saved membership has species labels only, so symmetric broad-rank gains cannot be reconstructed here.','Counts describe the selected cohort, not sensitivity, specificity, infection, or prevalence.']}
a.output.mkdir(parents=True,exist_ok=False);(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'counts':summary,'reference_above_species':broad},indent=2))
