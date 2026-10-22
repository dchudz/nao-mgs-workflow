"""Join measured mate evidence into pair-level descriptive comparisons, including unchanged NR calls."""
import argparse,collections,csv,hashlib,itertools,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--nr-calls',type=Path,required=True);p.add_argument('--old-panel',type=Path,required=True);p.add_argument('--nodes',type=Path,required=True);a=p.parse_args();d=a.directory;out=d/'pair-analysis';out.mkdir(exist_ok=False)
parents={}
for line in a.nodes.open():
 f=line.split('|');parents[int(f[0])]=int(f[1])
cache={}
def ancestry(t):
 if t not in cache:
  seen=set();path=[];original=t
  while t in parents and t not in seen:
   seen.add(t);path.append(t)
   if parents[t]==t:cache[original]=path;break
   t=parents[t]
  else:cache[original]=None
 return cache[t]
def lca(taxa):
 paths=[ancestry(t) for t in taxa]
 if not paths or any(p is None for p in paths):return 0
 common=set(paths[0]).intersection(*(set(p) for p in paths[1:]));return next(t for t in paths[0] if t in common)
def classify(taxa):
 unknown={t for t in taxa if ancestry(t) is None};known=taxa-unknown;viral={t for t in known if 10239 in ancestry(t)};nonviral=known-viral;artificial={t for t in known if 81077 in ancestry(t)};cellular={t for t in known if any(x in ancestry(t) for x in [131567,2,2157,2759])};natural=known-artificial
 # The raw pair status preserves the frozen NR classifier's precedence. Unknown counts remain explicit.
 status='no_alignment' if not taxa else 'mixed_viral_nonviral' if viral and nonviral else 'unknown_taxonomy' if unknown else 'nonviral_only' if nonviral else 'viral_only'
 nstatus='unknown_taxonomy' if unknown else 'artificial_only' if artificial and not natural else 'no_alignment' if not natural else 'mixed_viral_natural_nonviral' if viral and natural-viral else 'viral_only' if viral else 'natural_nonviral_only'
 return dict(status=status,natural_status=nstatus,taxid_lca=lca(taxa),known_natural_lca=lca(natural),viral_taxa=len(viral),artificial_taxa=len(artificial),cellular_taxa=len(cellular),unknown_taxa=len(unknown),near_best_taxids=';'.join(map(str,sorted(taxa))))
selection={x['seq_id']:x for x in csv.DictReader((d/'selection.tsv').open(),delimiter='\t')};assert len(selection)==232
queries={line[1:].split()[0] for line in (d/'queries.fna').read_text().splitlines() if line.startswith('>')};assert queries=={q+'/'+str(m) for q in selection for m in [1,2]}
calls={};methods=set();thresholds=set()
for row in csv.DictReader((d/'analysis/per-mate.tsv').open(),delimiter='\t'):
 key=(row['method'],row['threshold'],row['query']);assert key not in calls and row['query'] in queries;calls[key]={int(t) for t in row['near_best_taxids'].split(';') if t};methods.add(row['method']);thresholds.add(row['threshold'])
assert len(calls)==len(methods)*len(thresholds)*len(queries)
nr={}
for line in a.nr_calls.open():
 q,score,taxa=line.rstrip('\n').split('\t');assert q in queries and q not in nr;nr[q]={int(t) for t in taxa.split(';')}
old={x['query_pair']:x for x in csv.DictReader((a.old_panel/'comparison.tsv').open(),delimiter='\t') if x['task']=='protein_NR'};assert old.keys()==selection.keys()
rows=[]
for method,threshold in [*itertools.product(sorted(methods),sorted(thresholds)),('protein_NR','nr_native_E1e-5_min25')]:
 for pair,sel in selection.items():
  mates=[nr.get(pair+'/'+str(m),set()) if method=='protein_NR' else calls[method,threshold,pair+'/'+str(m)] for m in [1,2]];taxa=set().union(*mates)
  row=dict(pair=pair,sample=sel['sample'],stratum=sel['stratum'],method=method,threshold=threshold,aligned_mates=sum(bool(x) for x in mates),**classify(taxa));rows.append(row)
  if method=='protein_NR':assert all(str(row[k])==old[pair][k] for k in ['status','taxid_lca','aligned_mates']),pair
with (out/'per-pair.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
lookup={(x['method'],x['threshold'],x['pair']):x for x in rows};groups={};changes=[];transitions={};strata=sorted({x['stratum'] for x in rows})
for method,threshold in sorted({(x['method'],x['threshold']) for x in rows}):
 for stratum in strata:
  group=[x for x in rows if x['method']==method and x['threshold']==threshold and x['stratum']==stratum];groups[method+'/'+threshold+'/'+stratum]={'pairs':len(group),'status':dict(collections.Counter(x['status'] for x in group)),'natural_status':dict(collections.Counter(x['natural_status'] for x in group))}
comparisons=[(l,t,r,t) for t in sorted(thresholds) for l,r in itertools.combinations(sorted(methods),2)]+[(m,'strict_E1e-5','protein_NR','nr_native_E1e-5_min25') for m in sorted(methods)]
for left,lt,right,rt in comparisons:
 for stratum in strata:
  counts=collections.Counter();n=0
  for pair,sel in selection.items():
   if sel['stratum']!=stratum:continue
   l=lookup[left,lt,pair];r=lookup[right,rt,pair];n+=1;counts[l['status']+' -> '+r['status']]+=1
   if any(l[k]!=r[k] for k in ['status','taxid_lca','natural_status','known_natural_lca']):changes.append(dict(pair=pair,stratum=stratum,left=left,left_threshold=lt,right=right,right_threshold=rt,**{'left_'+k:l[k] for k in ['status','taxid_lca','natural_status','known_natural_lca']},**{'right_'+k:r[k] for k in ['status','taxid_lca','natural_status','known_natural_lca']}))
  transitions[left+'/'+lt+' -> '+right+'/'+rt+'/'+stratum]={'pairs':n,'status':dict(counts)}
with (out/'disagreements.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(changes[0]) if changes else ['pair'],delimiter='\t');w.writeheader();w.writerows(changes)
result={'pairs':232,'mates':464,'groups':groups,'transitions':transitions,'NR_frozen_pair_verdict_parity_verified':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'analysis/per-mate.tsv',a.nr_calls,a.old_panel/'comparison.tsv',a.nodes]},'limits':['Frozen sample-balanced selected panel, not prevalence or known-truth accuracy.','Core_nt mate taxon sets use production subject deduplication,10 dense score ranks and90percent near-best. Pair status uses the frozen NR precedence; unknown counts remain explicit.','NR keeps its original E1e-5,minimum25AA and90percent score filter, with no10-rank restriction. References, programs and filters differ.','NR pair verdicts are checked exactly against the completed original NR panel; no earlier verdict is replaced.','Known natural LCA omits artificial and unknown taxa; natural_status and unknown_taxa preserve that uncertainty.']};(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['groups','transitions','input_hashes']}))
