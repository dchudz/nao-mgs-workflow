"""Separate artificial constructs from cellular competitors without changing primary NR classifications."""
import collections,csv,importlib.util,json
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/combined-validation'
spec=importlib.util.spec_from_file_location('nr_summary',r/'protein-nr-summary.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
parents={}
for line in open(r/'protein/current-nr/nodes.dmp'):
 f=line.split('|');parents[int(f[0])]=int(f[1])
cache={}
def category(t):
 if t in cache:return cache[t]
 original=t;seen=set();found='other_known'
 if t not in parents:found='unknown'
 while t in parents and t not in seen:
  if t==10239:found='viral';break
  if t==81077:found='artificial';break
  if t in [131567,2,2157,2759]:found='cellular';break
  seen.add(t);t=parents[t]
 if t not in parents:found='unknown'
 cache[original]=found;return found
from nr_evidence_inputs import near_best_calls
calls=near_best_calls(d,mod)
rows=[]
for source,scope in [('validated-pairs.tsv','protein_candidates'),('native-control-validation.tsv','native_controls')]:
 for x in csv.DictReader(open(d/source),delimiter='\t'):
  q=x['query_pair'];taxa=set().union(*(calls.get(q+'/'+str(m),(set(),0))[0] for m in [1,2]));counts=collections.Counter(category(t) for t in taxa)
  assert counts['viral']==int(x['viral_taxa']) and counts['unknown']==int(x['unknown_taxa'])
  assert sum(counts[k] for k in ['artificial','cellular','other_known'])==int(x['nonviral_taxa'])
  rows.append(dict(query_pair=q,scope=scope,primary_nr_status=x['nr_status'],viral_taxa=counts['viral'],artificial_taxa=counts['artificial'],cellular_taxa=counts['cellular'],other_known_taxa=counts['other_known'],unknown_taxa=counts['unknown'],viral_and_artificial_only=counts['viral']>0 and counts['artificial']>0 and counts['cellular']==counts['other_known']==counts['unknown']==0))
with open(d/'competitor-clades.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={'primary_classifications_unchanged':True,'scope':{},'interpretation':'Known hits outside the viral clade include artificial/engineered sequences, not only cellular organisms. This descriptive annotation uses the same eligible near-best taxa as the primary summary. Viral plus artificial-only evidence remains in the original mixed category; it is not promoted to a validated call. Unknown or unclassified evidence is not known-truth negative evidence.'}
for scope in ['protein_candidates','native_controls']:
 selected=[x for x in rows if x['scope']==scope];summary['scope'][scope]={'pairs':len(selected),'viral_and_artificial_only':sum(x['viral_and_artificial_only'] for x in selected),'any_cellular_competitor':sum(x['cellular_taxa']>0 for x in selected),'any_other_known_outside_viral_taxonomy':sum(x['other_known_taxa']>0 for x in selected),'any_unknown_taxonomy':sum(x['unknown_taxa']>0 for x in selected)}
(d/'competitor-clades.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
