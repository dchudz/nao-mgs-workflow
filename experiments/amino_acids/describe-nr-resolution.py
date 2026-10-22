"""Describe NR taxonomic resolution without promoting broad viral LCAs to species support."""
import argparse,collections,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('evidence');p.add_argument('taxonomy');p.add_argument('output');a=p.parse_args();tax=Path(a.taxonomy)
parents={};ranks={}
for line in open(tax/'nodes.dmp'):
 f=[x.strip() for x in line.split('|')];parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2]
def species(t):
 seen=set()
 while t in parents and t not in seen:
  if ranks[t]=='species':return t
  seen.add(t);t=parents[t]
 return None
rows=[];counts=collections.defaultdict(collections.Counter);wanted=set()
for x in csv.DictReader(open(a.evidence),delimiter='\t'):
 if x['nr_status']!='viral_only':continue
 v=int(x['viral_reference_taxid']);n=int(x['nr_taxid_lca']);vs=species(v);ns=species(n)
 relation='same_species' if vs is not None and vs==ns else 'different_species' if vs is not None and ns is not None else 'NR_above_species' if ns is None and n in parents else 'unresolved'
 key='also_native' if x['exact_pair_also_native']=='True' else 'absent_from_native'
 counts[key].update({'pairs':1,'rank/'+ranks.get(n,'unknown'):1,'species_relation/'+relation:1})
 rows.append(dict(query_pair=x['query_pair'],sample=x['sample'],exact_pair_also_native=x['exact_pair_also_native'],viral_reference_taxid=v,nr_taxid_lca=n,nr_rank=ranks.get(n,'unknown'),reference_species=vs,nr_species=ns,species_relation=relation))
 wanted.update(t for t in [v,n,vs,ns] if t is not None)
names={}
for line in open(tax/'names.dmp'):
 f=[x.strip() for x in line.split('|')]
 if f[3]=='scientific name' and int(f[0]) in wanted:names[int(f[0])]=f[1]
for x in rows:
 for key in ['viral_reference_taxid','nr_taxid_lca','reference_species','nr_species']:x[key+'_name']=names.get(x[key],'')
out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
with open(out/'resolution.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'groups':{k:dict(v) for k,v in counts.items()},'interpretation':'NR viral-only evidence can have an LCA above species. Same-species means the NR LCA and viral-reference LCA resolve to the same species in the NR snapshot taxonomy; it is still sequence evidence rather than known-truth validation. No primary evidence tier is changed.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
