"""Verify report counts against independent saved tables and cluster assignments."""
import collections,csv,gzip,hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parent/'anchor-analysis';d=r/'result-v1'
def rows(p):
 with gzip.open(p,'rt') if str(p).endswith('.gz') else p.open() as f:yield from csv.DictReader(f,delimiter='\t')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
s=json.loads((d/'summary.json').read_text());e={x['query_pair']:x for x in rows(d/'pair-evidence.tsv')};ind=json.loads((r/'independent-saved-table-counts.json').read_text());calls=list(rows(d/'species-calls.tsv'))
extra={q:x for q,x in e.items() if x['full_exact_pair_also_native']=='False'}
assert len(extra)==ind['protein_only']==10620
assert dict(collections.Counter(x['raw_status'] for x in extra.values()))==ind['protein_only_status']
assert sum(x['exclude_constructs_status']=='viral_only' for x in extra.values())==876
assert sum(x['raw_status']=='viral_only' for x in e.values())==2615
checked=0
for m in s['clusters']:
 arm,depth,policy=m['arm'],m['depth'],m['policy'];folder=d/f'{arm}-{depth}-{policy}'
 assignments={x['seq_id']:x['sim_dup_exemplar'] for x in rows(folder/'output.tsv.gz')};supports=collections.defaultdict(set);additional=collections.defaultdict(set)
 for x in rows(folder/'protein-clusters.tsv.gz'):
  assert assignments[x['canonical_pair']]==x['similarity_exemplar']
  if x['eligible_both_concordant']=='True':
   sp=e[x['query_pair']][policy+'_species'];supports[sp].add(x['similarity_exemplar'])
   if x['native_absent']=='True':additional[sp].add(x['similarity_exemplar'])
 for x in calls:
  if (x['arm'],x['depth'],x['policy'])!=(arm,depth,policy):continue
  n=int(x['native_clusters']);p=len(supports[x['species']]);assert p==int(x['protein_clusters']);assert len(additional[x['species']])==int(x['protein_only_clusters'])
  assert (x['protein_call']=='True')==(p>=3);assert (x['gained_protein_call']=='True')==(p>=3 and n<3);checked+=1
inputs={}
for rp in ['nr-restore-receipt.json','full-restore-receipt.json']:
 for name,v in json.loads((r/rp).read_text())['verified_files'].items():
  assert sha(r/name)==v['sha256'];inputs[name]=v['sha256']
inputs['reference/total-virus-db-annotated.tsv.gz']=sha(r/'reference/total-virus-db-annotated.tsv.gz')
result={'verified_species_rows':checked,'independent_candidate_counts_match':True,'all_used_restored_files_reverified':inputs,'result_files':{str(p.relative_to(d)):sha(p) for p in d.rglob('*') if p.is_file()},'source_sha256':sha(Path(__file__))}
(r/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'verified_species_rows':checked,'files':len(result['result_files'])}))
