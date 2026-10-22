"""Compare protein evidence to native calls using original IDs and the same complete trimmed pairs."""
import argparse,collections,csv,gzip,json,sys
from pathlib import Path
from union_pairs import pairs
p=argparse.ArgumentParser();p.add_argument('--metadata',required=True);p.add_argument('--members',required=True);p.add_argument('--viral',required=True);p.add_argument('--nr');p.add_argument('--native',required=True);p.add_argument('--fastp',required=True);p.add_argument('--prefix',default='');p.add_argument('--output',required=True);a=p.parse_args()
def table(path):
 with (gzip.open(path,'rt') if path.endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))
native={x['seq_id'] for x in table(a.native)};native_keys=set();seen=set()
for key,raw in pairs(a.fastp):
 identifier=key[0].decode().lstrip('@')
 if identifier in native:native_keys.add((key[1].decode(),key[3].decode()));seen.add(identifier)
assert seen==native,('Native IDs absent from original FASTP',len(native-seen))
meta={x['query_pair']:x for x in table(a.metadata) if x['query_pair'].startswith(a.prefix)}
members=collections.defaultdict(set)
for x in table(a.members):
 if x['query_pair'] in meta:members[x['query_pair']].add(x['seq_id'])
viral={x['query_pair']:x for x in table(a.viral) if x['query_pair'] in meta};nr={x['query_pair']:x for x in table(a.nr)} if a.nr else {}
if a.nr:assert viral.keys()<=nr.keys(),'Incomplete broad-reference validation'
rows=[]
for q,call in viral.items():
 m=meta[q];ids=members[q];assert len(ids)==int(m['occurrences']),(q,len(ids),m['occurrences'])
 key=(m['query_seq'],m['query_seq_rev']);n=nr.get(q,{})
 rows.append(dict(query_pair=q,original_pairs=len(ids),pairs_also_native=len(ids&native),pairs_added_to_native=len(ids-native),exact_trimmed_pair_also_native=key in native_keys,protein_taxid_lca=call['protein_taxid_lca'],nr_status=n.get('nr_status','pending'),taxonomic_relation=n.get('taxonomic_relation','pending')))
with open(a.output,'w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(rows)
candidate_ids=set().union(*(members[q] for q in viral)) if viral else set()
candidate_keys={(meta[q]['query_seq'],meta[q]['query_seq_rev']) for q in viral}
summary={'native_pairs':len(native),'native_exact_complete_trimmed_pairs':len(native_keys),'protein_viral_reference_candidates':len(rows),'protein_candidate_original_pairs':sum(x['original_pairs'] for x in rows),'candidate_original_pairs_added_to_native':sum(x['pairs_added_to_native'] for x in rows),'native_pairs_absent_from_protein_candidates':len(native-candidate_ids),'native_exact_trimmed_pairs_absent_from_protein_candidates':len(native_keys-candidate_keys),'candidate_exact_trimmed_pairs_absent_from_native':sum(not x['exact_trimmed_pair_also_native'] for x in rows),'broad_reference_status':dict(collections.Counter(x['nr_status'] for x in rows)),'viral_only_NR_exact_pairs_absent_from_native':None if not a.nr else sum(x['nr_status']=='viral_only' and not x['exact_trimmed_pair_also_native'] for x in rows),'interpretation':'Candidates are not validated detections. Common exact-pair comparison uses complete original FASTP pairs on both sides; native TSVs may omit an unaligned mate.'}
Path(a.output+'.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
