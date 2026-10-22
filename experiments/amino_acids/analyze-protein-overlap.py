"""Separate new original reads from changed complete-pair sequences; optionally join raw-gate attrition."""
import argparse,collections,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('pairs');p.add_argument('output');p.add_argument('--gate');a=p.parse_args()
gate={x['query_pair']:x for x in csv.DictReader(open(a.gate),delimiter='\t')} if a.gate else None
groups=collections.defaultdict(collections.Counter);seen=set()
for x in csv.DictReader(open(a.pairs),delimiter='\t'):
 q=x['query_pair'];assert q not in seen;seen.add(q)
 shared=x.get('exact_pair_also_native',x.get('exact_trimmed_pair_also_native'))=='True'
 original=int(x['original_pairs']);added=int(x.get('added_original_pairs',x.get('pairs_added_to_native')));assert 0<=added<=original
 category='exact_pair_also_native' if shared else 'new_exact_pair_has_native_original_member' if added<original else 'new_exact_pair_all_original_members_absent_from_native'
 for key in ['all',x['nr_status']]:
  c=groups[key];c['candidate_exact_pairs']+=1;c[category]+=1;c['original_members']+=original;c['additional_original_members']+=added
  if gate is not None:
   g=gate[q];assert int(g['original_members'])==original
   c['repeat_gate/'+g['gate_status']]+=1;c['repeat_gate/'+category+'/'+g['gate_status']]+=1
if gate is not None:assert seen==set(gate),'Mismatched candidate scopes'
summary={'groups':{k:dict(v) for k,v in groups.items()},'repeat_gate_included':gate is not None,'interpretation':['A new exact pair with an original member already assigned by the native path can reflect trimming or correction differences. It is not an additional original read.','Repeat-gate retention describes raw support from existing members; it does not rerun downstream processing or establish final pipeline parity.']}
Path(a.output).write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
