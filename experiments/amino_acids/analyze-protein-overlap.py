"""Separate new original reads from changed complete-pair sequences; optionally join raw-gate attrition."""
import argparse,collections,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('pairs');p.add_argument('output');p.add_argument('--gate');p.add_argument('--screen-membership');p.add_argument('--native-arm',choices=['nt','nt20']);a=p.parse_args()
gate={x['query_pair']:x for x in csv.DictReader(open(a.gate),delimiter='\t')} if a.gate else None
screen=None
if a.screen_membership:
 assert a.native_arm,'Specify the native arm for raw-screen membership'
 screen_rows=list(csv.DictReader(open(a.screen_membership),delimiter='\t'));screen={x['query_pair']:x for x in screen_rows};assert len(screen)==len(screen_rows)
groups=collections.defaultdict(collections.Counter);seen=set()
for x in csv.DictReader(open(a.pairs),delimiter='\t'):
 q=x['query_pair'];assert q not in seen;seen.add(q)
 shared=x.get('exact_pair_also_native',x.get('exact_trimmed_pair_also_native'))=='True'
 original=int(x['original_pairs']);added=int(x.get('added_original_pairs',x.get('pairs_added_to_native')));assert 0<=added<=original
 category='exact_pair_also_native' if shared else 'new_exact_pair_has_native_original_member' if added<original else 'new_exact_pair_all_original_members_absent_from_native'
 for key in ['all',x['nr_status']]:
  c=groups[key];c['candidate_exact_pairs']+=1;c[category]+=1;c['original_members']+=original;c['additional_original_members']+=added
  if screen is not None:
   raw=screen[q];assert int(raw['original_members'])==original
   retained=int(raw[a.native_arm+'_retained_original_members']);assert 0<=retained<=original
   status=raw[a.native_arm+'_screen_status'];assert status==('none_retained' if retained==0 else 'all_retained' if retained==original else 'some_retained')
   c['initial_screen/'+category+'/'+status]+=1
   c['initial_screen_original_members/'+category+'/retained']+=retained
   c['initial_screen_original_members/'+category+'/excluded']+=original-retained
  if gate is not None:
   g=gate[q];assert int(g['original_members'])==original
   c['repeat_gate/'+g['gate_status']]+=1;c['repeat_gate/'+category+'/'+g['gate_status']]+=1
if screen is not None:assert seen==set(screen),'Mismatched raw-screen candidate scope'
if gate is not None:assert seen==set(gate),'Mismatched candidate scopes'
summary={'groups':{k:dict(v) for k,v in groups.items()},'repeat_gate_included':gate is not None,'initial_screen_included':screen is not None,'initial_screen_native_arm':a.native_arm,'interpretation':['A new exact pair with an original member already assigned by the native path can reflect trimming or correction differences. It is not an additional original read.','For queries whose original members are all absent from final native assignments, raw retained members passed the initial nucleotide screen but were absent later; raw excluded members never appeared in that screen output. This does not identify the responsible later operation.','Repeat-gate retention describes raw support from existing members; it does not rerun downstream processing or establish final pipeline parity.']}
Path(a.output).write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
