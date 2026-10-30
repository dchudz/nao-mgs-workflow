"""Measure query-threading work and exact alignment parity with the frozen baseline."""
import argparse,collections,gzip,json,re
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('baseline');a=p.parse_args()
d=Path(a.root);baseline=Path(a.baseline)
rows=json.loads((d/'results.json').read_text())
assert len(rows)==3 and [x['mode'] for x in rows]==[1,0,1]

def hits(path):
 with gzip.open(path,'rt') as f:return collections.Counter(line.rstrip('\n') for line in f)

def difference(left,right):
 only_left=left-right;only_right=right-left
 queries=sorted({line.split('\t',1)[0] for line in only_left.keys()|only_right.keys()})
 return {'alignment_multiset_equal':left==right,'left_rows':sum(left.values()),'right_rows':sum(right.values()),'rows_only_left':sum(only_left.values()),'rows_only_right':sum(only_right.values()),'differing_query_count':len(queries),'differing_query_ids':queries}

reference=hits(baseline/'hits.tsv.gz');outputs=[]
for x in rows:
 run=d/f"run{x['run']}-mode{x['mode']}";timing=(run/'blast.time').read_text()
 def number(label):
  m=re.search(re.escape(label)+r':\s*([\d.]+)',timing);assert m,label;return float(m[1])
 assert number('Exit status')==0
 x['cpu_seconds']=number('User time (seconds)')+number('System time (seconds)')
 x['peak_rss_kib']=int(number('Maximum resident set size (kbytes)'))
 x['log_text']=(run/'run.log').read_text()
 observed=hits(run/'hits.tsv.gz');outputs.append(observed)
 x['baseline_parity']=difference(reference,observed)
 changes={side:[{'multiplicity':count,'alignment_columns':line.split('\t')} for line,count in sorted(delta.items())] for side,delta in [('baseline_only',reference-observed),('trial_only',observed-reference)]}
 (run/'alignment-differences.json').write_text(json.dumps(changes,indent=2))
summary={'runs':rows,'baseline':str(baseline),'mode1_repeat_parity':difference(outputs[0],outputs[2]),'notes':['Exact output multisets include every reported column and duplicate alignment multiplicity.','Shared warm-cache worker; timings are not isolated service latency.','All new trials have an explicit 400 GiB memory ceiling; the preceding cap20 baseline did not.','No assignment parity claim is inferred from alignment parity failure.']}
(d/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps([{k:v for k,v in x.items() if k not in ('command','log_text','baseline_parity')}|{'baseline_parity':{k:v for k,v in x['baseline_parity'].items() if k!='differing_query_ids'}} for x in rows]))
print(json.dumps(summary['mode1_repeat_parity']))
