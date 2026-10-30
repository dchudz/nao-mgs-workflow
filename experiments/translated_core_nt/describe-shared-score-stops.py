"""Measure score differences with and without stops; do not infer their cause."""
import argparse,collections,csv,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--right-score-column',default='fast_bitscore');a=p.parse_args()
matched=collections.defaultdict(collections.Counter);presence=collections.defaultdict(collections.Counter);total=0
with gzip.open(a.input,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  assert len(x['qseq'])==len(x['sseq'])
  n=sum(q==s=='*' for q,s in zip(x['qseq'],x['sseq']));delta=round(float(x[a.right_score_column])-float(x['tblastx_bitscore']),3)
  matched[n][delta]+=1;presence['contains_stop' if '*' in x['qseq']+x['sseq'] else 'no_stop_in_either_alignment'][delta]+=1;total+=1
def stats(v):return dict(HSPs=sum(v.values()),changed=sum(n for d,n in v.items() if d!=0),delta_range=[min(v),max(v)],most_common_bit_score_deltas=v.most_common(8))
result=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),right_score_column=a.right_score_column,shared_HSPs=total,alignment_stop_presence={k:stats(v) for k,v in presence.items()},matched_stop_count={str(k):stats(v) for k,v in sorted(matched.items())},limits=['Only already returned byte-identical HSP sequences, coordinates, frames and lengths are compared.','An alignment without a stop can still come from a query frame containing stops elsewhere; this measurement does not establish the scoring mechanism.','Removing stop-containing HSPs alone would leave the observed score differences in the stop-free subset. No filter was changed.','Score differences and threshold crossings are not established biological sensitivity gains. No new search was run.'])
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'shared_HSPs':total,'alignment_stop_presence':result['alignment_stop_presence']},indent=2))
