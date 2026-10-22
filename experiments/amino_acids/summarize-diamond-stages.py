"""Extract recorded stage durations only from a successfully completed DIAMOND log."""
import argparse,collections,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('log');p.add_argument('output');a=p.parse_args();text=Path(a.log).read_text()
assert re.search(r'Exit status:\s*0\s*$',text) and 'Command terminated by signal' not in text,'A completed search is required'
counts=collections.Counter();seconds=collections.Counter();chunks=[]
for line in text.splitlines():
 m=re.match(r'^(.*?)\.\.\.\s+\[([0-9.]+)s\]',line)
 if m:counts[m[1]]+=1;seconds[m[1]]+=float(m[2])
 m=re.match(r'Processing query block (\d+), reference block (\d+)/(\d+), shape (\d+)/(\d+), index chunk (\d+)/(\d+)\.',line)
 if m:chunks.append(tuple(map(int,m.groups())))
assert chunks and len(chunks)==len(set(chunks)), 'Missing or duplicate chunk records'
dimensions={(q,B,S,C) for q,_,B,_,S,_,C in chunks}
expected={(q,b,B,s,S,c,C) for q,B,S,C in dimensions for b in range(1,B+1) for s in range(1,S+1) for c in range(1,C+1)}
assert set(chunks)==expected,'Incomplete query/reference/shape/chunk coverage'
for label in ['Building reference seed array','Building query seed array','Computing hash join','Searching alignments']:assert counts[label]==len(chunks),(label,counts[label],len(chunks))
m=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)',text);assert m;elapsed=0
for part in m[1].split(':'):elapsed=elapsed*60+float(part)
groups={'reference_loading':['Loading reference sequences'],'reference_masking':['Masking reference'],'reference_histograms_and_seed_arrays':['Building reference histograms','Building reference seed array'],'query_histograms_and_seed_arrays':['Building query histograms','Building query seed array'],'hash_join':['Computing hash join'],'alignment_search_and_computation':['Searching alignments','Computing alignments']}
summary=dict(wall_seconds=elapsed,logged_chunks=len(chunks),groups={name:{'logged_seconds':sum(seconds[label] for label in labels),'labels':labels} for name,labels in groups.items()},stages={label:{'calls':counts[label],'logged_seconds':seconds[label]} for label in counts},interpretation='DIAMOND-reported elapsed stage durations from the completed native-BLAST-format NR search. Histogram and seed-array construction recur during this search and do not create a persistent reusable DIAMOND database. These are reported stage durations, not CPU time or a separate benchmark. Total wall time also includes remaining setup, masking, cleanup and output work.')
Path(a.output).write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='stages'}))
