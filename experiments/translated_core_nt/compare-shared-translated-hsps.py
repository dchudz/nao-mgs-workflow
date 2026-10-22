"""Separate retrieval changes from score changes on byte-identical translated alignments."""
import argparse,collections,csv,gzip,hashlib,json
from pathlib import Path
FIELDS='qseqid sseqid staxids bitscore evalue length pident qcovs qcovhsp qlen slen qstart qend sstart send qframe sframe qseq sseq'.split()
KEY='qseqid sseqid length qlen slen qstart qend sstart send qframe sframe qseq sseq'.split()
p=argparse.ArgumentParser();p.add_argument('pilot',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=False)
methods=['tblastx','query-tblastn-fast'];tables={};duplicates={}
for method in methods:
 values={};n=0
 for row in csv.DictReader((a.pilot/method/'alignments.tsv').open(),fieldnames=FIELDS,delimiter='\t'):
  key=tuple(row[k] for k in KEY);n+=1
  if key in values:assert all(values[key][k]==row[k] for k in ['bitscore','evalue','staxids'])
  values[key]=row
 tables[method]=values;duplicates[method]=n-len(values)
l,r=(tables[x] for x in methods);common=l.keys()&r.keys();queries=[x[1:].split()[0] for x in (a.pilot/'queries.fna').read_text().splitlines() if x.startswith('>')];assert len(queries)==43
counts={q:collections.Counter() for q in queries};cross=set();deltas=[];stop_cross=0
for method in methods:
 for key in tables[method]:counts[key[0]][method+'_unique_HSPs']+=1
rows=[]
for key in sorted(common):
 x,y=l[key],r[key];assert x['staxids']==y['staxids'];q=key[0];counts[q]['shared_identical_HSPs']+=1;delta=float(y['bitscore'])-float(x['bitscore']);deltas.append(delta);changed=delta!=0
 if changed:counts[q]['shared_HSPs_with_different_bitscore']+=1
 left_pass=float(x['evalue'])<=1e-5;right_pass=float(y['evalue'])<=1e-5
 if left_pass!=right_pass:
  counts[q]['same_HSP_crosses_strict_E_cutoff']+=1;cross.add(q);stop_cross+=int('*' in x['qseq']+x['sseq'])
 rows.append(dict(zip(KEY,key),tblastx_bitscore=x['bitscore'],fast_bitscore=y['bitscore'],tblastx_evalue=x['evalue'],fast_evalue=y['evalue'],tblastx_passes_strict=left_pass,fast_passes_strict=right_pass,has_stop='*' in x['qseq']+x['sseq']))
with gzip.open(a.output/'shared-identical-alignments.tsv.gz','wt') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
per_query=[dict(query=q,**counts[q]) for q in queries];fields=['query']+sorted({k for row in per_query for k in row if k!='query'})
with (a.output/'per-query.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(per_query)
result={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_hashes':{method:hashlib.sha256((a.pilot/method/'alignments.tsv').read_bytes()).hexdigest() for method in methods},'query_count':len(queries),'unique_HSPs':{m:len(tables[m]) for m in methods},'duplicate_identical_HSP_rows':duplicates,'shared_identical_HSPs':len(common),'only_tblastx_HSPs':len(l.keys()-r.keys()),'only_fast_HSPs':len(r.keys()-l.keys()),'shared_HSPs_with_different_bitscores':sum(x!=0 for x in deltas),'fast_minus_tblastx_bitscore_range':[min(deltas),max(deltas)],'same_HSP_crosses_strict_E_cutoff':sum(x['same_HSP_crosses_strict_E_cutoff'] for x in counts.values()),'cutoff_crossing_HSPs_with_stops':stop_cross,'queries_with_identical_HSP_cutoff_crossing':sorted(cross),'limits':['Matching requires identical query, subject, coordinates, frames, lengths and both aligned amino-acid strings.','A threshold crossing on the same HSP is a change in reported significance, not newly retrieved sequence evidence.','Unshared HSPs may reflect heuristic retrieval, output caps or frame pooling; this comparison does not identify one unique cause.','Scores and E-values are not assumed interchangeable between programs.']};(a.output/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
