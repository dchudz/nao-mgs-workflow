"""Attribute existing prepared-read losses to reported alignment availability and fixed length cutoff."""
import argparse,collections,csv,json,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--losses',required=True);p.add_argument('--hits',required=True);p.add_argument('--output',required=True);a=p.parse_args()
selected=[x for x in csv.DictReader(open(a.losses),delimiter='\t') if x['loss_boundary']=='prepared_without_eligible_viral_protein_alignment']
members=collections.defaultdict(list)
for x in selected:assert x['protein_query_pair'];members[x['protein_query_pair']].append(x['seq_id'])
wanted={q.encode():q for q in members};counts=collections.Counter();mates=collections.defaultdict(set);best={};maximum=collections.defaultdict(int);scanned=0;start=time.monotonic()
with open(a.hits,'rb') as f:
 for line in f:
  scanned+=1;q=line.split(b'\t',1)[0];pair=q.rsplit(b'/',1)[0]
  if pair not in wanted:continue
  q=q.decode();name=wanted[pair];mate=q.rsplit('/',1)[1];assert mate in ('1','2')
  fields=line.decode().rstrip('\n').split('\t');assert len(fields)==12
  length=int(fields[4]);assert length<25,(name,'Eligible alignment found for a supposedly absent candidate')
  counts[name]+=1;mates[name].add(mate);maximum[name]=max(maximum[name],length)
  row=dict(query=q,subject=fields[1],bitscore=float(fields[2]),evalue=float(fields[3]),alignment_aa=length,identity_pct=float(fields[5]),query_nt=int(fields[6]),query_start=int(fields[7]),query_end=int(fields[8]),query_frame=int(fields[11]))
  if q not in best or row['bitscore']>best[q]['bitscore']:best[q]=row
rows=[]
for q,ids in sorted(members.items()):
 rows.append(dict(query_pair=q,original_lost_pairs=len(ids),original_read_IDs=';'.join(sorted(ids)),observed_boundary='reported_alignment_below_25_aa_only' if counts[q] else 'no_reported_viral_alignment',reported_alignment_rows=counts[q],aligned_mates=len(mates[q]),maximum_alignment_aa=maximum[q]))
out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
for name,data in [('query-loss-boundaries.tsv',rows),('best-subthreshold-alignments.tsv',list(best.values()))]:
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['query'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'selected_original_pairs':len(selected),'selected_prepared_queries':len(members),'raw_alignment_rows_scanned':scanned,'wall_seconds':time.monotonic()-start,'boundaries':{},'interpretation':'Existing DIAMOND output was inspected without rerunning search or changing primary classifications. No reported alignment means no row in the saved search output under its frozen parameters; it does not establish absence of homology. Subthreshold alignments had already passed search E-value/coverage reporting but were excluded by the frozen25-amino-acid summary rule. No NR verdict is extrapolated.'}
for boundary in sorted({x['observed_boundary'] for x in rows}):
 z=[x for x in rows if x['observed_boundary']==boundary];summary['boundaries'][boundary]={'prepared_queries':len(z),'original_pairs':sum(x['original_lost_pairs'] for x in z)}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
