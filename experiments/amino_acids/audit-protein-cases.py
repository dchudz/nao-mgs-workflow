"""Select fixed-hash illustrative cases, retaining eligible evidence from each competing clade."""
import argparse,collections,csv,gzip,hashlib,heapq,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--evidence',required=True);p.add_argument('--protein',required=True);p.add_argument('--nr',required=True);p.add_argument('--taxonomy',required=True);p.add_argument('--output',required=True);a=p.parse_args()
protein=Path(a.protein);out=Path(a.output);out.mkdir(parents=True,exist_ok=True);groups=collections.defaultdict(list)
for x in csv.DictReader(open(a.evidence),delimiter='\t'):
 shared=x['exact_pair_also_native']=='True';status=x['nr_status']
 if not shared and status=='viral_only':group='added_viral_only_reference_identity_below80' if float(x['maximum_identity_pct'])<80 else 'added_viral_only_reference_identity_at_least80'
 elif not shared and status in ['mixed_viral_nonviral','nonviral_only','no_alignment']:group='added_'+status
 elif shared and status in ['mixed_viral_nonviral','nonviral_only']:group='also_native_'+status
 else:continue
 groups[group].append(x)
selected={};selection=[]
for group,rows in sorted(groups.items()):
 for x in sorted(rows,key=lambda x:hashlib.sha256(('NR-case-audit:'+x['query_pair']).encode()).digest())[:3]:
  x=dict(x);x['audit_stratum']=group;x['stratum_total']=len(rows);selected[x['query_pair']]=x;selection.append(x)
for x in selected.values():x['example_original_read_ids']=[]
with gzip.open(protein/'all.members.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['query_pair'] in selected:
   z=selected[x['query_pair']]['example_original_read_ids'];z.append(x['seq_id']);z[:]=heapq.nsmallest(3,z)
with gzip.open(protein/'all.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['query_pair'] in selected:
   z=selected[x['query_pair']]
   for key in ['query_seq','query_seq_rev']:
    seq=x[key];z[key+'_length']=len(seq);z[key+'_max_homopolymer']=max((len(m.group()) for m in re.finditer(r'(A+|C+|G+|T+)',seq.upper())),default=0)
parents={}
for line in open(Path(a.taxonomy)/'nodes.dmp'):
 x=line.split('|');parents[int(x[0])]=int(x[1])
cache={}
def category(t):
 if t in cache:return cache[t]
 original=t;seen=set();kind='other_known' if t in parents else 'unknown'
 while t in parents and t not in seen:
  if t==10239:kind='viral';break
  if t==81077:kind='artificial';break
  if t in [131567,2,2157,2759]:kind='cellular';break
  seen.add(t);t=parents[t]
 cache[original]=kind;return kind
hits=collections.defaultdict(list);best=collections.defaultdict(float)
for row in csv.reader(open(a.nr),delimiter='\t'):
 q,subject,score,evalue,length,identity,taxids=row;pair=q.rsplit('/',1)[0]
 if pair not in selected or int(length)<25:continue
 score=float(score);best[q]=max(best[q],score);taxa={int(t) for t in taxids.split(';') if t.isdigit()} or {0}
 for kind in {category(t) for t in taxa}:
  z=hits[q,kind];z.append((score,subject,evalue,int(length),float(identity),taxids));z.sort(key=lambda x:(-x[0],x[1]));del z[2:]
rows=[];wanted=set()
for (q,kind),values in sorted(hits.items()):
 for score,subject,evalue,length,identity,taxids in values:
  if score<.9*best[q]:continue
  taxa={int(t) for t in taxids.split(';') if t.isdigit()} or {0};wanted.update(taxa)
  rows.append(dict(query=q,competitor_clade=kind,subject=subject,bitscore=score,evalue=evalue,alignment_aa=length,identity_pct=identity,taxids=taxids))
names={}
for line in open(Path(a.taxonomy)/'names.dmp'):
 x=[f.strip() for f in line.split('|')]
 if x[3]=='scientific name' and int(x[0]) in wanted:names[int(x[0])]=x[1]
for x in rows:x['taxon_names']='; '.join(names.get(int(t),'unknown') for t in x['taxids'].split(';') if t.isdigit())
for x in selection:x['example_original_read_ids']=';'.join(x['example_original_read_ids'])
for name,data in [('cases.tsv',selection),('case-alignments.tsv',rows)]:
 with open(out/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'cases':len(selection),'strata':{k:len(v) for k,v in groups.items()},'selection':'Up to three SHA256-ranked cases per descriptive outcome stratum. Representatives are illustrations, not an accuracy sample.','evidence':'Up to two highest-scoring eligible alignments per mate and competitor clade; each retained row reaches90% of the best eligible mate score and25alignedaminoacids. Raw full search output is retained separately. Identity threshold refers to the viral-reference match, not NR.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
