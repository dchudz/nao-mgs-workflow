"""Pool category-level validation counts and reconcile to the existing arm audit."""
import argparse,collections,csv,hashlib,json
from pathlib import Path

def rows(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
p=argparse.ArgumentParser()
for key in ('study','input','tables'):p.add_argument('--'+key,type=Path,required=True)
a=p.parse_args();cases=[r for r in rows(a.study/'results/aggregate/runs.tsv') if r['status']=='complete'];pooled=collections.defaultdict(collections.Counter);input_hashes={}
fields=('rows','joint_exemplars','attempted','aligned','no_alignment','not_sampled','other_status_rows','rows_with_taxonomic_distance','agreement_rows','sum_taxonomic_distance')
for case in cases:
 d=a.input/case['label'];path=d/'categories.tsv';proof=json.loads((d/'complete.json').read_text());assert sha(path)==proof['output_sha256'];input_hashes[str(path)]=sha(path)
 for r in rows(path):
  key=(int(case['depth_million_pairs']),r['arm'],r['category'])
  for field in fields:pooled[key][field]+=float(r[field]) if field=='sum_taxonomic_distance' else int(r[field])
output=[]
for (depth,arm,category),c in sorted(pooled.items()):
 n=c['rows_with_taxonomic_distance'];attempted=c['attempted'];total=c['rows']
 output.append(dict(depth_million_pairs=depth,units=sum(int(r['depth_million_pairs'])==depth for r in cases),arm=arm,category=category,**c,attempted_fraction=attempted/total if total else '',no_alignment_rate=c['no_alignment']/attempted if attempted else '',agreement_rate=c['agreement_rows']/n if n else '',mean_taxonomic_distance=c['sum_taxonomic_distance']/n if n else ''))
for r in rows(a.tables/'validation-audit.tsv'):
 group=[x for x in output if x['depth_million_pairs']==int(r['depth_million_pairs']) and x['arm']==r['arm']]
 assert group and group[0]['units']==int(r['units'])
 for key in ('rows','attempted','aligned','no_alignment','not_sampled','other_status_rows','rows_with_taxonomic_distance','agreement_rows'):
  assert sum(x[key] for x in group)==int(r[key]),(key,r['arm'])
 n=int(r['rows_with_taxonomic_distance'])
 if n:assert abs(sum(x['sum_taxonomic_distance'] for x in group)/n-float(r['mean_taxonomic_distance']))<1e-9
path=a.tables/'validation-by-category.tsv'
with path.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t');w.writeheader();w.writerows(output)
(a.tables/'validation-category-provenance.json').write_text(json.dumps(dict(accepted_comparisons=len(cases),analysis_sha256=sha(Path(__file__)),input_hashes=input_hashes,pooled_audit_sha256=sha(a.tables/'validation-audit.tsv'),output_sha256=sha(path),reconciled_to_pooled_audit=True),indent=2)+'\n')
print(json.dumps({'accepted_comparisons':len(cases),'category_rows':len(output),'reconciled_to_pooled_audit':True}))
