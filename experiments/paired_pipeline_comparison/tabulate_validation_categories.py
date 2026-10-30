"""Aggregate existing validation annotations by the symmetric reporting category.

No sequences are exported, no taxonomic assignments are changed, and no searches run.
"""
import argparse,collections,csv,gzip,hashlib,json,math
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def table(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:
  yield from csv.DictReader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
def aggregate(rows,taxonomy,categories):
 counts=collections.defaultdict(collections.Counter)
 for row in rows:
  sp=taxonomy[row['aligner_taxid_lca']]['taxid_species']
  category=categories.get(sp,'below_threshold_both') if sp not in ('','NA','0','None') else 'above_species_or_unresolved'
  c=counts[category];c['rows']+=1;c[row['validation_status']]+=1
  c['joint_exemplars']+=row['seq_id']==row['prim_align_dup_exemplar']==row['sim_dup_exemplar']
  try:distance=float(row['validation_distance_aligner'])
  except (ValueError,TypeError):continue
  if math.isfinite(distance):
   assert distance>=0;c['rows_with_taxonomic_distance']+=1;c['sum_taxonomic_distance']+=distance;c['agreement_rows']+=distance==0
 return counts

def main():
 p=argparse.ArgumentParser()
 for name in ('root','reference','audit','out'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 taxonomy={r['taxid']:r for r in table(a.reference)};reference_sha=sha(a.reference);completed=[]
 for proof in sorted((a.root/'comparisons').glob('*/completed-evidence-backup.json')):
  d=proof.parent;dest=a.out/d.name
  if (dest/'complete.json').exists():completed.append(d.name);continue
  audit=a.audit/d.name
  if not (audit/'complete.json').exists():continue
  meta=json.loads((audit/'complete.json').read_text());assert meta['reference_sha256']==reference_sha
  species=audit/'species.tsv';assert sha(species)==meta['outputs']['species.tsv']
  categories={r['species_taxid']:r['threshold_3'] for r in table(species)}
  expected=json.loads((d/'comparison-report/validation-status-counts.json').read_text());result=[];inputs={str(species):sha(species)}
  for arm in ('production','protein_variant'):
   paths=list((d/(arm+'-downstream/output/results_downstream')).glob('*_validation_hits.tsv.gz'));assert len(paths)==1
   path=paths[0];assert sha(path)==meta['input_tables'][arm]['sha256'];inputs[str(path)]=sha(path)
   counts=aggregate(table(path),taxonomy,categories)
   assert sum(x['rows'] for x in counts.values())==expected[arm]['rows']
   assert sum(x['joint_exemplars'] for x in counts.values())==expected[arm]['joint_exemplars']
   for status,n in expected[arm]['validation_status'].items():assert sum(x[status] for x in counts.values())==n
   for category,c in sorted(counts.items()):
    attempted=c['aligned']+c['no_alignment'];n=c['rows_with_taxonomic_distance']
    result.append(dict(group=d.name,arm=arm,category=category,rows=c['rows'],joint_exemplars=c['joint_exemplars'],attempted=attempted,aligned=c['aligned'],no_alignment=c['no_alignment'],not_sampled=c['not_sampled'],other_status_rows=c['rows']-attempted-c['not_sampled'],no_alignment_rate=c['no_alignment']/attempted if attempted else '',rows_with_taxonomic_distance=n,agreement_rows=c['agreement_rows'],agreement_rate=c['agreement_rows']/n if n else '',sum_taxonomic_distance=c['sum_taxonomic_distance'],mean_taxonomic_distance=c['sum_taxonomic_distance']/n if n else ''))
  dest.mkdir(exist_ok=True);output=dest/'categories.tsv'
  with output.open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
  (dest/'complete.json').write_text(json.dumps(dict(label=d.name,threshold=3,analysis_sha256=sha(Path(__file__)),reference_sha256=reference_sha,input_hashes=inputs,output_sha256=sha(output),scope='All existing assigned rows grouped by the species reporting category; no query export or new searches.'),indent=2)+'\n')
  completed.append(d.name);print(json.dumps({'label':d.name,'category_rows':len(result)}),flush=True)
 print(json.dumps({'completed':len(completed)}),flush=True)
if __name__=='__main__':main()
