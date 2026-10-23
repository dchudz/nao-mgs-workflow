"""Describe all strict pair transitions and reconstruct every newly viral-only pair."""
import argparse,collections,csv,hashlib,json
from pathlib import Path
csv.field_size_limit(100000000)
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();d=a.directory;out=d/'transition-analysis';out.mkdir(exist_ok=False)
paths=[d/'pair-analysis/per-pair.tsv',d/'analysis/codon-audit-tblastx-word4-vs-blastn.json',d/'analysis/taxonomy-description/taxa.tsv']
rows=list(csv.DictReader(paths[0].open(),delimiter='\t'));calls={(x['pair'],x['method']):x for x in rows if x['threshold']=='strict_E1e-5' or x['method']=='protein_NR'}
assert len(calls)==232*4
names={x['taxid']:x for x in csv.DictReader(paths[2].open(),delimiter='\t')}
audit=json.loads(paths[1].read_text());changes=[];gains=[];losses=[];counts=collections.Counter()
for pair in sorted({x['pair'] for x in rows}):
 left=calls[pair,'blastn'];right=calls[pair,'tblastx-word4'];nr=calls[pair,'protein_NR']
 counts[left['stratum']+'/'+left['natural_status']+' -> '+right['natural_status']]+=1
 if left['natural_status']!=right['natural_status'] or (left['status']=='viral_only')!=(right['status']=='viral_only'):
  changes.append(dict(pair=pair,stratum=left['stratum'],left_status=left['status'],right_status=right['status'],left_natural_status=left['natural_status'],right_natural_status=right['natural_status'],translated_aligned_mates=int(right['aligned_mates']),translated_artificial_taxa=int(right['artificial_taxa']),translated_cellular_taxa=int(right['cellular_taxa']),translated_unknown_taxa=int(right['unknown_taxa'])))
 if left['status']=='viral_only' and right['status']!='viral_only':losses.append(changes[-1])
 if left['natural_status']=='no_alignment' and right['natural_status']=='viral_only':
  hits=[x for x in audit['rows'] if x['query'].rsplit('/',1)[0]==pair]
  taxon=names.get(right['known_natural_lca'],{})
  gains.append(dict(pair=pair,sample=left['sample'],stratum=left['stratum'],translated_aligned_mates=int(right['aligned_mates']),translated_natural_lca=right['known_natural_lca'],translated_name=taxon.get('name',''),translated_rank=taxon.get('rank',''),NR_status=nr['status'],NR_LCA=nr['taxid_lca'],reconstructed_mates=hits))
result=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},threshold='strict_E1e-5',left='blastn',right='tblastx-word4',pair_count=232,natural_status_transitions=dict(counts),changes=changes,newly_viral_only_pairs=gains,raw_viral_only_losses=losses,limits=['All gains and losses are included; neither higher positive counts nor narrower assignments establish accuracy.','The native core_nt search settings and per-read classifications remain unchanged.','Pair evidence combines only the two mates of that pair; a single aligned mate is explicitly identified.','This selected one-pair-per-sample panel cannot establish new three-exemplar sample/species calls.','NR is preserved contextual evidence from a different reference and filtering scheme, with the separately measured batching limitation.'])
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'newly_viral_only_pairs':len(gains),'raw_viral_only_losses':len(losses),'raw_losses_remaining_naturally_viral':sum(x['right_natural_status']=='viral_only' for x in losses),'gains':[{'pair':x['pair'],'name':x['translated_name'],'rank':x['translated_rank'],'aligned_mates':x['translated_aligned_mates']} for x in gains]},indent=2))
