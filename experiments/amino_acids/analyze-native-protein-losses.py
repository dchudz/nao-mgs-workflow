"""Locate original native-read losses at measured protein preparation and alignment boundaries."""
import argparse,collections,csv,gzip,hashlib,json,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
p=argparse.ArgumentParser();p.add_argument('--nr-controls');p.add_argument('--raw-screen');a=p.parse_args();d=r/'cohort-storage/protein200';out=d/'native-losses';out.mkdir(exist_ok=True)
def table(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:yield from csv.DictReader(f,delimiter='\t')
order=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200]
batches={s:m.parent for m in (r/'cohort-storage/batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']}
assert set(order)<=batches.keys()
native={};aa12=set()
for sample in order:
 for arm in ['nt','aa12']:
  paths=list((batches[sample]/'output/marked').glob(f'*_{sample}_{arm}_duplicate_reads.tsv.gz'));assert len(paths)==1
  for x in table(paths[0]):
   if arm=='nt':assert x['seq_id'] not in native;native[x['seq_id']]=dict(sample=sample,taxid=x['aligner_taxid_lca'])
   else:aa12.add(x['seq_id'])
viral={x['query_pair'] for x in table(d/'viral-pairs.tsv')};prepared={};candidate_ids=set()
for x in table(d/'all.members.tsv.gz'):
 q=x['seq_id']
 if q in native:
  assert q not in prepared;prepared[q]=x['query_pair']
  if x['query_pair'] in viral:candidate_ids.add(q)
lost=native.keys()-candidate_ids;rows={q:dict(seq_id=q,**native[q],loss_boundary='prepared_without_eligible_viral_protein_alignment' if q in prepared else 'absent_from_prepared_protein_queries',protein_query_pair=prepared.get(q,''),native_AA12_assignment_preserved=q in aa12) for q in sorted(lost)}
controls={(x['sample'],x['exact_pair_sha256']):x['query_pair'] for x in table(r/'cohort-storage/nr-native-controls/selected.tsv')}
seen=set()
for sample in order:
 for key,raw in pairs(batches[sample]/'output/fastp'/f'{sample}_nt_fastp.fastq.gz'):
  q=key[0].decode().lstrip('@')
  if q not in rows:continue
  assert q not in seen;seen.add(q);z=rows[q];z['native_mate1_nt']=len(key[1]);z['native_mate2_nt']=len(key[3]);z['native_exact_pair_sha256']=hashlib.sha256(key[1]+b'\0'+key[3]).hexdigest();z['nr_control_query_pair']=controls.get((sample,z['native_exact_pair_sha256']),'')
assert seen==lost
control_evidence={x['query_pair']:x for x in table(a.nr_controls)} if a.nr_controls else {}
for z in rows.values():
 q=z['nr_control_query_pair'];z['nr_control_status']=control_evidence[q]['nr_status'] if q and a.nr_controls else 'pending' if q else 'not_selected';z['nr_control_taxid_lca']=control_evidence[q]['nr_taxid_lca'] if q and a.nr_controls else ''
raw_screen={x['seq_id']:x for x in table(a.raw_screen)} if a.raw_screen else {}
if a.raw_screen:
 assert raw_screen.keys()==rows.keys()
 for q,z in rows.items():
  assert raw_screen[q]['loss_boundary']==z['loss_boundary']
  value=raw_screen[q]['present_in_original_AA12_screen'];assert value in ('True','False')
  z['present_in_original_AA12_screen']=value=='True'
  if z['loss_boundary']=='prepared_without_eligible_viral_protein_alignment':assert z['present_in_original_AA12_screen']
with (out/'lost-original-pairs.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(next(iter(rows.values()))) if rows else ['seq_id'],delimiter='\t');w.writeheader();w.writerows(rows.values())
summary=dict(samples=len(order),native_original_pairs=len(native),native_original_pairs_present_in_protein_candidates=len(candidate_ids),lost_original_pairs=len(rows),distinct_native_exact_pairs_among_lost_original_members=len({(x['sample'],x['native_exact_pair_sha256']) for x in rows.values()}),loss_boundaries={},nr_controls_complete=bool(a.nr_controls),interpretation='Every loss is an original NT24 read pair absent from all viral-protein candidates. Preparation uses AA12-screen FASTP output followed by unchanged human/other contaminant filters. Absence from prepared queries does not by itself identify which upstream operation removed it. Prepared queries without an eligible viral-protein call passed that boundary but lack an alignment satisfying the frozen protein criteria. Native AA12 assignment preservation is an additional observed comparison, not a presumed causal explanation. Linked NR evidence belongs to the explicitly identified exact-sequence control query. Control-status counts use unique control queries, not extrapolated verdicts for original reads. Other losses remain unvalidated.')
for boundary in sorted({x['loss_boundary'] for x in rows.values()}):
 selected=[x for x in rows.values() if x['loss_boundary']==boundary];summary['loss_boundaries'][boundary]=dict(original_pairs=len(selected),native_AA12_assignment_preserved=sum(x['native_AA12_assignment_preserved'] for x in selected),nr_control_exact_pairs=len({x['nr_control_query_pair'] for x in selected if x['nr_control_query_pair']}),nr_control_status_exact_pairs=dict(collections.Counter({x['nr_control_query_pair']:x['nr_control_status'] for x in selected if x['nr_control_query_pair']}.values())))
summary['raw_screen_membership_measured']=bool(a.raw_screen)
if a.raw_screen:
 for boundary,z in summary['loss_boundaries'].items():
  selected=[x for x in rows.values() if x['loss_boundary']==boundary]
  z['present_in_original_AA12_screen']=sum(x['present_in_original_AA12_screen'] for x in selected)
  z['absent_from_original_AA12_screen']=sum(not x['present_in_original_AA12_screen'] for x in selected)
assert summary['lost_original_pairs']==json.loads((d/'native-comparison/summary.json').read_text())['tiers']['viral_reference_candidate']['lost_original_pairs']
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
