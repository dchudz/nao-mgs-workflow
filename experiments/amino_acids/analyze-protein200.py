"""Compare frozen first-200 protein calls with the original NT24 path on common exact pair keys."""
import argparse,collections,csv,gzip,json,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'))
from union_pairs import pairs
p=argparse.ArgumentParser();p.add_argument('--native-arm',choices=['nt','nt22','nt20','aa12','union'],default='nt');p.add_argument('--protein',required=True);p.add_argument('--nr');p.add_argument('--samples-json',help='Explicit ordered sample IDs; default remains frozen first200');p.add_argument('--output',required=True);a=p.parse_args()
d=Path(a.protein);out=Path(a.output);out.mkdir(exist_ok=True)
def table(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:yield from csv.DictReader(f,delimiter='\t')
order=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200]
if a.samples_json:
 order=json.loads(Path(a.samples_json).read_text())
 assert order and all(isinstance(x,str) for x in order) and len(order)==len(set(order)), 'Invalid ordered sample IDs'
batches={s:m.parent for m in (r/'cohort-storage/batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']}
assert set(order)<=set(batches)
viral={x['query_pair']:x for x in table(d/'viral-pairs.tsv')}
assert {q.split('_',1)[0] for q in viral} <= set(order), 'Protein calls outside requested samples'
nr={x['query_pair']:x for x in table(a.nr) if x['query_pair'] in viral} if a.nr else {}
if a.nr:assert viral.keys()==nr.keys(), 'Incomplete broad reference validation'
meta=collections.defaultdict(dict)
for x in table(d/'all.tsv.gz'):
 q=x['query_pair']
 if q in viral:meta[q.split('_',1)[0]][q]=x
members=collections.defaultdict(set)
for x in table(d/'all.members.tsv.gz'):
 if x['query_pair'] in viral:members[x['query_pair']].add(x['seq_id'])
parents={};ranks={}
for line in open(r/'reference/results/taxonomy-nodes.dmp'):
 f=[x.strip() for x in line.split('|')];parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2]
def species(t):
 seen=set()
 while t in parents and t not in seen:
  if ranks[t]=='species':return t
  seen.add(t);t=parents[t]
 return None
with gzip.open(r/'reference/results/total-virus-db-annotated.tsv.gz','rt') as f:annotations={int(x['taxid']):x for x in csv.DictReader(f,delimiter='\t') if x['rank']=='species'}
def human(t):return annotations.get(t,{}).get('infection_status_human')=='1'
sample_rows=[];pair_rows=[];species_rows=[];species_changes=[]
for sample in order:
 batch=batches[sample];files=list((batch/'output/marked').glob(f'*_{sample}_{a.native_arm}_duplicate_reads.tsv.gz'));assert len(files)==1
 native={x['seq_id']:x for x in table(files[0])};native_keys=set();native_species=collections.defaultdict(set);seen=set()
 for key,raw in pairs(batch/'output/fastp'/f'{sample}_{a.native_arm}_fastp.fastq.gz'):
  q=key[0].decode().lstrip('@')
  if q in native:
   seqs=(key[1].decode(),key[3].decode());native_keys.add(seqs);seen.add(q)
   t=species(int(native[q]['aligner_taxid_lca']))
   if t is not None:native_species[t].add(seqs)
 assert seen==native.keys(),(sample,'Missing original FASTP pair')
 tier_keys=collections.defaultdict(set);tier_ids=collections.defaultdict(set);tier_species=collections.defaultdict(lambda:collections.defaultdict(set));statuses=collections.Counter()
 for q,m in meta[sample].items():
  v=viral[q];n=nr.get(q,{});ids=members[q];seqs=(m['query_seq'],m['query_seq_rev']);assert len(ids)==int(m['occurrences'])
  status=n.get('nr_status','pending');statuses[status]+=1
  tiers=['viral_reference_candidate']
  if status=='viral_only':
   tiers.append('nr_viral_only')
   if n['aligned_mates']=='2':tiers.append('nr_viral_only_both_mates')
   if n['taxonomic_relation'] in ['exact','ancestor_related']:tiers.append('nr_viral_only_related')
  for tier in tiers:
   tier_keys[tier].add(seqs);tier_ids[tier].update(ids);t=species(int(v['protein_taxid_lca']))
   if t is not None:tier_species[tier][t].add(seqs)
  pair_rows.append(dict(sample=sample,query_pair=q,exact_pair_also_native=seqs in native_keys,original_pairs=len(ids),added_original_pairs=len(ids-native.keys()),viral_reference_taxid=v['protein_taxid_lca'],viral_reference_aligned_mates=v['aligned_mates'],maximum_identity_pct=v['maximum_identity_pct'],maximum_alignment_aa=v['maximum_alignment_aa'],nr_status=status,nr_aligned_mates=n.get('aligned_mates',''),nr_taxid_lca=n.get('nr_taxid_lca',''),taxonomic_relation=n.get('taxonomic_relation','pending')))
 for tier in ['viral_reference_candidate']+(['nr_viral_only','nr_viral_only_both_mates','nr_viral_only_related'] if a.nr else []):
  keys=tier_keys[tier];ids=tier_ids[tier]
  sample_rows.append(dict(sample=sample,tier=tier,native_original_pairs=len(native),native_exact_pairs=len(native_keys),protein_exact_pairs=len(keys),added_exact_pairs=len(keys-native_keys),lost_exact_pairs=len(native_keys-keys),shared_exact_pairs=len(keys&native_keys),protein_original_pairs=len(ids),added_original_pairs=len(ids-native.keys()),lost_original_pairs=len(native.keys()-ids)))
  for minimum in [1,2,3]:
   base={t for t,s in native_species.items() if len(s)>=minimum};called={t for t,s in tier_species[tier].items() if len(s)>=minimum}
   species_rows.append(dict(sample=sample,tier=tier,minimum_exact_pairs=minimum,native_species=len(base),protein_species=len(called),gained_species=len(called-base),lost_species=len(base-called),native_human_species=sum(map(human,base)),protein_human_species=sum(map(human,called)),gained_human_species=sum(map(human,called-base)),lost_human_species=sum(map(human,base-called))))
   for t in sorted(base^called):species_changes.append(dict(sample=sample,tier=tier,minimum_exact_pairs=minimum,change='gained' if t in called else 'lost',taxid=t,name=annotations.get(t,{}).get('name',''),human_annotation=annotations.get(t,{}).get('infection_status_human','NA'),native_exact_pairs=len(native_species.get(t,set())),protein_exact_pairs=len(tier_species[tier].get(t,set()))))
for name,rows in [('sample-metrics.tsv',sample_rows),('pair-evidence.tsv',pair_rows),('species-metrics.tsv',species_rows),('species-changes.tsv',species_changes)]:
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['sample'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'samples':len(order),'native_arm':a.native_arm,'nr_complete':bool(a.nr),'tiers':{},'nr_status':dict(collections.Counter(x['nr_status'] for x in pair_rows)),'limits':['Broad protein evidence is not known-truth sensitivity or specificity.','Both arms use exact complete trimmed nucleotide-pair keys for this comparison; the native production duplicate policy also uses alignment coordinates and sequence similarity.','Species summaries use the frozen viral-reference taxonomy, including the original protein LCA, with exact-pair support; they do not establish infection or novel-virus discovery.','NR tiers and both-mate/related subsets are descriptive evidence summaries, not tuned calling thresholds.']}
for tier in sorted({x['tier'] for x in sample_rows}):
 rows=[x for x in sample_rows if x['tier']==tier];summary['tiers'][tier]={k:sum(x[k] for x in rows) for k in rows[0] if k not in ['sample','tier']}
 summary['tiers'][tier]['samples_with_gain']=sum(x['added_exact_pairs']>0 for x in rows);summary['tiers'][tier]['samples_with_loss']=sum(x['lost_exact_pairs']>0 for x in rows)
species_summary={'native_arm':a.native_arm,'samples':len(order),'counts':{},'interpretation':'Fixed viral-reference species names and human-infection annotations. Support uses exact complete-pair keys on both sides at descriptive thresholds1,2,3; it is not a diagnostic cutoff or an infection claim.'}
for tier in sorted({x['tier'] for x in species_rows}):
 for minimum in [1,2,3]:
  selected=[x for x in species_rows if x['tier']==tier and x['minimum_exact_pairs']==minimum]
  species_summary['counts'][f'{tier}/minimum_{minimum}']={k:sum(x[k] for x in selected) for k in selected[0] if k not in ['sample','tier','minimum_exact_pairs']}
(out/'species-summary.json').write_text(json.dumps(species_summary,indent=2))
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
