"""Freeze a sample-balanced protein-only diagnostic before broad-reference outcomes."""
import collections,csv,datetime,hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parent;d=r/'results/protein-nt-panel';d.mkdir(exist_ok=True)
a={x['query_pair']:x for x in csv.DictReader((r/'results/protein200/pair-evidence.tsv').open(),delimiter='\t')}
b={x['query_pair']:x for x in csv.DictReader((r/'results/protein200-nt20/pair-evidence.tsv').open(),delimiter='\t')}
assert a.keys()==b.keys() and all(x['nr_status']=='pending' for x in a.values())
groups=collections.defaultdict(lambda:collections.defaultdict(list))
for q,x in a.items():
 if x['exact_pair_also_native']==b[q]['exact_pair_also_native']=='False':group='added_beyond_NT24_and_NT20'
 elif x['exact_pair_also_native']==b[q]['exact_pair_also_native']=='True':group='shared_with_NT24_and_NT20'
 else:continue
 groups[group][x['sample']].append(q)
def rank(q):return hashlib.sha256(('protein-nt-panel:'+q).encode()).digest()
rows=[]
for group,samples in sorted(groups.items()):
 chosen=sorted(samples,key=lambda s:rank(group+':'+s))
 if group.startswith('shared'):chosen=chosen[:32]
 for s in chosen:
  q=min(samples[s],key=rank);x=a[q]
  rows.append(dict(query_pair=q,sample=s,stratum=group,sample_stratum_candidates=len(samples[s]),viral_reference_taxid=x['viral_reference_taxid'],viral_reference_maximum_identity_pct=x['maximum_identity_pct']))
assert len(rows)==232 and len({x['query_pair'] for x in rows})==232
with (d/'selection.tsv').open('x') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
protocol=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),secondary=True,selection='One SHA256-ranked candidate absent from both NT24 and NT20 per each of the fixed first200 samples; one shared candidate per32 SHA256-ranked sample controls. Selection made while every broad-NR outcome was pending.',selected_pairs=len(rows),selected_mates=464,strata=dict(collections.Counter(x['stratum'] for x in rows)),selection_sha256=hashlib.sha256((d/'selection.tsv').read_bytes()).hexdigest(),query_representation='Exactly the same complete trimmed mate sequences already included in the frozen full NR query pool. No query or protein search is repeated or altered.',searches=['megablast','blastn -task blastn'],reference='Existing comprehensive core_nt snapshot, unchanged.',resources='Existing primary worker only, after the sampling-cap benchmark completes,8threads with explicit8192DockerCPUshares; no added instance. Six-hour per-search resource bound and overall $500 account cap.',limits=['One candidate per sample is not a read-weighted accuracy estimate.','Shared controls were selected independently of any NR outcome.','NR and core_nt differ in reference representation and coverage; this comparison cannot isolate algorithm alone.','Both nucleotide tasks use identical separate mates; original production validation merges/joins mates.','Search tasks change both seed and scoring defaults.','The primary cohort, thresholds and original NR search remain unchanged.'])
(d/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');print(json.dumps(protocol))
