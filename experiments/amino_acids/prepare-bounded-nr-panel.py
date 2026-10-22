"""Freeze a small broad-reference panel before inspecting NR outcomes."""
import collections,csv,datetime,hashlib,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/bounded-validation';d.mkdir(exist_ok=True)
source=r/'protein/combined-validation'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp','s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/protein200-native-comparison/pair-evidence.tsv',str(d/'candidate-evidence.tsv'),'--only-show-errors'],check=True)
groups=collections.defaultdict(list)
for x in csv.DictReader(open(d/'candidate-evidence.tsv'),delimiter='\t'):
 identity=float(x['maximum_identity_pct']);bin='below80' if identity<80 else '80to95' if identity<95 else 'atleast95'
 groups[x['sample'],x['exact_pair_also_native'],bin].append(x)
selected=[]
for key,rows in sorted(groups.items()):
 chosen=sorted(rows,key=lambda x:hashlib.sha256(('bounded-nr-panel:'+x['query_pair']).encode()).digest())[:2]
 for x in chosen:selected.append(dict(x,identity_bin=key[2],stratum_available_pairs=len(rows)))
with open(d/'selected.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t');w.writeheader();w.writerows(selected)
with open(d/'selected-viral-pairs.tsv','w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['query_pair','protein_taxid_lca']);w.writerows((x['query_pair'],x['viral_reference_taxid']) for x in selected)
controls=list(csv.DictReader(open(source/'native-controls.tsv'),delimiter='\t'))
seeds=list(csv.DictReader(open(source/'seed-selection.tsv'),delimiter='\t'))
wanted={x['query_pair'] for x in selected}|{x['query_pair'] for x in controls}|{'SEED_'+x['seq_id'] for x in seeds};seen=set();queries=0
with open(source/'combined-queries.fna') as f,open(d/'queries.fna','w') as out:
 while True:
  h=f.readline()
  if not h:break
  seq=f.readline();assert h.startswith('>') and seq and not seq.startswith('>')
  q=h[1:].strip();pair=q.rsplit('/',1)[0]
  if pair in wanted:out.write(h+seq);seen.add(pair);queries+=1
assert wanted==seen
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_search_unchanged':True,'purpose':'A separate bounded evidence panel if the full817744-query search cannot finish before the deadline. No verdict is extrapolated to unsampled candidates.','selection':'Up to2 SHA-ranked candidate pairs per first200 sample, crossed by existing native exact-pair overlap and viral-alignment maximum identity bins below80,80to95,atleast95. Include every already-frozen native control and BLAST diagnostic pair. No NR outcome used.','candidate_pairs':len(selected),'native_control_pairs':len(controls),'seed_diagnostic_pairs':len(seeds),'pairs':len(wanted),'query_reads':queries,'query_sha256':hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),'strata':'Selection is balanced across strata, not an unbiased cohort prevalence estimate; report exact sampled denominators and do not extrapolate total validated gains.','query_order':'Filtered from the original complete query pool, preserving retained order.','reference':'Same full comprehensive NR snapshot, frozen sensitivity and4index chunks.'}
(d/'protocol.json').write_text(json.dumps(protocol,indent=2));print(json.dumps(protocol),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/bounded-nr-panel/','--recursive','--only-show-errors'],check=True)
