"""Reuse the outcome-independent protein diagnostic panel without selecting new reads."""
import csv,datetime,hashlib,json,shutil
from pathlib import Path
r=Path('/home/ec2-user/experiment/cohort-storage');src=r/'protein-nt-panel';dst=r/'translated-core-nt/protein232';pilot=r/'translated-core-nt/pilot24'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((src/'protocol.json').read_text()); provenance=json.loads((src/'query-provenance.json').read_text())
assert sha(src/'queries.fna')==provenance['query_sha256']=='a98e68ff4763d7702a0755cbf79c551df48ea2e9fe01cc833e039d9c3f1cc794'
assert sha(src/'selection.tsv')==protocol['selection_sha256']=='bf7510b756f248fee0da885043b4de3a52910983d428282ef8a42e304a6c6bdf'
rows=list(csv.DictReader((src/'selection.tsv').open(),delimiter='\t'));assert len(rows)==232
queries=[x[1:].split()[0] for x in (src/'queries.fna').read_text().splitlines() if x.startswith('>')];assert len(queries)==len(set(queries))==464
assert {x.rsplit('/',1)[0] for x in queries}=={x['query_pair'] for x in rows}
dst.mkdir(exist_ok=False)
for name in ['queries.fna','protocol.json','query-provenance.json','selection.tsv']:
 shutil.copy2(src/name,dst/('original-'+name if name!='queries.fna' else name))
for row in rows:row['seq_id']=row['query_pair'];row['available_mates']=2
with (dst/'selection.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for name in ['reference-inventory.json','reference-info.txt']:shutil.copy2(pilot/name,dst/name)
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pairs':232,'queries':464,'query_sha256':sha(dst/'queries.fna'),'selection_sha256':sha(dst/'selection.tsv'),'reference_inventory_sha256':sha(dst/'reference-inventory.json'),'source_query_sha256':sha(src/'queries.fna'),'source_selection_sha256':sha(src/'selection.tsv'),'selection_rule':protocol['selection'],'limits':['Frozen one-pair-per-sample diagnostic panel, not a read-weighted accuracy estimate.','Original separate trimmed mates, not production mate merging.','Original selection predates broad-NR outcomes; choosing to reuse this panel occurs after those outcomes.','New nucleotide controls match16threads and full19field outputs; older8thread controls retained separately.','Translated search uses word4 threshold18 after the matched43mate pilot; defaultTBLASTX not scaled.','Evalues may change when retrieved HSP sets change; shared bitscore parity does not assert Evalue parity.']}
(dst/'input-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))
