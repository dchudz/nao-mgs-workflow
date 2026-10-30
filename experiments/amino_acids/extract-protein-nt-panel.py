"""Copy selected mates byte-for-byte from the frozen NR query pool."""
import csv,hashlib,json
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/protein-nt-panel';source=r/'protein/combined-validation/combined-queries.fna'
protocol=json.loads((d/'protocol.json').read_text());assert hashlib.sha256((d/'selection.tsv').read_bytes()).hexdigest()==protocol['selection_sha256']
wanted={x['query_pair']+'/'+str(m) for x in csv.DictReader((d/'selection.tsv').open(),delimiter='\t') for m in [1,2]};assert len(wanted)==464
seen=set();h=hashlib.sha256();keep=False
with source.open('rb') as src,(d/'queries.fna').open('xb') as dst:
 for line in src:
  h.update(line)
  if line.startswith(b'>'):
   q=line[1:].split()[0].decode();keep=q in wanted
   if keep:assert q not in seen;seen.add(q)
  if keep:dst.write(line)
assert seen==wanted
assert h.hexdigest()=='6aa28bedbaf612c3a7c5a895b24a3be4f18d2c8264a439846427d943b289c597'
(d/'query-provenance.json').write_text(json.dumps(dict(source_sha256=h.hexdigest(),query_sha256=hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),queries=len(seen)),indent=2))
print((d/'query-provenance.json').read_text())
