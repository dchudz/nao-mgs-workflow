"""Select already-measured NR taxon sets; no new search or verdict extrapolation."""
import datetime,hashlib,json
from pathlib import Path
r=Path('/home/ec2-user/experiment/protein');assert (r/'seed-nr-standalone/complete.json').exists(),'Run after the measured search'
d=r/'protein232-nr-evidence';d.mkdir(exist_ok=False);q=r/'protein232-queries.fna';assert hashlib.sha256(q.read_bytes()).hexdigest()=='a98e68ff4763d7702a0755cbf79c551df48ea2e9fe01cc833e039d9c3f1cc794'
queries={line[1:].split()[0] for line in q.read_bytes().splitlines() if line.startswith(b'>')};assert len(queries)==464
source=r/'nr-query-call-full/query-calls.tsv';h=hashlib.sha256();seen=set();output=d/'query-calls.tsv'
with source.open('rb') as inp,output.open('wb') as out:
 for line in inp:
  h.update(line);name=line.split(b'\t',1)[0]
  if name in queries:assert name not in seen;seen.add(name);out.write(line)
assert h.hexdigest()=='c76cf68eed3cc6f2bdaf8f6788b2e4fe6a1a2c4f9c76853e4bf51114cd4e3f42'
result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queries':464,'measured_eligible_queries':len(seen),'missing_query_IDs':sorted(x.decode() for x in queries-seen),'source_sha256':h.hexdigest(),'output_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'interpretation':'Exact subset of previously verified NR near-best taxon cache: E<=1e-5, at least25AA, within90percent of best eligible score per mate. Absent queries remain in the denominator.'};(d/'complete.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
