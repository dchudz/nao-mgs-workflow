"""Retrieve one already-returned target; no new homology search or verdict."""
import datetime,hashlib,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/seed-nr-forced-double-index-trace/output-audit-v2'
audit=json.loads((d/'complete.json').read_text());assert audit['comparisons']['pooled_double']['alignment_multiset_equal']
index=json.loads((d/'human-target-location.json').read_text());assert len(index['targets'])==1
target=index['targets'][0];assert target['global_ordinal']==815174796
output=d/'human-target-815174796.faa';assert not output.exists()
cmd=[str(r/'protein/tools/diamond'),'getseq','--db',target['database'],'--seq',str(target['diamond_getseq_number']),'--out',str(output)]
started=time.monotonic()
with (d/'human-target-retrieval.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=300)
text=output.read_text();assert text.count('>')==1 and text.startswith('>')
header,*lines=text.splitlines();sequence=''.join(lines);assert sequence and set(sequence)<=set('ABCDEFGHIJKLMNOPQRSTUVWXYZ*')
receipt=dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=cmd,wall_seconds=time.monotonic()-started,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),target=target,header=header,amino_acids=len(sequence),sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),fasta_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),limits=['Exact target retrieval after the comprehensive NR comparison, not restricted-reference validation.','All four recovered human HSPs refer to this ordinal. A human taxonomic label does not establish the reads originate from humans.','Numeric identity is tied to the frozen NR snapshot. No accession or gene name is inferred.'])
(d/'human-target-retrieval.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
