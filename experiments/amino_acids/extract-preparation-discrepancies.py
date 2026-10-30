"""Preserve exact AA12 FASTP records for native assignments missing from protein preparation."""
import csv,hashlib,json,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
source=r/'cohort-storage/protein200/native-losses/lost-original-pairs.tsv';d=r/'cohort-storage/preparation-discrepancies';d.mkdir(exist_ok=True)
selected={x['seq_id']:x for x in csv.DictReader(source.open(),delimiter='\t') if x['loss_boundary']=='absent_from_prepared_protein_queries' and x['native_AA12_assignment_preserved']=='True'};assert len(selected)==8
batches={s:m.parent for m in (r/'cohort-storage/batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']};seen=set();records=[]
with (d/'input.fastq').open('xb') as out:
 for sample in sorted({x['sample'] for x in selected.values()}):
  path=batches[sample]/'output/fastp'/f'{sample}_aa12_fastp.fastq.gz'
  for key,raw in pairs(path):
   q=key[0].decode().lstrip('@')
   if q in selected:
    assert q not in seen;seen.add(q);out.write(raw);lines=raw.splitlines();records.append(dict(seq_id=q,sample=sample,source=str(path.relative_to(r)),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),record_sha256=hashlib.sha256(raw).hexdigest(),mate1_header=lines[0].decode(),mate2_header=lines[4].decode(),mate1_nt=len(key[1]),mate2_nt=len(key[3])))
assert seen==selected.keys();(d/'selected.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
