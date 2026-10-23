"""Describe score differences on canonical-AA HSPs; do not recalibrate searches."""
import argparse,collections,csv,gzip,hashlib,json,math,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('shared_hsps',type=Path);p.add_argument('matrix',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
matrix_text=a.matrix.read_text().split('/*A*/',1)[1].split('}',1)[0]
values=list(map(int,re.findall(r'-?\d+',re.sub(r'/\*.*?\*/','',matrix_text,flags=re.S))))
assert len(values)==26*26
alphabet='ARNDCQEGHILKMFPSTWYVBJZX*';canonical=set(alphabet[:20])
matrix={(q,s):values[i*26+j] for i,q in enumerate(alphabet) for j,s in enumerate(alphabet)}
assert all(matrix[q,s]==matrix[s,q] for q in canonical for s in canonical)
counts=collections.Counter();examples=[];deviations=[]
with gzip.open(a.shared_hsps,'rt') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        is_canonical=set(row['qseq']+row['sseq'])<=canonical
        native=float(row['tblastx_bitscore']);fast=float(row['fast_bitscore'])
        counts[('canonical' if is_canonical else 'contains_noncanonical')+('/different' if native!=fast else '/equal')]+=1
        if not is_canonical:continue
        assert len(row['qseq'])==len(row['sseq'])==int(row['length'])
        score=sum(matrix[q,s] for q,s in zip(row['qseq'],row['sseq']))
        standard=(.3176*score-math.log(.134))/math.log(2)
        deviations.append((abs(standard-native),abs(standard-fast)))
        if native!=fast and len(examples)<5:
            examples.append({k:row[k] for k in ['qseqid','qframe','sseqid','qseq','sseq','tblastx_bitscore','fast_bitscore']})
assert sum(counts.values())==28285
result=dict(counts=dict(counts),examples=examples,unmasked_matrix_probe_not_score_reconstruction=dict(canonical_HSPs=len(deviations),assumed_lambda=.3176,assumed_K=.134,native_max_absolute_bit_deviation=max(x[0] for x in deviations),fast_max_absolute_bit_deviation=max(x[1] for x in deviations),native_above_half_bit_deviation=sum(x[0]>.51 for x in deviations),fast_above_half_bit_deviation=sum(x[1]>.51 for x in deviations)),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.shared_hsps,a.matrix]},limits=['Characterization of existing matched HSPs only; no new search or altered threshold.','Canonical means every reported aligned character is one of the20standard amino acids. This does not assert absence of masking elsewhere in the query.','The unmasked matrix probe uses the pinned DIAMOND BLOSUM62 table and standard ungapped constants. It does not reconstruct either NCBI search scoring context or masking, and fails to reproduce their emitted scores.','The frozen output schema lacks raw scores and per-context statistical parameters. No causal mechanism or score-correction formula is established.','Removing stops or ambiguous characters alone cannot establish score parity on the observed canonical alignments.'])
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result['counts']))
