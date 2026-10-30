"""Exercise adapter output through the unchanged production duplicate executable."""
import csv,gzip,hashlib,json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from diamond_adapter import FIELDS,paired_row
from test_diamond_adapter import hit,PARENTS

binary=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
rows=[]
for name,hs,quality in [
    ('same_low',[[hit()],[hit(query='p/2',frame=-1)]],'?'),
    ('same_high',[[hit()],[hit(query='p/2',frame=-1)]],'I'),
    ('single',[[hit()],[]],'I'),
    ('split_za',[[hit(subject='z~0~30')],[hit(subject='a~3~30',query='p/2')]],'I'),
    ('split_az',[[hit(subject='a~3~30')],[hit(subject='z~0~30',query='p/2')]],'?'),
]:
    rows.append(paired_row(name,'synthetic',('A'*21,'C'*21),(quality*21,quality*21),hs,PARENTS,{81077,99},set()))
source=out/'input.tsv.gz';dest=out/'marked.tsv.gz';meta=out/'groups.tsv.gz'
with gzip.open(source,'wt') as f:
    w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)
cmd=[str(binary),'-i',str(source),'-o',str(dest),'-m',str(meta),'-d','1','-n','1']
result=subprocess.run(cmd,capture_output=True,text=True,check=True)
(out/'stderr.log').write_text(result.stderr)
with gzip.open(dest,'rt') as f:observed={r['seq_id']:r for r in csv.DictReader(f,delimiter='\t')}
assert set(observed)=={r['seq_id'] for r in rows}
assert observed['same_low']['prim_align_dup_exemplar']=='same_high'
assert observed['same_high']['prim_align_dup_exemplar']=='same_high'
assert observed['single']['prim_align_dup_exemplar']=='single'
assert observed['split_az']['prim_align_dup_exemplar']=='split_za'
for row in rows:
    for name,value in row.items():assert observed[row['seq_id']][name]==str(value),(row['seq_id'],name)
receipt={'verified':True,'command':cmd,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
    'invariants':['Original columns and all five pair identities preserved','Higher original quality chooses the exemplar',
                  'One-mate evidence does not collapse into two-mate evidence','Split reference/mate order canonicalizes in production code'],
    'scope':'Synthetic adapter-to-production-coordinate-marker check; not a real DIAMOND alignment or baseline replay.'}
(out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
