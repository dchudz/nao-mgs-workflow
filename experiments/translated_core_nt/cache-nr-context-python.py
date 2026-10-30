"""Use the unchanged Python reader for a small local subset when macOS cannot build the Linux compactor."""
import argparse, csv, hashlib, importlib.util, json, platform, time
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('classifier',type=Path);a=p.parse_args()
spec=importlib.util.spec_from_file_location('original_nr_classifier',a.classifier)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
manifest=json.loads((a.evidence/'download-receipt.json').read_text())
known={x['name']:x for x in manifest['files']}
for label,name,expected_rows in [('pooled','pooled-seed-alignments.tsv',5942196),('standalone','standalone-alignments.tsv',5678862)]:
    source=a.evidence/name;d=a.evidence/(label+'-cache');d.mkdir(exist_ok=False)
    h=hashlib.sha256()
    with source.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    assert h.hexdigest()==known[name]['sha256']
    count=[0]
    def rows(reader):
        for row in reader:
            count[0]+=1
            yield row
    start=time.monotonic()
    with source.open() as f:calls=module.near_best_query_taxa(rows(csv.reader(f,delimiter='\t')))
    assert count[0]==expected_rows
    target=d/'query-calls.tsv'
    with target.open('w') as f:
        for q,(taxa,score) in sorted(calls.items()):f.write(q+'\t'+str(score)+'\t'+';'.join(map(str,sorted(taxa)))+'\n')
    result=dict(input_sha256=h.hexdigest(),input_rows=count[0],query_rows=len(calls),wall_seconds=time.monotonic()-start,output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),classifier_sha256=hashlib.sha256(a.classifier.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),platform=platform.platform(),reader='Unchanged original near_best_query_taxa Python function',reason='The macOS C++ standard library does not provide the floating-point from_chars overload used by the Linux compactor. No compactor code was changed; existing Python aggregation is sufficient for these two small files.')
    (d/'complete.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
