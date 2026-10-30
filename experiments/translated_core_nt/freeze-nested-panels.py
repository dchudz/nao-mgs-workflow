"""Freeze nested query-count panels without reading any alignment output."""
import csv, hashlib, importlib.util, json, shutil
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
root=m.ROOT/'cohort-storage/translated-core-nt';pilot=root/'pilot24';source=m.ROOT/'cohort-storage/seed-diagnostic200'
original=list(csv.DictReader((pilot/'selection.tsv').open(),delimiter='\t'));parent=list(csv.DictReader((source/'selection.tsv').open(),delimiter='\t'));seqs=m.fasta(source/'queries.fna');used={x['seq_id'] for x in original};large=list(original)
flags=['baseline_control','AA12_added','NT20_added']
for flag in flags:
    present=[x for x in large if x['stratum']==flag];samples={x['sample'] for x in present};count=len(present)
    candidates=sorted((x for x in parent if x[flag]=='True' and x['seq_id'] not in used),key=lambda x:hashlib.sha256(('same-reference-expansion-v1:'+flag+':'+x['seq_id']).encode()).digest())
    for row in candidates:
        if row['sample'] in samples:continue
        large.append(dict(row,stratum=flag));samples.add(row['sample']);used.add(row['seq_id']);count+=1
        if count==32:break
    assert count==32
small=[next(x for x in original if x['stratum']==flag) for flag in flags]
for name,rows in [('small3',small),('expanded96',large)]:
    d=root/name;d.mkdir(exist_ok=False)
    with (d/'selection.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(original[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    queries={q:s for row in rows for q,s in seqs.items() if q.rsplit('/',1)[0]==row['seq_id']}
    (d/'queries.fna').write_text(''.join('>'+q+'\n'+s+'\n' for q,s in queries.items()))
    for file in ['reference-inventory.json','reference-info.txt']:shutil.copy2(pilot/file,d/file)
    for q,s in m.fasta(pilot/'queries.fna').items():
        if q in queries:assert queries[q]==s
    if name=='expanded96':assert set(m.fasta(pilot/'queries.fna')).issubset(queries)
    else:assert set(queries).issubset(m.fasta(pilot/'queries.fna'))
    manifest=dict(created_utc=m.now(),pairs=len(rows),queries=len(queries),query_sha256=m.sha(d/'queries.fna'),selection_sha256=m.sha(d/'selection.tsv'),reference_inventory_sha256=m.sha(d/'reference-inventory.json'),source_query_sha256=m.sha(source/'queries.fna'),source_selection_sha256=m.sha(source/'selection.tsv'),pilot_query_sha256=m.sha(pilot/'queries.fna'),selection_rule='small3: first pilot pair in each stratum. expanded96: retain all24pilot pairs; append24per stratum, unique samples within each stratum, hash-ranked from same parent panel, excluding already selected pairs. No alignment output read.',limits=['Nested diagnostic panels, not prevalence estimates','Changing query batches can change search behavior; check shared-read output and assignment parity','Execution requires completed pilot runtime review; freezing does not launch searches'])
    m.save(d/'input-manifest.json',manifest);print(json.dumps(manifest))
