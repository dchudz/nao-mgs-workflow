"""Incomplete, inconsistent, or unpreserved results must not enter summaries."""
import json
from pathlib import Path
import shutil
import tempfile
from aggregate_reports import collect, read, summarize

p=Path(__file__).parent
selection=read(p/'selection.json');orders=read(p/'validation-order.json')
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);case=root/'SRR38294967_5m';case.mkdir()
    original=p/'results'/case.name
    for name in ('complete.json','synthetic-control-gate.json','bowtie-output-integrity.json','completed-evidence-backup.json'):
        shutil.copyfile(original/name,case/name)
    report=case/'comparison-report';(report/'standard').mkdir(parents=True)
    for name in ('accounting.json','integrity.json','validation-status-counts.json'):
        shutil.copyfile(original/'comparison-report'/name,report/name)
    shutil.copyfile(original/'comparison-report/standard/viral_read_status.tsv',report/'standard/viral_read_status.tsv')
    rows=collect(root,selection,orders)
    assert sum(r['status']=='complete' for r in rows)==1
    assert rows[0]['production_no_alignment']==0
    assert all(r['status']=='not_collected' for r in rows[1:])
    for path,key,value in ((report/'accounting.json','shared_output_pair_ids',1),
                           (case/'completed-evidence-backup.json','readback_verified',False),
                           (report/'integrity.json','passed',False)):
        text=path.read_text();data=json.loads(text);data[key]=value;path.write_text(json.dumps(data))
        rows=collect(root,selection,orders)
        assert rows[0]['status']=='blocked'
        assert not summarize(rows)['strata'][0].get('completed_run_count')
        path.write_text(text)
    (report/'accounting.json').unlink()
    assert collect(root,selection,orders)[0]['status']=='incomplete'
print('Report acceptance and missing-data checks passed')
