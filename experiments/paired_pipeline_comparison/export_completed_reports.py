"""Export only reports from completed, privately preserved comparison units."""
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

root=Path('/home/ec2-user/paired-study-20261024')
target=root/('report-export-'+sys.argv[1]+'.tar.gz')
labels=[]
with tarfile.open(target,'w:gz') as tar:
    for proof_path in sorted((root/'comparisons').glob('*/completed-evidence-backup.json')):
        d=proof_path.parent
        assert json.loads(proof_path.read_text())['readback_verified']
        assert json.loads((d/'comparison-report/accounting.json').read_text())['completed']
        labels.append(d.name)
        sample=d.name.rsplit('_',1)[0]
        inputs=root/'inputs'/sample
        acquisition_path=inputs/'verified.json'
        acquisition=json.loads(acquisition_path.read_text())
        downloads=[json.loads(p.read_text()) for p in (root/'batch-logs').glob('*/'+sample+'-download.json')]
        summary={key:acquisition[key] for key in ('sample','population_pairs','sampling','full_quality','source_sha256',
            'seed_sha256','files','seconds_excluding_initial_download','vdb_validation_seconds','completed_utc')}
        summary.update(source_bytes=acquisition['source_metadata']['ContentLength'],
            acquisition_receipt_sha256=hashlib.sha256(acquisition_path.read_bytes()).hexdigest(),
            download_receipts=downloads,download_time_recorded_in_batch=bool(downloads),
            download_seconds=sum(r['seconds'] for r in downloads) if downloads else None,
            backup_record=json.loads((inputs/'backup-readback.json').read_text()),
            interpretation='Acquisition covers the complete source archive. It is shared by shallow/deep comparisons and both arms; do not charge it once per arm or once per depth. Missing download timing is unknown, not zero.')
        payload=(json.dumps(summary,indent=2)+'\n').encode()
        item=tarfile.TarInfo(d.name+'/input-acquisition.json');item.size=len(payload)
        tar.addfile(item,io.BytesIO(payload))
        for path in sorted(d.rglob('*')):
            if not path.is_file() or path.is_symlink():continue
            relative=path.relative_to(d)
            if 'work' in relative.parts or '.nextflow' in relative.parts:continue
            if any(part.startswith('.') for part in relative.parts):continue
            if path.suffix=='.json' or path.name=='trace.tsv' or ('comparison-report' in relative.parts and path.suffix in {'.tsv','.log'}):
                tar.add(path,arcname=str(path.relative_to(root/'comparisons')))
print(json.dumps({'archive':str(target),'labels':labels,'bytes':target.stat().st_size,
                  'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))
