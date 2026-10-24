"""Export only reports from completed, privately preserved comparison units."""
import hashlib
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
        for path in sorted(d.rglob('*')):
            if not path.is_file() or path.is_symlink():continue
            relative=path.relative_to(d)
            if 'work' in relative.parts or '.nextflow' in relative.parts:continue
            if any(part.startswith('.') for part in relative.parts):continue
            if path.suffix=='.json' or path.name=='trace.tsv' or ('comparison-report' in relative.parts and path.suffix in {'.tsv','.log'}):
                tar.add(path,arcname=str(path.relative_to(root/'comparisons')))
print(json.dumps({'archive':str(target),'labels':labels,'bytes':target.stat().st_size,
                  'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))
