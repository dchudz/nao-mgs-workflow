"""Require full output-byte parity on the saved real-HSP integration fixture."""
import gzip,hashlib,json,subprocess,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/parallel-adapter';D.mkdir(exist_ok=True)
implementation_sha=hashlib.sha256((R/'analysis-tools/parallel_diamond_adapter.py').read_bytes()).hexdigest()
prior=json.loads((R/'integration/adapter-bounded/verification.json').read_text());records=[]
for workers in [1,32]:
    for mode,key,expected in [('extract','bounded_extract',prior['extract']['decompressed_sha256']),('final','bounded_final',prior['final']['decompressed_sha256'])]:
        command=prior[key]['command'].copy();command[1]=str(R/'analysis-tools/parallel_diamond_adapter.py')
        output=D/f'{mode}-{workers}.gz';command[command.index('--output')+1]=str(output);command+=['--workers',str(workers)]
        start=time.monotonic()
        with (D/f'{mode}-{workers}.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
        elapsed=time.monotonic()-start
        with gzip.open(output,'rb') as f:sha=hashlib.sha256(f.read()).hexdigest()
        assert sha==expected,(mode,workers,sha,expected)
        record={'mode':mode,'workers':workers,'seconds':elapsed,'decompressed_sha256':sha,'matches_frozen_output':True,'command':command};records.append(record);print(json.dumps(record),flush=True)
assert implementation_sha==hashlib.sha256((R/'analysis-tools/parallel_diamond_adapter.py').read_bytes()).hexdigest()
receipt={'verified':True,'implementation_sha256':implementation_sha,'frozen_adapter_sha256':hashlib.sha256((R/'experiments/amino_acids/diamond_adapter.py').read_bytes()).hexdigest(),'eligible_HSPs':prior['eligible_HSPs'],'scope':'Original real-HSP fixture and original FASTQ qualities; candidate FASTQ and final complete hit rows are byte-identical to the frozen adapter for both worker counts.','results':records}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
