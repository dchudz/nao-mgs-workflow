"""Reproduce the saved protein FASTA from the exact current Bowtie2 genomes."""
import hashlib,json,subprocess,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/translation-reproduction';D.mkdir(exist_ok=True)
source=R/'analysis-tools/translate-alignment-reference.original.py'
expected='60f07af878602f9cd6d596616afbfa20746401b5aa33b1393f2364aab0c53a46'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
assert sha(source)==expected
root="r=Path('/home/ec2-user/experiment')";code=source.read_text();assert code.count(root)==1
script=D/'translate-isolated.py';script.write_text(code.replace(root,f'r=Path({str(D)!r})'))
(D/'reference/results').mkdir(parents=True,exist_ok=True);(D/'protein').mkdir(exist_ok=True)
genomes=R/'reference/results/virus-genomes-masked.fasta.gz';link=D/'reference/results'/genomes.name
if not link.exists():link.symlink_to(genomes)
assert link.resolve()==genomes
python=D/'.venv/bin/python';start=time.monotonic()
with (D/'translation.log').open('w') as log:subprocess.run([str(python),str(script)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
record=json.loads((D/'protein/alignment-translation.json').read_text())
actual=sha(D/'protein/viral-alignment.faa');original=sha(R/'protein/viral-alignment.faa')
assert actual==original=='e9861f81d94d567905610e544d906db2e3bc75a79a24d220247e6121a6b18966'
assert actual==record['sha256']
receipt={'verified':True,'original_builder_sha256':expected,'isolated_script_sha256':sha(script),
 'only_source_edit':'Output root changed to an isolated reproduction directory; input is a symlink to the current frozen Bowtie2 DUST-masked FASTA.',
 'genome_fasta_gzip_sha256':sha(genomes),'protein_fasta_sha256':actual,
 'restored_diamond_index_sha256':sha(R/'protein/viral-alignment.dmnd'),
 'wall_seconds':time.monotonic()-start,'translation':record,
 'dependencies':subprocess.check_output([str(python),'-m','pip','freeze'],text=True).splitlines(),
 'scope':'Complete six-frame protein FASTA regenerated from the current native alignment genomes. Exact bytes match the archived FASTA used to build the restored DIAMOND index; no scientific index or pipeline input is replaced.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
