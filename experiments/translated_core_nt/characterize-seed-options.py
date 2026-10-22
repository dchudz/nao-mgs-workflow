"""Record accepted tblastx seed settings with the fixed synthetic reference."""
import importlib.util,json,subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
d=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';out=d/'seed-option-characterization.json';assert not out.exists();rows=[]
for size,threshold in [(3,None),(4,None),(4,18),(5,20),(6,21)]:
    args=['tblastx','-word_size',str(size)]
    if threshold is not None:args+=['-threshold',str(threshold)]
    command=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE]+args+['-db',str(d/'db'),'-query',str(d/'queries.fna'),'-outfmt','6 qseqid sseqid pident length evalue']
    result=subprocess.run(command,capture_output=True,text=True,timeout=30)
    rows.append(dict(word_size=size,threshold=threshold,command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,synonymous80aa_exact_hit=any(x.startswith('synonymous\treference\t100.000\t80\t') for x in result.stdout.splitlines())))
assert all(x['returncode']==0 and x['synonymous80aa_exact_hit'] for x in rows[:3])
assert all(x['returncode']!=0 for x in rows[3:])
m.save(out,dict(created_utc=m.now(),rows=rows,excluded_from_real_data=True));print(json.dumps([{k:v for k,v in x.items() if k not in ['command','stdout']} for x in rows]))
