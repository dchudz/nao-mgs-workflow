import importlib.util,json,subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('translation',Path(__file__).with_name('characterize-query-translation.py'));t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);m=t.m
f=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=f/'score-statistics';d.mkdir(exist_ok=False);p=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE];report={}
for name in ['tblastx','query-tblastn-fast']:
 cfg=json.loads(Path(__file__).with_name(name+'.json').read_text());query=f/('query-only/queries.faa' if name.startswith('query') else 'queries.fna');cmd=p+[cfg['executable']]+cfg['arguments']+['-query',str(query),'-db',str(f/'db'),'-evalue','10','-max_hsps','5','-num_alignments','250']
 for fmt in ['0','6 qseqid sseqid qframe sframe score bitscore evalue length pident']:
  result=subprocess.run(cmd+['-outfmt',fmt],capture_output=True,text=True,check=True,timeout=60);(d/(name+('.text' if fmt=='0' else '.scores.tsv'))).write_text(result.stdout);(d/(name+'.stderr')).write_text(result.stderr)
 report[name]=dict(command=cmd)
m.save(d/'protocol.json',dict(created_utc=m.now(),methods=report,source_sha256=m.sha(__file__),purpose='Inspect raw score and statistical parameters; do not assume a universal Evalue correction between different programs'))
