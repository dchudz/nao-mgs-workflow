"""Freeze a local locus-inspection panel; never use it for comprehensive validation."""
import collections,csv,gzip,hashlib,json
from pathlib import Path
r=Path('/home/ec2-user/experiment');src=r/'protein/species-gain-cases-compacted-v2';out=r/'protein/common-reference-case-audit-v2';out.mkdir(exist_ok=False)
table=lambda n:list(csv.DictReader((src/n).open(),delimiter='\t'))
calls={x['sample']:x for x in table('calls.tsv') if x['human_annotation']=='1'}
pairs={x['query_pair']:x for x in table('pairs.tsv') if x['sample'] in calls}
raw=[x for x in csv.reader((src/'selected-viral-alignments.tsv').open(),delimiter='\t') if x[0].rsplit('/',1)[0] in pairs and int(x[4])>=25]
best=collections.defaultdict(float)
for x in raw:best[x[0]]=max(best[x[0]],float(x[2]))
support=collections.defaultdict(lambda:collections.defaultdict(set));scores=collections.defaultdict(lambda:collections.defaultdict(dict))
for x in raw:
 if float(x[2])<.9*best[x[0]]:continue
 sample=pairs[x[0].rsplit('/',1)[0]]['sample'];genome=x[1].rsplit('~',2)[0];support[sample][genome].add(x[0]);scores[sample][genome][x[0]]=max(float(x[2]),scores[sample][genome].get(x[0],0))
chosen={s:min(gs,key=lambda g:(-len(gs[g]),-sum(scores[s][g].values()),g)) for s,gs in support.items()}
def fasta(path):
 name=None;seq=[]
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq)
    name=line[1:].split()[0];seq=[]
   else:seq.append(line.strip())
  if name is not None:yield name,''.join(seq)
wanted=set(chosen.values());refs={k:v for k,v in fasta(r/'reference/results/virus-genomes-masked.fasta.gz') if k in wanted};assert refs.keys()==wanted
queries={k:v for k,v in fasta(r/'protein/combined-validation/combined-queries.fna') if k.rsplit('/',1)[0] in pairs};assert len(queries)==2*len(pairs)
meta=[]
for s,g in sorted(chosen.items()):
 d=out/s;d.mkdir();(d/'reference.fna').write_text('>'+g+'\n'+refs[g]+'\n');qs={q:v for q,v in queries.items() if pairs[q.rsplit('/',1)[0]]['sample']==s};(d/'queries.fna').write_text(''.join('>'+q+'\n'+v+'\n' for q,v in sorted(qs.items())))
 meta.append({'sample':s,'species':calls[s]['species'],'reference':g,'reference_nt':len(refs[g]),'query_mates':len(qs),'existing_near_best_query_support':len(support[s][g]),'queries_sha256':hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),'reference_sha256':hashlib.sha256((d/'reference.fna').read_bytes()).hexdigest()})
(out/'protocol.json').write_text(json.dumps({'cases':meta,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference_choice':'Most audited mates with eligible near-best hits, then total best bitscore, then lexical accession. Chosen only from preserved viral-reference hits.','selection':'Previously frozen audit pairs for all four human-annotated concordant gains; twelve of eighteen similarity supports for Betacoronavirus and all pairs for other three calls.','limits':['Local coordinate-inspection only, never a replacement validation database or a source of changed calls.','Small-database E-values are not comparable to comprehensive-search E-values.','Failure to align to this one reference does not negate the broad-reference evidence.','Overlapping intervals are not a molecular-duplicate determination.']},indent=2));print(json.dumps(meta))
