"""Measure subject-deduplication effects on the preserved earlier NT-panel outputs."""
import collections,csv,hashlib,json
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/protein-nt-panel';out=r/'cohort-storage/translated-core-nt/old-panel-filter-audit';out.mkdir(exist_ok=False)
def kept(values,deduplicate):
 if deduplicate:
  best={}
  for x in values:
   previous=best.get(x[1])
   if previous is None or (float(x[3]),int(x[5]))>(float(previous[3]),int(previous[5])):best[x[1]]=x
  values=list(best.values())
 scores=sorted({float(x[3]) for x in values},reverse=True)
 if not scores:return set(),None
 cutoff=max(.9*scores[0],scores[min(9,len(scores)-1)])
 return set().union(*({int(t) for t in x[2].split(';') if t.isdigit()} or {0} for x in values if float(x[3])>=cutoff)),cutoff
fixture=[['q','a','1',str(score),'0','100','100','100'] for score in range(100,90,-1)]+[['q','b','2','90','0','100','100','100']]
assert kept(fixture,False)[0]=={1} and kept(fixture,True)[0]=={1,2}
queries=[x[1:].split()[0] for x in (d/'queries.fna').read_text().splitlines() if x.startswith('>')];summary={};details=[]
for method in ['megablast','blastn']:
 path=d/(method+'.tsv');hits=collections.defaultdict(list)
 for x in csv.reader(path.open(),delimiter='\t'):assert len(x)==8;hits[x[0]].append(x)
 old={};new={};changed=[];cutoffs=0
 for q in queries:
  old[q],lc=kept(hits[q],False);new[q],rc=kept(hits[q],True);cutoffs+=lc!=rc
  if old[q]!=new[q]:changed.append(q);details.append({'method':method,'query':q,'old_taxids':';'.join(map(str,sorted(old[q]))),'production_order_taxids':';'.join(map(str,sorted(new[q]))),'old_minimum_bitscore':lc,'production_order_minimum_bitscore':rc})
 pairs={q.rsplit('/',1)[0] for q in queries};pair_changes=[]
 for q in pairs:
  if set().union(*(old.get(q+'/'+str(m),set()) for m in [1,2]))!=set().union(*(new.get(q+'/'+str(m),set()) for m in [1,2])):pair_changes.append(q)
 summary[method]={'queries':len(queries),'pairs':len(pairs),'changed_score_cutoffs':cutoffs,'changed_mate_taxon_sets':len(changed),'changed_pair_taxon_sets':len(pair_changes),'changed_pair_IDs':sorted(pair_changes),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
with (out/'differences.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['method','query','old_taxids','production_order_taxids','old_minimum_bitscore','production_order_minimum_bitscore'],delimiter='\t');w.writeheader();w.writerows(details)
result={'summary':summary,'mechanism_fixture_verified':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'production_filter_order_source_sha256':hashlib.sha256((r/'subworkflows/local/blastFasta/main.nf').read_bytes()).hexdigest(),'limits':['Existing raw E10 outputs only; no search repeated and no historical result overwritten.','Taxon-set parity is stronger than LCA parity; changed sets would require taxonomy-level analysis before reporting changed verdicts.','The earlier summary ranked all HSPs; production first retains one highest-score/longest HSP per query and subject.']};(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
