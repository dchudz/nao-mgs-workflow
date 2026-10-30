from pathlib import Path
import csv,gzip,json,collections,functools
r=Path('/home/ec2-user/experiment');d=r/'e2e/SRR37076045';parents={};ranks={}
for line in (r/'reference/results/taxonomy-nodes.dmp').open():
 fields=[f.strip() for f in line.split('|')];parents[int(fields[0])]=int(fields[1]);ranks[int(fields[0])]=fields[2]
@functools.lru_cache(None)
def species(v):
 if v in (None,'NA',''):return None
 t=int(v);seen=set()
 while t in parents and t not in seen:
  if ranks[t]=='species':return t
  seen.add(t);t=parents[t]
 return None
with gzip.open(r/'reference/results/total-virus-db-annotated.tsv.gz','rt') as f:annotations={int(x['taxid']):x for x in csv.DictReader(f,delimiter='\t')}
sets={};support={}
for arm in ['production','validation_control']:
 candidate=collections.Counter();validated=collections.Counter();artificial=collections.Counter()
 with gzip.open(d/(arm+'-downstream/output/results_downstream')/f'{d.name}_validation_hits.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   sid=row['seq_id']
   if row['prim_align_dup_exemplar']!=sid or row['sim_dup_exemplar']!=sid:continue
   t=species(row['aligner_taxid_lca'])
   if t is None:continue
   candidate[t]+=1
   if row['validation_status']=='aligned' and species(row['validation_staxid_lca'])==t:
    validated[t]+=1
    if row.get('validation_staxid_lca_artificial') not in (None,'','NA'):artificial[t]+=1
 sets[arm]={'candidate':{t for t,c in candidate.items() if c>=3},'validated':{t for t,c in validated.items() if c>=3}}
 support[arm]={'candidate':candidate,'validated':validated,'artificial':artificial}
assert sets['production']['candidate']==sets['validation_control']['candidate']
changes=[]
for t in sorted(sets['production']['validated']^sets['validation_control']['validated']):
 note=annotations.get(t,{})
 changes.append({'species_taxid':t,'name':note.get('name',''),'human_annotation':note.get('infection_status_human','NA'),'change':'gained' if t in sets['validation_control']['validated'] else 'lost','production_validated_exemplars':support['production']['validated'][t],'control_validated_exemplars':support['validation_control']['validated'][t]})
record={'sample':d.name,'candidate_species_calls':len(sets['production']['candidate']),'validated_species_calls':{a:len(x['validated']) for a,x in sets.items()},'changes':changes,'rule':'At least three joint alignment-and-similarity duplicate exemplars, each with its own aligned validation natural LCA concordant at species level. Production merged/joined queries; no extrapolation to unsampled pairs.','scope':'Native alignment arms only; expanded-screen and protein arms pending.'}
(d/'native-call-comparison.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
