"""Characterize record selection using only generated strings and invented labels."""
import collections,csv,gzip,hashlib,json,os,random,subprocess,time
from pathlib import Path
ROOT=Path('/Users/david/Documents/ChatGPT/pipeline/.tmp/duplicate-characterization-20261024')
BIN=Path('/Users/david/Documents/ChatGPT/pipeline/.tmp/anchor-tools/source/target/release')
FIELDS=['seq_id','prim_align_genome_id_all','prim_align_ref_start_unclipped','prim_align_ref_start_unclipped_rev','prim_align_ref_end_unclipped','prim_align_ref_end_unclipped_rev','prim_align_query_rc','prim_align_query_rc_rev','query_seq','query_seq_rev','query_qual','query_qual_rev','synthetic_category']
PROTOCOL={'seed':284104,'data':'Generated uniform random strings; fictional reference keys and categories X/Y only. No biological sample/reference used.','source_commit':os.environ.get('CHARACTERIZATION_SOURCE_COMMIT'),'groups':[3,30],'random_replicates_per_cell':30,'random_removal_probabilities':[0.1,0.5,0.9],'quality_models':['independent','X_higher','Y_higher'],'curated_sizes':[1,3,10,100,1000],'curated_cases':['category_swap','total_increase','total_decrease'],'interpretation':'Synthetic frequencies apply only to this declared generator, not to experimental reads, biological prevalence, or the published 53 calls. Empty input has zero selected records by definition; tool empty-input behavior is not tested.','production_code_modified':False}
(ROOT/'protocol.json').write_text(json.dumps(PROTOCOL,indent=2)+'\n')
def row(identifier,key,seqs,q,label):
 return dict(zip(FIELDS,[identifier,key,'100','300','260','460','False','True',*seqs,chr(q+33)*160,chr(q+33)*160,label]))
def strings(rng):return [''.join(rng.choice('ACGT') for _ in range(160)) for _ in range(2)]
def evaluate(name,rows):
 p=ROOT/'cases'/name;p.mkdir(parents=True,exist_ok=False)
 f=p/'input.tsv.gz'
 with gzip.open(f,'wt') as h:
  w=csv.DictWriter(h,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda r:r['seq_id']))
 if not rows:return {'records':0,'joint':0,'X':0,'Y':0,'empty_defined':True}
 cmds=[[str(BIN/'mark_duplicates'),'-i',str(f),'-o',str(p/'coordinate.tsv.gz'),'-m',str(p/'groups.tsv'),'--num-threads','1'],[str(BIN/'mark_duplicates_similarity'),'-i',str(p/'coordinate.tsv.gz'),'-o',str(p/'final.tsv.gz')]]
 for i,c in enumerate(cmds):
  with (p/f'{i}.log').open('w') as h:subprocess.run(c,stdout=h,stderr=subprocess.STDOUT,check=True,timeout=60,env={**os.environ,'RAYON_NUM_THREADS':'1'})
 with gzip.open(p/'final.tsv.gz','rt') as h:
  out=[r for r in csv.DictReader(h,delimiter='\t') if r['seq_id']==r['prim_align_dup_exemplar']==r['sim_dup_exemplar']]
 c=collections.Counter(r['synthetic_category'] for r in out)
 return {'records':len(rows),'joint':len(out),'X':c['X'],'Y':c['Y'],'empty_defined':False}
results=[];started=time.monotonic();trial=0
for quality in PROTOCOL['quality_models']:
 for groups in PROTOCOL['groups']:
  for probability in PROTOCOL['random_removal_probabilities']:
   for repeat in range(PROTOCOL['random_replicates_per_cell']):
    rng=random.Random(PROTOCOL['seed']+trial);rows=[]
    for g in range(groups):
     s=strings(rng);qx,qy=(rng.randrange(20,41),rng.randrange(20,41)) if quality=='independent' else ((40,20) if quality=='X_higher' else (20,40))
     for label,q in [('X',qx),('Y',qy)]:rows.append(row(f'{g:05d}_{label}',f'toy_{g}',s,q,label))
    kept=[r for r in rows if rng.random()>=probability]
    before=evaluate(f'r{trial:04d}-before',rows);after=evaluate(f'r{trial:04d}-after',kept)
    assert before['joint']==groups
    results.append({'kind':'random_removal','quality':quality,'groups':groups,'probability':probability,'repeat':repeat,'before':before,'after':after})
    trial+=1
for kind in PROTOCOL['curated_cases']:
 for groups in PROTOCOL['curated_sizes']:
  rng=random.Random(PROTOCOL['seed']+trial);rows=[]
  for g in range(groups):
   s,t=strings(rng),strings(rng)
   if kind=='category_swap':rows.extend([row(f'{g:05d}_drop',f'toy_{g}',s,40,'Y'),row(f'{g:05d}_keep',f'toy_{g}',s,20,'X')])
   else:rows.extend([row(f'{g:05d}_drop',f'toy_{g}_A',s,40,'X'),row(f'{g:05d}_keep_A',f'toy_{g}_A',t,20,'X'),row(f'{g:05d}_keep_B',f'toy_{g}_B',s if kind=='total_increase' else t,30,'X')])
  before=evaluate(f'c{trial:04d}-before',rows);after=evaluate(f'c{trial:04d}-after',[r for r in rows if not r['seq_id'].endswith('_drop')])
  results.append({'kind':kind,'groups':groups,'before':before,'after':after});trial+=1
summary=[]
for quality in PROTOCOL['quality_models']:
 for groups in PROTOCOL['groups']:
  for probability in PROTOCOL['random_removal_probabilities']:
   rs=[r for r in results if r['kind']=='random_removal' and (r['quality'],r['groups'],r['probability'])==(quality,groups,probability)]
   changes=[r['after']['X']-r['before']['X'] for r in rs]
   summary.append({'quality':quality,'groups':groups,'removal_probability':probability,'trials':len(rs),'X_increased':sum(d>0 for d in changes),'X_decreased':sum(d<0 for d in changes),'X_unchanged':sum(d==0 for d in changes),'X_delta_min':min(changes),'X_delta_max':max(changes),'X_delta_mean':sum(changes)/len(changes),'threshold_gained':sum(r['before']['X']<3<=r['after']['X'] for r in rs),'threshold_lost':sum(r['after']['X']<3<=r['before']['X'] for r in rs)})
record={'protocol':PROTOCOL,'wall_seconds':time.monotonic()-started,'trials':len(results),'random_cells':summary,'curated':[r for r in results if r['kind']!='random_removal'],'binary_sha256':{n:hashlib.sha256((BIN/n).read_bytes()).hexdigest() for n in ['mark_duplicates','mark_duplicates_similarity']},'limits':['No empirical frequency estimate for the real cohort.','Category support and total joint-exemplar counts are different statistics.','The curated coordinate/sequence-discordant examples demonstrate possibility, not realistic abundance.','No mutation to selection rules or biological reanalysis was performed.']}
(ROOT/'trials.json').write_text(json.dumps(results,indent=2)+'\n');(ROOT/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'trials':len(results),'wall_seconds':record['wall_seconds'],'curated':record['curated']}))
