"""Analyze preserved full-depth anchor evidence; never launch a sequence search.

Raw and construct-excluded policies are separate. Exclusion removes taxid 81077
and descendants only from the original eligible near-best set; it does not
retrieve lower-scoring hits. Each depth/policy has a separately formed joint
native/protein similarity pool, with the original real quality tie-breaker.
"""
import argparse,collections,csv,functools,gzip,hashlib,json,re,subprocess,sys,time
from pathlib import Path

def table(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
        yield from csv.DictReader(f,delimiter='\t')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def write_table(p,rows):
    with (gzip.open(p,'wt') if str(p).endswith('.gz') else open(p,'w')) as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def read_pairs(p):
    with gzip.open(p,'rt') as f:
        while True:
            a=[f.readline().rstrip('\n') for _ in range(4)]
            if not a[0]:break
            b=[f.readline().rstrip('\n') for _ in range(4)]
            assert all(a+b) and a[0].startswith('@') and b[0].startswith('@')
            name=lambda h:re.sub(r'/[12]$','',h.split()[0].lstrip('@'))
            assert name(a[0])==name(b[0])
            assert len(a[1])==len(a[3]) and len(b[1])==len(b[3])
            yield name(a[0]),a[1],a[3],b[1],b[3]
def ordinal(s):
    m=re.fullmatch(r'SRR39694979\.(\d+)',s)
    assert m,('Unrecognized original read ID',s)
    return int(m[1])
def main(a):
    r=a.root;o=a.output;o.mkdir(exist_ok=False);start=time.monotonic()
    parents={};ranks={}
    for line in (r/'protein/current-nr/nodes.dmp').open():
        f=[x.strip() for x in line.split('|')];parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2]
    @functools.lru_cache(None)
    def ancestors(t):
        if t not in parents:return ()
        out=[];seen=set()
        while t not in seen:
            if t not in parents:return ()
            out.append(t);seen.add(t)
            if t==parents[t]:break
            t=parents[t]
        return tuple(out)
    def species(t):return next((v for v in ancestors(t) if ranks[v]=='species'),None)
    def lca(taxa):
        paths=[ancestors(t) for t in sorted(taxa)]
        if not paths or any(not p for p in paths):return 0
        shared=set(paths[0]).intersection(*paths[1:]);return next((t for t in paths[0] if t in shared),0)
    def classify(taxa,mates):
        unknown={t for t in taxa if not ancestors(t)};known=taxa-unknown
        viral={t for t in known if 10239 in ancestors(t)};nonviral=known-viral
        status='no_alignment' if not mates else 'mixed_viral_nonviral' if viral and nonviral else 'unknown_taxonomy' if unknown else 'nonviral_only' if nonviral else 'viral_only'
        return dict(status=status,lca=0 if unknown else lca(known),mates=mates,viral=len(viral),nonviral=len(nonviral),unknown=len(unknown))
    viral={x['query_pair']:x for x in table(r/'full/protein/viral-pairs.tsv')};assert len(viral)==22359
    old={x['query_pair']:x for x in table(r/'full/protein/native-comparison.tsv')};assert old.keys()==viral.keys()
    saved={x['query_pair']:x for x in table(r/'protein/compacted-validation-v2/validated-pairs.tsv') if x['query_pair'] in viral};assert saved.keys()==viral.keys()
    meta={x['query_pair']:x for x in table(r/'full/protein/queries.tsv.gz') if x['query_pair'] in viral};assert meta.keys()==viral.keys()
    members=collections.defaultdict(set);by_id={}
    for x in table(r/'full/protein/queries.fna.members.tsv.gz'):
        q=x['query_pair']
        if q in viral:
            name=x['seq_id'];ordinal(name);assert name not in by_id;by_id[name]=q;members[q].add(name)
    assert all(len(members[q])==int(meta[q]['occurrences']) for q in viral)
    # The full query cache has already been checksum-verified before this script.
    calls={}
    with (r/'pooled-query-calls.tsv').open() as f:
        for line in f:
            if not line.startswith('FULL_'):continue
            q,score,taxids=line.rstrip('\n').split('\t')
            if q.rsplit('/',1)[0] in viral:
                assert q not in calls;calls[q]={int(t) for t in taxids.split(';')}
    evidence={};pair_rows=[]
    for q,v in viral.items():
        policies={}
        for policy in ['raw','exclude_constructs']:
            mt=[calls[q+'/'+str(m)] for m in [1,2] if q+'/'+str(m) in calls]
            if policy=='exclude_constructs':mt=[{t for t in ts if 81077 not in ancestors(t)} for ts in mt];mt=[ts for ts in mt if ts]
            result=classify(set().union(*mt) if mt else set(),len(mt));result['species']=species(result['lca']);result['concordant']=result['species'] is not None and result['species']==species(int(v['protein_taxid_lca']));policies[policy]=result
        raw=policies['raw'];s=saved[q]
        assert (raw['status'],raw['lca'],raw['mates'],raw['viral'],raw['nonviral'],raw['unknown'])==(s['nr_status'],int(s['nr_taxid_lca']),int(s['aligned_mates']),int(s['viral_taxa']),int(s['nonviral_taxa']),int(s['unknown_taxa'])),q
        evidence[q]=policies
        row=dict(query_pair=q,original_pairs=len(members[q]),first5m_original_pairs=sum(ordinal(n)<=5000000 for n in members[q]),full_exact_pair_also_native=old[q]['exact_trimmed_pair_also_native'],viral_reference_taxid=v['protein_taxid_lca'])
        row.update({policy+'_'+k:val for policy,z in policies.items() for k,val in z.items()});pair_rows.append(row)
    write_table(o/'pair-evidence.tsv',pair_rows)
    # Recover real qualities independently at each read depth, checking all memberships.
    depth_limits={'full':247633776,'first5m':5000000};quality={d:{} for d in depth_limits};seen=set()
    source=r/'full/protein/contaminants/OTHER/FULL_SRR39694979_other_bowtie2_unmapped.fastq.gz'
    for name,s1,q1,s2,q2 in read_pairs(source):
        if name not in by_id:continue
        assert name not in seen;seen.add(name);q=by_id[name];assert (s1,s2)==(meta[q]['query_seq'],meta[q]['query_seq_rev'])
        candidate=(sum(map(ord,q1+q2)),q1,q2,'@'+name)
        for d,limit in depth_limits.items():
            if ordinal(name)<=limit and candidate>quality[d].get(q,(-1,)) :quality[d][q]=candidate
    assert seen==by_id.keys();assert quality['full'].keys()==viral.keys()
    annotations={int(x['taxid']):x for x in table(r/'reference/total-virus-db-annotated.tsv.gz') if x['rank']=='species'}
    names={}
    for line in (r/'protein/current-nr/names.dmp').open():
        f=[x.strip() for x in line.split('|')]
        if f[3]=='scientific name' and (int(f[0]) in annotations or ranks.get(int(f[0]))=='species'):names[int(f[0])]=f[1]
    call_rows=[];metrics=[];all_counts={}
    for arm in a.arms:
        paths=list((r/'full/comparison/output/marked').glob(f'*_SRR39694979_{arm}_duplicate_reads.tsv.gz'));assert len(paths)==1
        native={x['seq_id']:x for x in table(paths[0])};native_quality={d:{} for d in depth_limits};native_species={d:collections.defaultdict(set) for d in depth_limits};native_ids=set()
        for name,s1,q1,s2,q2 in read_pairs(r/f'full/comparison/output/fastp/SRR39694979_{arm}_fastp.fastq.gz'):
            if name not in native:continue
            assert name not in native_ids;native_ids.add(name);seq=(s1,s2);candidate=(sum(map(ord,q1+q2)),q1,q2,'@'+name);sp=species(int(native[name]['aligner_taxid_lca']))
            for d,limit in depth_limits.items():
                if ordinal(name)>limit:continue
                if candidate>native_quality[d].get(seq,(-1,)):native_quality[d][seq]=candidate
                if sp:native_species[d][seq].add(sp)
        assert native_ids==native.keys()
        if arm=='nt':
            assert len(native_quality['full'])==11759
            assert all(((meta[q]['query_seq'],meta[q]['query_seq_rev']) in native_quality['full'])==(old[q]['exact_trimmed_pair_also_native']=='True') for q in viral)
        for depth in depth_limits:
            qs=set(quality[depth]);keys={q:(meta[q]['query_seq'],meta[q]['query_seq_rev']) for q in qs};assert len(set(keys.values()))==len(qs)
            absent={q for q in qs if keys[q] not in native_quality[depth]}
            for policy in ['raw','exclude_constructs']:
                eligible={q for q in qs if evidence[q][policy]['status']=='viral_only'}
                both={q for q in eligible if evidence[q][policy]['mates']==2};concordant={q for q in both if evidence[q][policy]['concordant']}
                counters=dict(collections.Counter(evidence[q][policy]['status'] for q in qs))
                counts=dict(candidate_pairs=len(qs),protein_only_pairs=len(absent),statuses=counters,viral_only=len(eligible),protein_only_viral_only=len(eligible&absent),protein_only_viral_only_both_mates=len(both&absent),protein_only_concordant_both_mates=len(concordant&absent))
                all_counts[f'{arm}/{depth}/{policy}']=counts
                pool=dict(native_quality[depth])
                for q in eligible:
                    if quality[depth][q]>pool.get(keys[q],(-1,)):pool[keys[q]]=quality[depth][q]
                ids={s:'p'+hashlib.sha256((s[0]+'\n'+s[1]).encode()).hexdigest() for s in pool};assert len(set(ids.values()))==len(ids)
                dest=o/f'{arm}-{depth}-{policy}';dest.mkdir()
                rows=[dict(seq_id=ids[s],query_seq=s[0],query_seq_rev=s[1],query_qual=pool[s][1],query_qual_rev=pool[s][2],prim_align_dup_exemplar=ids[s]) for s in sorted(pool,key=ids.get)]
                write_table(dest/'input.tsv.gz',rows)
                command=[str(a.binary.resolve()),'-i',str(dest/'input.tsv.gz'),'-o',str(dest/'output.tsv.gz')]
                with (dest/'stderr.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=log,timeout=600)
                cluster={x['seq_id']:x['sim_dup_exemplar'] for x in table(dest/'output.tsv.gz')};assert cluster.keys()==set(ids.values()) and set(cluster.values())<=cluster.keys()
                base=collections.defaultdict(set);prot=collections.defaultdict(set);extra=collections.defaultdict(set)
                for s,taxa in native_species[depth].items():
                    for t in taxa:base[t].add(cluster[ids[s]])
                for q in concordant:
                    t=evidence[q][policy]['species'];prot[t].add(cluster[ids[keys[q]]])
                    if q in absent:extra[t].add(cluster[ids[keys[q]]])
                for t in sorted(base.keys()|prot.keys()):
                    b,p=base[t],prot[t]
                    call_rows.append(dict(arm=arm,depth=depth,policy=policy,species=t,name=names.get(t,''),human_annotation=annotations.get(t,{}).get('infection_status_human','NA'),native_clusters=len(b),protein_clusters=len(p),protein_only_clusters=len(extra[t]),union_clusters=len(b|p),native_call=len(b)>=3,protein_call=len(p)>=3,gained_protein_call=len(p)>=3 and len(b)<3,union_call=len(b|p)>=3,gained_union_call=len(b|p)>=3 and len(b)<3))
                mapping=[dict(query_pair=q,canonical_pair=ids[keys[q]],similarity_exemplar=cluster[ids[keys[q]]],eligible_both_concordant=q in concordant,native_absent=q in absent) for q in sorted(eligible)]
                if mapping:write_table(dest/'protein-clusters.tsv.gz',mapping)
                metrics.append(dict(arm=arm,depth=depth,policy=policy,pool_exact_pairs=len(pool),pool_clusters=len(set(cluster.values())),command=command))
                print(json.dumps(dict(arm=arm,depth=depth,policy=policy,**counts)),flush=True)
    write_table(o/'species-calls.tsv',call_rows)
    result=dict(sample='SRR39694979',depth_limits=depth_limits,raw_NR_parity_pairs=len(saved),counts=all_counts,clusters=metrics,seconds=time.monotonic()-start,source_sha256=sha(__file__),binary_sha256=sha(a.binary),inputs={str(p.relative_to(r)):sha(p) for p in [r/'full/protein/viral-pairs.tsv',r/'full/protein/queries.tsv.gz',r/'full/protein/queries.fna.members.tsv.gz',r/'protein/compacted-validation-v2/validated-pairs.tsv',r/'pooled-query-calls.tsv',r/'protein/current-nr/nodes.dmp',r/'reference/total-virus-db-annotated.tsv.gz']},limits=['No new searches. Original pooled double-indexed search evidence, with raw classification parity required for all anchor pairs.','Construct exclusion removes taxid 81077 descendants within the original near-best window; does not recover lower-scoring hits.','Both mates means two NR-aligned mates whose combined eligible taxa are viral only.','Species concordance uses the same frozen taxonomy for NR and viral reference.','Clusters use the original pinned similarity algorithm, sequence-hash order, real original qualities; no alignment-coordinate duplicate stage.','Each depth/policy reclusters its joint baseline plus all viral-only candidates before both-mate/species restrictions; verdicts never transferred between reads.','First5m is a post-hoc subset by original SRA spot number, not a standalone shallow search. Anchor was not a heldout cohort member.','Native baseline not symmetrically NR validated. Human annotation is not evidence of infection.'])
    (o/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--arms',nargs='+',default=['nt','nt20']);main(p.parse_args())
