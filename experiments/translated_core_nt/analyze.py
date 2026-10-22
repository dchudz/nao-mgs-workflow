"""Describe measured same-reference alignments, retaining every queried mate."""
import argparse, collections, csv, hashlib, importlib.util, itertools, json, subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def write(path, rows, fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def main(d, methods=None, output_name="analysis"):
    methods = methods or ["megablast", "blastn", "tblastx"]
    assert len(methods) == len(set(methods)) and len(methods) >= 2
    manifest=json.loads((d/'input-manifest.json').read_text());assert m.sha(d/'queries.fna')==manifest['query_sha256']
    queries=m.fasta(d/'queries.fna');selection={x['seq_id']:x for x in csv.DictReader((d/'selection.tsv').open(),delimiter='\t')}
    allrows=[];resources={};details=[]
    for method in methods:assert (d/method/'complete.json').exists(),'Incomplete method cannot be treated as no hit: '+method
    parents={}
    for line in (m.ROOT/'reference/results/taxonomy-nodes.dmp').open():
        x=line.split('|');parents[int(x[0])]=int(x[1])
    cache={}
    def lineage(t):
        if t in cache:return cache[t]
        seen=set();path=[];original=t
        while t in parents and t not in seen:
            seen.add(t);path.append(t)
            if parents[t]==t:cache[original]=path;return path
            t=parents[t]
        cache[original]=[];return []
    def status(taxa):
        if not taxa:return 'no_alignment'
        if any(not lineage(t) for t in taxa):return 'unknown_taxonomy'
        viral=any(10239 in lineage(t) for t in taxa);nonviral=any(10239 not in lineage(t) for t in taxa)
        return 'mixed_viral_nonviral' if viral and nonviral else 'viral_only' if viral else 'nonviral_only'
    out=d/output_name;out.mkdir(exist_ok=False)
    python_image=subprocess.check_output(['docker','image','inspect','public.ecr.aws/q0n1c7g8/nao-mgs-workflow/python:e7c8524f75c0f132','--format','{{index .RepoDigests 0}}'],text=True).strip()
    lca_script=m.ROOT/'modules/local/lcaTsv/resources/usr/bin/lca_tsv.py'
    for method in methods:
        raw=[]
        for values in csv.reader((d/method/'alignments.tsv').open(),delimiter='\t'):
            assert len(values)==len(m.FIELDS)
            x=dict(zip(m.FIELDS,values));assert x['qseqid'] in queries
            for field in ['bitscore','evalue','pident','qcovs','qcovhsp']:x[field]=float(x[field])
            for field in ['length','qlen','slen','qstart','qend','sstart','send','qframe','sframe']:x[field]=int(x[field])
            x['coordinate_query_coverage_pct']=100*(abs(x['qend']-x['qstart'])+1)/x['qlen']
            assert 0 < x['coordinate_query_coverage_pct'] <= 100
            raw.append(x)
        resource=json.loads((d/method/'resources.json').read_text());assert resource['status']=='complete'
        for line in (d/method/'time.txt').read_text().splitlines():
            for text,key in [('User time (seconds):','user_seconds'),('System time (seconds):','system_seconds'),('Maximum resident set size (kbytes):','max_RSS_KiB')]:
                if text in line:resource[key]=float(line.split(text,1)[1])
        resource['raw_alignment_rows']=len(raw)
        resource['reported_HSP_coverage_disagreements']=sum(abs(x['coordinate_query_coverage_pct']-x['qcovhsp'])>1 for x in raw)
        resources[method]=resource
        for threshold in [10,1e-5]:
            label='default_E10' if threshold==10 else 'strict_E1e-5';best={}
            for x in raw:
                if x['evalue']>threshold:continue
                key=(x['qseqid'],x['sseqid']);old=best.get(key)
                if old is None or (x['bitscore'],x['length'])>(old['bitscore'],old['length']):best[key]=x
            byquery=collections.defaultdict(list)
            for x in best.values():byquery[x['qseqid']].append(x)
            retained={};taxonomy=[]
            for q,values in byquery.items():
                scores=sorted({x['bitscore'] for x in values},reverse=True);minimum=max(.9*scores[0],scores[min(9,len(scores)-1)])
                retained[q]=[x for x in values if x['bitscore']>=minimum]
                for x in retained[q]:
                    details.append(dict(method=method,threshold=label,**x))
                    taxa=[int(t) for t in x['staxids'].split(';') if t.isdigit()]
                    for t in taxa:
                        if lineage(t):taxonomy.append(dict(qseqid=q,staxid=t,bitscore=x['bitscore']))
            taxonomy.sort(key=lambda x:(x['qseqid'],-x['bitscore'],x['staxid']))
            input_path=out/(method+'-'+label+'-taxonomy.tsv');write(input_path,taxonomy,['qseqid','staxid','bitscore'])
            lca_path=out/(method+'-'+label+'-pipeline-lca.tsv')
            cmd=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',python_image,'python',str(lca_script),'-i',str(input_path),'-o',str(lca_path),'-d',str(m.ROOT/'reference/results/taxonomy-nodes.dmp'),'-n',str(m.ROOT/'reference/results/taxonomy-names.dmp'),'-g','qseqid','-t','staxid','-s','bitscore','-a','81077']
            with (out/(method+'-'+label+'-lca.log')).open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
            lcas={x['qseqid']:x for x in csv.DictReader(lca_path.open(),delimiter='\t')}
            for q in queries:
                xs=retained.get(q,[]);taxa={int(t) for x in xs for t in x['staxids'].split(';') if t.isdigit()}
                if xs and not taxa:taxa={0}
                original=selection[q.rsplit('/',1)[0]];call=lcas.get(q,{})
                allrows.append(dict(query=q,pair=q.rsplit('/',1)[0],sample=original['sample'],stratum=original['stratum'],method=method,threshold=label,status=status(taxa),pipeline_natural_lca=call.get('staxid_lca','NA'),pipeline_combined_lca=call.get('staxid_lca_combined','NA'),near_best_taxids=';'.join(map(str,sorted(taxa))),retained_subjects=len(xs),best_evalue=min((x['evalue'] for x in xs),default='NA'),best_bitscore=max((x['bitscore'] for x in xs),default='NA')))
    assert len(allrows)==len(queries)*len(methods)*2
    write(out/'per-mate.tsv',allrows,list(allrows[0]));write(out/'near-best-alignments.tsv',details,['method','threshold']+m.FIELDS+['coordinate_query_coverage_pct'])
    lookup={(x['query'],x['method'],x['threshold']):x for x in allrows};disagreements=[];groups={}
    for threshold in ['default_E10','strict_E1e-5']:
        for left,right in itertools.combinations(methods,2):
            changed=[]
            for q in queries:
                a,b=lookup[q,left,threshold],lookup[q,right,threshold]
                if (a['status'],a['pipeline_natural_lca'])!=(b['status'],b['pipeline_natural_lca']):
                    changed.append(q);disagreements.append(dict(query=q,pair=a['pair'],stratum=a['stratum'],threshold=threshold,left=left,right=right,left_status=a['status'],right_status=b['status'],left_lca=a['pipeline_natural_lca'],right_lca=b['pipeline_natural_lca']))
            groups[left+' -> '+right+'/'+threshold]=dict(disagreeing_mates=len(changed),disagreeing_pairs=len({q.rsplit('/',1)[0] for q in changed}))
    write(out/'disagreements.tsv',disagreements,['query','pair','stratum','threshold','left','right','left_status','right_status','left_lca','right_lca'])
    summary=dict(pairs=len(selection),mates=len(queries),resources=resources,disagreements=groups,status_counts={method+'/'+threshold:dict(collections.Counter(x['status'] for x in allrows if x['method']==method and x['threshold']==threshold)) for method in methods for threshold in ['default_E10','strict_E1e-5']},production_lca_source_sha256=m.sha(lca_script),python_image=python_image,notes=['Counts describe a selected diagnostic panel, not biological truth or prevalence','Natural and combined LCA computed by unchanged production lca_tsv.py','No-hit mates retained in denominator; all three searches must finish','Top10 dense score ranks and90percent near-best applied separately within each mate and method','Strict E1e-5 is a postfilter of returned E10 hits, not a separate search','Unknown taxonomy included in status; unknown entries omitted from production LCA input','Alignment length/identity units are amino acids for tblastx and nucleotides for the other methods','Separate-mate comparison is not a complete production-pipeline replay','Reported tblastx qcovhsp can disagree with nucleotide coordinate coverage; retain raw and derived values. Synthetic threshold-boundary checks support internal filtering, not display, on tested cases'])
    m.save(out/'summary.json',summary);print(json.dumps(summary))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--methods',nargs='+');p.add_argument('--output-name',default='analysis');args=p.parse_args();main(args.directory,args.methods,args.output_name)
