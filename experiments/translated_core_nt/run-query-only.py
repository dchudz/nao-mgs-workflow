"""Full-reference query-only translation experiment with explicit per-read normalization."""
import argparse,collections,csv,datetime,importlib.util,json,subprocess,time
from pathlib import Path
spec=importlib.util.spec_from_file_location('translation',Path(__file__).with_name('characterize-query-translation.py'));t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);m=t.m

def normalize_pool(raw,mapping,records):
    bypair=collections.defaultdict(list);raw_count=0
    for values in csv.reader(raw.open(),delimiter='\t'):
        assert len(values)==len(m.FIELDS);row=t.normalize(dict(zip(m.FIELDS,values)),mapping,records);raw_count+=1
        bypair[row['qseqid'],row['sseqid']].append(row)
    # Preserve all raw evidence separately; normalize the stated retrieval caps per original read.
    byquery=collections.defaultdict(list)
    for (query,subject),rows in bypair.items():
        rows.sort(key=lambda x:(-float(x['bitscore']),-int(x['length']),int(x['qstart']),int(x['sstart']),int(x['qframe']),int(x['sframe'])))
        byquery[query].append((subject,rows[:5]))
    normalized=[];dropped_subjects=0
    for query,subjects in sorted(byquery.items()):
        subjects.sort(key=lambda x:(-float(x[1][0]['bitscore']),-int(x[1][0]['length']),x[0]));dropped_subjects+=max(0,len(subjects)-250)
        for subject,rows in subjects[:250]:
            spans=sorted((min(int(x['qstart']),int(x['qend'])),max(int(x['qstart']),int(x['qend']))) for x in rows);merged=[]
            for left,right in spans:
                if merged and left<=merged[-1][1]+1:merged[-1][1]=max(merged[-1][1],right)
                else:merged.append([left,right])
            coverage=100*sum(right-left+1 for left,right in merged)/len(records[query])
            for row in rows:row['qcovs']=str(coverage);normalized.append(row)
    return normalized,dict(raw_frame_HSPs=raw_count,normalized_HSPs=len(normalized),dropped_subjects_after_pool_cap=dropped_subjects,original_queries=len(records),aligned_original_queries=len(byquery))

def main(d,commit,seconds):
    cfg=Path(__file__).with_name('query-tblastn-fast.json');method=json.loads(cfg.read_text());manifest=json.loads((d/'input-manifest.json').read_text())
    assert m.sha(d/'queries.fna')==manifest['query_sha256'] and m.sha(d/'selection.tsv')==manifest['selection_sha256']
    for x in json.loads((d/'reference-inventory.json').read_text()):
        p=m.DB.parent/x['name'];assert p.stat().st_size==x['bytes'] and p.stat().st_mtime_ns==x['mtime_ns']
        if 'sha256' in x:assert m.sha(p)==x['sha256']
    out=d/method['name'];out.mkdir(exist_ok=False);start=time.monotonic();started=m.now();records,mapping=t.translate_queries(d/'queries.fna',out/'queries.faa');assert len(mapping)==6*len(records);m.save(out/'query-map.json',mapping)
    container='same-ref-'+d.name+'-'+method['name']
    command=[method['executable']]+method['arguments']+['-db',str(m.DB),'-query',str(out/'queries.faa'),'-out',str(out/'frame-alignments.tsv'),'-max_hsps','5','-num_alignments','250','-qcov_hsp_perc','30','-evalue','10','-num_threads','16','-outfmt',m.OUTFMT]
    docker=['docker','run','--rm','--name',container,'--network','none','--cpus','16','--memory','450g','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE,'/usr/bin/time','-v','-o',str(out/'time.txt')]+command
    cutoff=datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00');limit=min(seconds,(cutoff-datetime.datetime.now(datetime.timezone.utc)).total_seconds());assert limit>60
    m.save(out/'protocol.json',dict(started_utc=started,source_commit=commit,source_sha256=m.sha(__file__),translation_source_sha256=m.sha(Path(__file__).with_name('characterize-query-translation.py')),method=method,input_manifest=manifest,translated_query_sha256=m.sha(out/'queries.faa'),translated_queries=len(mapping),command=docker,maximum_seconds=limit,allocation_USD_per_hour_conservative=4.0,output_fields=m.FIELDS,normalization='Map every frame to original nucleotide coordinates; pool per original read; retain250best subjects and5best HSPs per subject; aggregate nucleotide qcovs over retained HSPs. Raw frame outputs and warnings retained. No verdict extrapolated between reads.'))
    status='failed';returncode=None
    try:
        with (out/'stderr.log').open('w') as log:result=subprocess.run(docker,stdout=log,stderr=log,timeout=limit)
        returncode=result.returncode;status='complete' if returncode==0 else 'failed'
    except subprocess.TimeoutExpired:status='time_limit'
    finally:
        subprocess.run(['docker','stop','-t','10',container],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        elapsed=time.monotonic()-start;m.save(out/'resources.json',dict(status=status,returncode=returncode,finished_utc=m.now(),wall_seconds=elapsed,allocated_cost_upper_estimate_USD=elapsed/3600*4))
    if status!='complete':raise SystemExit(2)
    rows,counts=normalize_pool(out/'frame-alignments.tsv',mapping,records)
    with (out/'alignments.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=m.FIELDS,delimiter='\t');w.writerows(rows)
    m.save(out/'normalization.json',counts);m.save(out/'complete.json',dict(completed_utc=m.now(),rows=len(rows),queries=len(records),output_sha256=m.sha(out/'alignments.tsv'),raw_output_sha256=m.sha(out/'frame-alignments.tsv')));print(json.dumps(counts))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--commit',required=True);p.add_argument('--seconds',type=int,default=21600);a=p.parse_args();main(a.directory,a.commit,a.seconds)
