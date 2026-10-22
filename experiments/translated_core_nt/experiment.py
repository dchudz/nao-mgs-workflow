"""Same-query, same-reference translated BLAST experiment; never edit production code."""
import argparse, csv, datetime, hashlib, json, subprocess, time
from pathlib import Path

IMAGE = 'public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
ROOT = Path('/home/ec2-user/experiment')
DB = ROOT / 'scratch/blast_db_190b4ae1/blast_db'
OUTFMT = '6 qseqid sseqid staxids bitscore evalue length pident qcovs qcovhsp qlen slen qstart qend sstart send qframe sframe qseq sseq'
FIELDS = OUTFMT.split()[1:]

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, obj): Path(p).write_text(json.dumps(obj, indent=2)+'\n')
def fasta(p):
    records = {}; name = None
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0]; assert name not in records
            records[name]=''
        else:
            assert name is not None; records[name]+=line.strip()
    assert records and all(records.values()); return records

def prepare(d, n):
    source=ROOT/'cohort-storage/seed-diagnostic200'
    rows=list(csv.DictReader((source/'selection.tsv').open(),delimiter='\t'))
    seqs=fasta(source/'queries.fna'); selected=[]; used=set()
    for flag in ['baseline_control','AA12_added','NT20_added']:
        candidates=sorted((x for x in rows if x[flag]=='True' and x['seq_id'] not in used), key=lambda x:hashlib.sha256(('same-reference-v1:'+flag+':'+x['seq_id']).encode()).digest())
        picked=[]; samples=set()
        for row in candidates:
            if row['sample'] not in samples:
                picked.append(dict(row, stratum=flag));samples.add(row['sample'])
            if len(picked)==n:break
        assert len(picked)==n
        selected.extend(picked);used.update(x['seq_id'] for x in picked)
    d.mkdir(parents=True,exist_ok=False)
    with (d/'selection.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t');w.writeheader();w.writerows(selected)
    q={k:v for row in selected for k,v in seqs.items() if k.rsplit('/',1)[0]==row['seq_id']}
    assert len(q)==sum(int(x['available_mates']) for x in selected)
    (d/'queries.fna').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in q.items()))
    dbfiles=sorted(DB.parent.glob('*'))
    inventory=[dict(name=p.name,bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,**({'sha256':sha(p)} if p.suffix in ['.nal','.njs'] else {})) for p in dbfiles if p.is_file()]
    save(d/'reference-inventory.json',inventory)
    info=subprocess.check_output(['docker','run','--rm','--network','none','-v',f'{ROOT}:{ROOT}:ro',IMAGE,'blastdbcmd','-db',str(DB),'-info'],text=True)
    (d/'reference-info.txt').write_text(info)
    save(d/'input-manifest.json',dict(created_utc=now(),pairs=len(selected),queries=len(q),per_stratum=n,query_sha256=sha(d/'queries.fna'),selection_sha256=sha(d/'selection.tsv'),reference_inventory_sha256=sha(d/'reference-inventory.json'),source_query_sha256=sha(source/'queries.fna'),source_selection_sha256=sha(source/'selection.tsv'),selection_rule='SHA256 same-reference-v1:flag:seq_id; disjoint strata baseline,AA12_added,NT20_added in that order; unique samples within each stratum; no BLAST or NR outcome used',limits=['Selected diagnostic panel, not population accuracy or prevalence','Identical separate available mate sequences; production validation merges/joins mates','Baseline and additions deliberately enriched; strata overlap in parent panel but each pair included once here']))
    print(json.dumps(json.loads((d/'input-manifest.json').read_text())),flush=True)

def run(d, config, commit, seconds):
    manifest=json.loads((d/'input-manifest.json').read_text());method=json.loads(config.read_text());name=method['name']
    assert sha(d/'queries.fna')==manifest['query_sha256'] and sha(d/'selection.tsv')==manifest['selection_sha256']
    for x in json.loads((d/'reference-inventory.json').read_text()):
        p=DB.parent/x['name'];assert p.stat().st_size==x['bytes'] and p.stat().st_mtime_ns==x['mtime_ns']
        if 'sha256' in x:assert sha(p)==x['sha256']
    output=d/name;output.mkdir(exist_ok=False)
    command=[method['executable']]+method['arguments']+['-db',str(DB),'-query',str(d/'queries.fna'),'-out',str(output/'alignments.tsv'),'-max_hsps','5','-num_alignments','250','-qcov_hsp_perc','30','-evalue','10','-num_threads','16','-outfmt',OUTFMT]
    container='same-ref-'+d.name+'-'+name
    docker=['docker','run','--rm','--name',container,'--network','none','--cpus','16','--memory','450g','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{ROOT}:{ROOT}',IMAGE,'/usr/bin/time','-v','-o',str(output/'time.txt')]+command
    cutoff=datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00')
    budget_seconds=min(seconds,(cutoff-datetime.datetime.now(datetime.timezone.utc)).total_seconds());assert budget_seconds>60
    protocol=dict(started_utc=now(),source_commit=commit,source_sha256=sha(__file__),method_config_sha256=sha(config),input_manifest=manifest,method=method,command=docker,maximum_seconds=budget_seconds,output_fields=FIELDS,allocation_USD_per_hour_conservative=4.0,cache_policy='Existing shared reference cache, not a cold-cache comparison',completion_required=True)
    save(output/'protocol.json',protocol); start=time.monotonic();status='failed';returncode=None
    try:
        with (output/'stderr.log').open('w') as log:
            result=subprocess.run(docker,stdout=log,stderr=log,timeout=budget_seconds)
        returncode=result.returncode;status='complete' if returncode==0 else 'failed'
    except subprocess.TimeoutExpired:status='time_limit'
    finally:
        subprocess.run(['docker','stop','-t','10',container],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        resource=dict(status=status,returncode=returncode,finished_utc=now(),wall_seconds=time.monotonic()-start)
        resource['allocated_cost_upper_estimate_USD']=resource['wall_seconds']/3600*4
        save(output/'resources.json',resource);print(json.dumps(dict(method=name,**resource)),flush=True)
    if status!='complete':raise SystemExit(2)
    query_ids=set(fasta(d/'queries.fna'));count=0
    with (output/'alignments.tsv').open() as f:
        for row in csv.reader(f,delimiter='\t'):
            assert len(row)==len(FIELDS) and row[0] in query_ids;count+=1
    save(output/'complete.json',dict(completed_utc=now(),rows=count,queries=len(query_ids),output_sha256=sha(output/'alignments.tsv')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('--directory',type=Path,required=True);p.add_argument('--per-stratum',type=int,default=8);p.add_argument('--config',type=Path);p.add_argument('--commit');p.add_argument('--seconds',type=int,default=21600);a=p.parse_args()
    if a.action=='prepare':prepare(a.directory,a.per_stratum)
    else:
        assert a.config and a.commit;run(a.directory,a.config,a.commit,a.seconds)
