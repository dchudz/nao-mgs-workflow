"""Parallel execution of the frozen adapter's unchanged per-query calculations.

Only complete query-pair ranges are distributed. The original FASTQ order,
qualities, near-best selection, geometry and production LCA functions remain.
"""
import argparse,contextlib,csv,json,multiprocessing,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
FROZEN=HERE if (HERE/'diamond_adapter.py').exists() else HERE.parent/'experiments/amino_acids'
sys.path.insert(0,str(FROZEN))
import diamond_adapter as serial

_STATE=None

def pair_key(line):
    return line.split(b'\t',1)[0].rsplit(b'/',1)[0]

def query_ranges(path,parts):
    size=Path(path).stat().st_size
    if not size:return []
    boundaries={0,size}
    with open(path,'rb') as f:
        for part in range(1,parts):
            f.seek(size*part//parts);f.readline()
            line=f.readline()
            if not line:continue
            key=pair_key(line)
            while True:
                boundary=f.tell();line=f.readline()
                if not line or pair_key(line)!=key:
                    boundaries.add(boundary);break
    ordered=sorted(boundaries)
    return list(zip(ordered,ordered[1:]))

def _fold(bounds):
    path,genomes,wanted,parents,artificial,unclassified,sample,extract=_STATE
    low,high=bounds;original_opener=serial.opener;original_reader=serial.read_hits
    @contextlib.contextmanager
    def bounded_opener(filename,mode='rt'):
        if str(filename)!=str(path):
            with original_opener(filename,mode) as f:yield f
            return
        with open(path,'rb') as f:
            f.seek(low)
            def lines():
                while f.tell()<high:
                    line=f.readline()
                    if not line:raise ValueError('Alignment file ended within a frozen range')
                    yield line.decode('utf-8')
            yield lines()
    first=last=None
    def tracked_reader(filename,mapping):
        nonlocal first,last
        for hit in original_reader(filename,mapping):
            key=serial.query_pair(hit)
            if last is not None and key<last:raise ValueError('DIAMOND evidence must be sorted by query ID')
            if first is None:first=key
            last=key;yield hit
    serial.opener=bounded_opener;serial.read_hits=tracked_reader
    try:
        if extract:value={serial.query_pair(hit) for hit in serial.read_hits(path,genomes)}
        else:
            value={}
            for query,hs in serial.grouped_hits(path,genomes,wanted):
                sequences,qualities=wanted[query]
                value[query]=serial.paired_row(query,sample,sequences,qualities,hs,parents,artificial,unclassified)
        return first,last,value
    finally:
        serial.opener=original_opener;serial.read_hits=original_reader

def reduce_hits(path,genomes,wanted,parents,artificial,unclassified,sample,workers,extract):
    global _STATE
    if workers<1:raise ValueError('workers must be positive')
    bounds=query_ranges(path,workers*4)
    result=set() if extract else {}
    if not bounds:return result
    _STATE=(str(path),genomes,wanted,parents,artificial,unclassified,sample,extract)
    previous=None
    # Fork shares immutable reference and query tables; only compact per-query
    # results return to the parent, never the full HSP stream.
    with multiprocessing.get_context('fork').Pool(min(workers,len(bounds))) as pool:
        for first,last,values in pool.imap(_fold,bounds,chunksize=1):
            if first is not None:
                if previous is not None and first<=previous:raise ValueError('DIAMOND evidence must be sorted by query ID without split pairs')
                previous=last
            if not extract and result.keys() & values.keys():raise ValueError('Query pair appeared in multiple ranges')
            result.update(values)
    _STATE=None
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ['hits','members','query-metadata','reads','sample','reference','genome-lengths','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--extract-only',action='store_true');p.add_argument('--workers',type=int,default=32);a=p.parse_args()
    ref=Path(a.reference)/'results';lengths=json.loads(Path(a.genome_lengths).read_text())
    with serial.opener(ref/'total-virus-db-annotated.tsv.gz') as f:viral={r['taxid'] for r in csv.DictReader(f,delimiter='\t')}
    with serial.opener(ref/'virus-genome-metadata-gid.tsv.gz') as f:
        genomes={r['genome_id']:(int(r['taxid'] if r['taxid'] in viral or r['species_taxid'] not in viral else r['species_taxid']),lengths[r['genome_id']]) for r in csv.DictReader(f,delimiter='\t')}
    membership=serial.members(a.members)
    with serial.opener(a.query_metadata) as f:query_sequences={r['query_pair']:(r['query_seq'],r['query_seq_rev']) for r in csv.DictReader(f,delimiter='\t')}
    representatives={};parents={};artificial=set();unclassified=set()
    if not a.extract_only:
        parents,children=serial.lca.parse_nodes_db(str(ref/'taxonomy-nodes.dmp'),81077)
        artificial=serial.lca.get_descendants({81077},children)
        unclassified=serial.lca.get_unclassified_taxids(serial.lca.parse_names_db(str(ref/'taxonomy-names.dmp')))
        for key,raw in serial.pairs(a.reads):
            sid,s1,q1,s2,q2=[x.decode() for x in key];sid=sid.removeprefix('@');query=membership[sid]
            if (s1,s2)!=query_sequences[query]:raise ValueError('Original read sequence differs from the searched query pair')
            representatives.setdefault(query,((s1,s2),(q1,q2)))
    reduced=reduce_hits(a.hits,genomes,representatives,parents,artificial,unclassified,a.sample,a.workers,a.extract_only)
    eligible=reduced if a.extract_only else reduced.keys();count=0;seen=set()
    with serial.opener(a.output,'wb' if a.extract_only else 'wt') as out:
        if not a.extract_only:writer=csv.DictWriter(out,fieldnames=serial.FIELDS,delimiter='\t');writer.writeheader()
        for key,raw in serial.pairs(a.reads):
            sid,s1,q1,s2,q2=[x.decode() for x in key];sid=sid.removeprefix('@')
            if sid in seen:raise ValueError('Repeated original read ID')
            seen.add(sid);query=membership[sid]
            if (s1,s2)!=query_sequences[query]:raise ValueError('Original read sequence differs from the searched query pair')
            if query not in eligible:continue
            if a.extract_only:out.write(raw)
            else:
                row=reduced[query].copy();row.update(seq_id=sid,query_qual=q1,query_qual_rev=q2)
                if len(s1)!=len(q1) or len(s2)!=len(q2):raise ValueError('Sequence/quality length mismatch')
                writer.writerow(row)
            count+=1
    Path(a.output+'.json').write_text(json.dumps({'pairs':count,'extract_only':a.extract_only,'workers':a.workers,
        'frozen_functions':str(FROZEN/'diamond_adapter.py'),'execution':'Complete query-pair byte ranges; original geometry, near-best and LCA functions; original FASTQ order and qualities.'},indent=2)+'\n')

if __name__=='__main__':main()
