"""Create explicitly paired, exact-sequence deduplicated protein-search queries."""
import argparse,collections,csv,gzip,json
from pathlib import Path
from union_pairs import pairs

def prepare(source,fasta,metadata):
    unique={};total=0
    for key,raw in pairs(source):
        token,s1,q1,s2,q2=key;seqs=(s1,s2);total+=1
        if seqs in unique:unique[seqs][1]+=1
        else:unique[seqs]=[token.decode(),1]
    with open(fasta,'x') as out, gzip.open(metadata,'wt') as meta:
        writer=csv.writer(meta,delimiter='\t');writer.writerow(['query_pair','example_seq_id','occurrences','query_seq','query_seq_rev'])
        for i,((a,b),(token,count)) in enumerate(unique.items()):
            name=f'p{i}';out.write(f'>{name}/1\n{a.decode()}\n>{name}/2\n{b.decode()}\n');writer.writerow([name,token.lstrip('@'),count,a.decode(),b.decode()])
    return {'total_pairs':total,'exact_unique_pairs':len(unique),'deduplication':'Exact trimmed nucleotide sequence pair; qualities and IDs do not define molecular uniqueness.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('fasta');p.add_argument('metadata');a=p.parse_args()
    result=prepare(a.source,a.fasta,a.metadata);Path(a.fasta+'.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
