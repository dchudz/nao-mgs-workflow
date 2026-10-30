"""Create paired protein queries with exact-sequence deduplication and complete membership."""
import argparse,csv,gzip,json,re
from pathlib import Path
from union_pairs import pairs

def prepare(source,fasta,metadata,prefix=''):
    if not re.fullmatch(r'[A-Za-z0-9_-]*',prefix):raise ValueError('Invalid query prefix')
    unique={};total=0
    for key,raw in pairs(source):
        token,s1,q1,s2,q2=key;seqs=(s1,s2);total+=1
        if seqs in unique:unique[seqs].append(token.decode().lstrip('@'))
        else:unique[seqs]=[token.decode().lstrip('@')]
    with open(fasta,'x') as out,gzip.open(metadata,'xt') as meta,gzip.open(str(fasta)+'.members.tsv.gz','xt') as membership:
        writer=csv.writer(meta,delimiter='\t');writer.writerow(['query_pair','example_seq_id','occurrences','query_seq','query_seq_rev'])
        member=csv.writer(membership,delimiter='\t');member.writerow(['query_pair','seq_id'])
        for i,((a,b),ids) in enumerate(unique.items()):
            name=f'{prefix}p{i}';out.write(f'>{name}/1\n{a.decode()}\n>{name}/2\n{b.decode()}\n');writer.writerow([name,ids[0],len(ids),a.decode(),b.decode()])
            member.writerows((name,identifier) for identifier in ids)
    assert sum(map(len,unique.values()))==total
    return {'total_pairs':total,'exact_unique_pairs':len(unique),'deduplication':'Exact trimmed nucleotide sequence pair; qualities and IDs do not define sequence uniqueness. Every original pair ID is retained in the membership table.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('fasta');p.add_argument('metadata');p.add_argument('--prefix',default='');a=p.parse_args()
    result=prepare(a.source,a.fasta,a.metadata,a.prefix);Path(a.fasta+'.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
