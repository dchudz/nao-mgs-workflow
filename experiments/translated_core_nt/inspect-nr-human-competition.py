"""Inspect the two observed lost-human-competitor cases in preserved raw results."""
import argparse, collections, csv, hashlib, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
pairs={'SEED_SRR33038061.678928042','SEED_SRR35987653.689121697'}
queries={pair+'/'+str(m) for pair in pairs for m in (1,2)}
allrows=[];summary=[]
a.output.mkdir(exist_ok=False)
for label,name in [('pooled','pooled-seed-alignments.tsv'),('standalone','standalone-alignments.tsv')]:
    scores={}
    for line in (a.evidence/(label+'-cache/query-calls.tsv')).open():
        q,score,taxa=line.rstrip('\n').split('\t')
        if q in queries:scores[q]=float(score)
    count=collections.Counter();human=collections.Counter();near=collections.Counter();maxhuman={};hsps=[]
    with (a.evidence/name).open() as f:
        for line in f:
            q=line.partition('\t')[0]
            if q not in queries:continue
            fields=line.rstrip('\n').split('\t');assert len(fields)==7
            count[q]+=1
            taxa={int(t) for t in fields[6].split(';') if t.isdigit()}
            if 9606 not in taxa:continue
            human[q]+=1;score=float(fields[2]);eligible=int(fields[4])>=25 and score>=.9*scores.get(q,float('inf'))
            near[q]+=eligible;maxhuman[q]=max(maxhuman.get(q,0),score)
            row=dict(context=label,query=q,subject=fields[1],bitscore=score,evalue=fields[3],alignment_AA=int(fields[4]),identity_pct=float(fields[5]),taxids=fields[6],eligible_near_best=eligible,best_eligible_query_score=scores.get(q,0));allrows.append(row)
    for q in sorted(queries):summary.append(dict(context=label,query=q,raw_alignment_rows=count[q],raw_human_alignment_rows=human[q],near_best_human_rows=near[q],maximum_raw_human_bitscore=maxhuman.get(q,0),best_eligible_query_score=scores.get(q,0)))
for name,rows in [('human-HSPs.tsv',allrows),('query-comparison.tsv',summary)]:
    with (a.output/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
result=dict(pairs=sorted(pairs),queries=summary,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_receipt_sha256=hashlib.sha256((a.evidence/'download-receipt.json').read_bytes()).hexdigest(),limits=['Cases selected after observing lost near-best human competition, not a prevalence estimate.','Human matches do not establish that a read originated from a human.','Raw output absence is distinguished from a retained match falling below the near-best threshold.','Only direct subject taxid9606 matches are tabulated; the preceding taxonomy analysis identified both cases using complete ancestry.','The fixed-mode repeat will test whether the automatic mode difference explains these specific cases.'])
(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))
