"""Annotate measured NR context disagreements without changing original classifications."""
import argparse, collections, csv, hashlib, importlib.util, json
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('evidence',type=Path)
p.add_argument('standalone',type=Path)
p.add_argument('classifier',type=Path)
p.add_argument('output',type=Path)
a=p.parse_args()
a.output.mkdir(exist_ok=False)
parents={}
for line in (a.evidence/'nodes.dmp').open():
    f=line.split('|');parents[int(f[0])]=int(f[1])
paths={}
def lineage(t):
    if t in paths:return paths[t]
    initial=t;seen=set();result=[]
    while t in parents and t not in seen:
        seen.add(t);result.append(t)
        if parents[t]==t:paths[initial]=result;return result
        t=parents[t]
    paths[initial]=None
    return None

spec=importlib.util.spec_from_file_location('original_nr_classifier',a.classifier)
original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
old_path=a.evidence/'pooled-seed-validation.tsv'
new_path=a.standalone/'seed-validation.tsv'
truth={label:{x['query_pair']:x for x in csv.DictReader(path.open(),delimiter='\t')} for label,path in [('pooled',old_path),('standalone',new_path)]}
assert truth['pooled'].keys()==truth['standalone'].keys() and len(truth['pooled'])==746
viral_calls={q:int(x['native_reference_taxid']) for q,x in truth['pooled'].items()}
observed={};allrows=[]
for label in ('pooled','standalone'):
    calls={}
    for line in (a.evidence/(label+'-cache')/'query-calls.tsv').open():
        q,score,taxa=line.rstrip('\n').split('\t');assert q not in calls
        calls[q]=({int(t) for t in taxa.split(';')},float(score))
    # Reuse the original pair classifier over the existing verified sufficient statistic.
    original.near_best_query_taxa=lambda unused:calls
    reproduced=original.summarize([],parents,viral_calls)
    for row in reproduced:
        expected=truth[label][row['query_pair']]
        for key in ('nr_status','taxonomic_relation'):
            assert row[key]==expected[key],(label,row['query_pair'],key)
        for key in ('nr_taxid_lca','aligned_mates','viral_taxa','nonviral_taxa','unknown_taxa'):
            assert row[key]==int(expected[key]),(label,row['query_pair'],key)
        assert row['best_bitscore']==float(expected['best_bitscore'])
    observed[label]={}
    for pair in viral_calls:
        taxa=set().union(*(calls.get(pair+'/'+str(m),(set(),0))[0] for m in (1,2)))
        known={t for t in taxa if lineage(t) is not None}
        groups={name:{t for t in known if root in lineage(t)} for name,root in [('viral',10239),('artificial',81077),('bacterial',2),('archaeal',2157),('eukaryotic',2759),('human',9606)]}
        groups['cellular']=groups['bacterial']|groups['archaeal']|groups['eukaryotic']
        groups['other_known_nonviral']=known-groups['viral']-groups['artificial']-groups['cellular']
        groups['unknown']=taxa-known
        natural_nonviral=known-groups['viral']-groups['artificial']
        natural_status='unknown_taxonomy' if groups['unknown'] else 'mixed_viral_natural_nonviral' if groups['viral'] and natural_nonviral else 'viral_only' if groups['viral'] else 'natural_nonviral_only' if natural_nonviral else 'artificial_only' if groups['artificial'] else 'no_alignment'
        row=dict(query_pair=pair,context=label,original_status=truth[label][pair]['nr_status'],natural_status=natural_status)
        for name,values in groups.items():row[name+'_taxa']=len(values);row[name+'_taxids']=';'.join(map(str,sorted(values)))
        observed[label][pair]=row;allrows.append(row)

changed=[]
for pair in viral_calls:
    left=observed['pooled'][pair];right=observed['standalone'][pair]
    if left['original_status']!=right['original_status'] or left['natural_status']!=right['natural_status']:
        changed.append(dict(query_pair=pair,**{'pooled_'+k:v for k,v in left.items() if k not in ['query_pair','context']},**{'standalone_'+k:v for k,v in right.items() if k not in ['query_pair','context']}))
for name,rows in [('per-pair-competitors.tsv',allrows),('status-change-competitors.tsv',changed)]:
    with (a.output/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['none'],delimiter='\t');w.writeheader();w.writerows(rows)
mixed_to_viral=[q for q in viral_calls if truth['pooled'][q]['nr_status']=='mixed_viral_nonviral' and truth['standalone'][q]['nr_status']=='viral_only']
counts=dict(
    original_classifications_reproduced=1492,
    pairs=746,
    mixed_to_viral_pairs=len(mixed_to_viral),
    mixed_to_viral_with_prior_artificial_only_nonviral=sum(observed['pooled'][q]['artificial_taxa']>0 and observed['pooled'][q]['cellular_taxa']==0 and observed['pooled'][q]['other_known_nonviral_taxa']==0 and observed['pooled'][q]['unknown_taxa']==0 for q in mixed_to_viral),
    mixed_to_viral_with_prior_cellular=sum(observed['pooled'][q]['cellular_taxa']>0 for q in mixed_to_viral),
    mixed_to_viral_with_prior_human=sum(observed['pooled'][q]['human_taxa']>0 for q in mixed_to_viral),
    mixed_to_viral_with_prior_other_known_nonviral=sum(observed['pooled'][q]['other_known_nonviral_taxa']>0 for q in mixed_to_viral),
    natural_status_changes=sum(observed['pooled'][q]['natural_status']!=observed['standalone'][q]['natural_status'] for q in viral_calls),
    original_status_transitions=dict(collections.Counter(truth['pooled'][q]['nr_status']+' -> '+truth['standalone'][q]['nr_status'] for q in viral_calls if truth['pooled'][q]['nr_status']!=truth['standalone'][q]['nr_status'])),
    changed_natural_status_transitions=dict(collections.Counter(observed['pooled'][q]['natural_status']+' -> '+observed['standalone'][q]['natural_status'] for q in viral_calls if observed['pooled'][q]['natural_status']!=observed['standalone'][q]['natural_status'])),
)
inputs=[a.evidence/'download-receipt.json',old_path,new_path,a.evidence/'pooled-cache/query-calls.tsv',a.evidence/'standalone-cache/query-calls.tsv',a.classifier]
result=dict(counts=counts,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_hashes={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},limits=[
    'Supplementary taxonomy annotation of completed original results; no new sequence search and no changed primary classifications.',
    'The original pair classifier exactly reproduces all1492 saved context-specific classifications and their score/taxon-count fields from the compact evidence.',
    'Cellular groups use bacterial, archaeal and eukaryotic ancestry. Artificial constructs and other known nonviral taxa are separate; cellular and artificial sets need not cover every known nonviral taxon.',
    'Counts of taxa describe the returned near-best set, not all competing sequences in NR or the true source organism.',
    'These746 pairs are a selected diagnostic panel, not a prevalence sample or the complete200-sample candidate set.',
    'Both query context and automatic search mode differ. The fixed-mode repeat remains the separate causal control.',
])
(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(counts,indent=2))
