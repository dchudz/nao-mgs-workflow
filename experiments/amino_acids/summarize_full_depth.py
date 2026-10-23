"""Describe candidate and validated calls from completed production downstream rows."""
import argparse,collections,csv,functools,gzip,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--reference',required=True,type=Path);p.add_argument('--analysis-commit');a=p.parse_args()
parents={};ranks={}
for line in (a.reference/'results/taxonomy-nodes.dmp').open():
    f=[x.strip() for x in line.split('|')];parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2]
@functools.lru_cache(maxsize=None)
def species(value):
    if value in (None,'NA',''):return None
    taxid=int(value);seen=set()
    while taxid in parents and taxid not in seen:
        if ranks[taxid]=='species':return taxid
        seen.add(taxid);taxid=parents[taxid]
    return None
with gzip.open(a.reference/'results/total-virus-db-annotated.tsv.gz','rt') as f:
    annotations={int(r['taxid']):r for r in csv.DictReader(f,delimiter='\t')}
arms=['production','validation_control','screen_control','protein_variant']
contrasts=[('validation_control','production'),('screen_control','validation_control'),('protein_variant','screen_control'),('protein_variant','production')]
out=a.root/'analysis';out.mkdir(exist_ok=True)
summary=[];changes=[];species_rows=[];pair_changes=[]
comparison_fields=['aligner_taxid_lca','aligner_taxid_lca_combined','aligner_taxid_lca_artificial',
                   'validation_staxid_lca','validation_staxid_lca_combined','validation_staxid_lca_artificial',
                   'validation_status','prim_align_dup_exemplar','sim_dup_exemplar']
evidence_fields=comparison_fields+['prim_align_genome_id_all','prim_align_ref_start','prim_align_ref_start_rev',
    'protein_bitscore','protein_bitscore_rev','protein_subject','protein_subject_rev','protein_query_frame','protein_query_frame_rev']
disagreement_file=gzip.open(out/'read-disagreements.tsv.gz','wt')
disagreement=csv.DictWriter(disagreement_file,fieldnames=['sample','variant','baseline','seq_id','change']+
    [side+'_'+field for side in ['baseline','variant'] for field in evidence_fields],delimiter='\t');disagreement.writeheader()
evidence_file=gzip.open(out/'call-evidence.tsv.gz','wt')
evidence=csv.DictWriter(evidence_file,fieldnames=['sample','arm','species_taxid','seq_id','species_concordant_validated']+evidence_fields,delimiter='\t');evidence.writeheader()
disagreement_count=0
samples=sorted(p.parent for p in a.root.glob('*/complete.json') if json.loads(p.read_text()).get('all_four_outputs_completed'))
for directory in samples:
    sample=directory.name;callsets={};pairsets={};readsets={}
    for arm in arms:
        path=directory/(arm+'-downstream')/'output/results_downstream'/f'{sample}_validation_hits.tsv.gz'
        metrics=collections.Counter();candidate=collections.defaultdict(set);validated=collections.defaultdict(set);all_pairs=set();compact={}
        with gzip.open(path,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                sid=r['seq_id'];assert sid not in all_pairs;all_pairs.add(sid)
                compact[sid]={k:r.get(k,'NA') for k in evidence_fields}
                metrics['rows']+=1;metrics['status_'+r['validation_status']]+=1
                if r['prim_align_dup_exemplar']!=sid or r['sim_dup_exemplar']!=sid:continue
                metrics['joint_duplicate_exemplars']+=1
                t=species(r['aligner_taxid_lca'])
                if t is None:
                    metrics['exemplars_without_species_resolution']+=1;continue
                candidate[t].add(sid)
                evidence.writerow({'sample':sample,'arm':arm,'species_taxid':t,'seq_id':sid,
                    'species_concordant_validated':r['validation_status']=='aligned' and species(r['validation_staxid_lca'])==t,**compact[sid]})
                if r['validation_status']=='aligned' and species(r['validation_staxid_lca'])==t:
                    validated[t].add(sid);metrics['species_concordant_validated_exemplars']+=1
        pairsets[arm]=all_pairs;readsets[arm]=compact
        for t in sorted(set(candidate)|set(validated)):
            note=annotations.get(t,{})
            species_rows.append({'sample':sample,'arm':arm,'species_taxid':t,'name':note.get('name',''),
                'human_annotation':note.get('infection_status_human','NA'),'candidate_exemplars':len(candidate[t]),
                'species_concordant_validated_exemplars':len(validated[t]),'candidate_call':len(candidate[t])>=3,
                'validated_call':len(validated[t])>=3})
        callsets[arm]={'candidate':{t for t,x in candidate.items() if len(x)>=3},'validated':{t for t,x in validated.items() if len(x)>=3}}
        summary.append({'sample':sample,'arm':arm,**metrics,'candidate_species_calls':len(callsets[arm]['candidate']),
                        'validated_species_calls':len(callsets[arm]['validated'])})
    for variant,baseline in contrasts:
        for sid in sorted(pairsets[variant]|pairsets[baseline]):
            left=readsets[baseline].get(sid);right=readsets[variant].get(sid)
            if left is None:reason='candidate_gained'
            elif right is None:reason='candidate_lost'
            else:
                changed=[k for k in comparison_fields if left[k]!=right[k]]
                if not changed:continue
                reason=','.join(changed)
            record={'sample':sample,'variant':variant,'baseline':baseline,'seq_id':sid,'change':reason}
            for side,values in [('baseline',left),('variant',right)]:
                record.update({side+'_'+k:values[k] if values is not None else 'ABSENT' for k in evidence_fields})
            disagreement.writerow(record);disagreement_count+=1
        pair_changes.append({'sample':sample,'variant':variant,'baseline':baseline,
            'candidate_pairs_gained':len(pairsets[variant]-pairsets[baseline]),'candidate_pairs_lost':len(pairsets[baseline]-pairsets[variant]),
            'candidate_pairs_shared':len(pairsets[variant]&pairsets[baseline])})
        for tier in ['candidate','validated']:
            left,right=callsets[variant][tier],callsets[baseline][tier]
            for t in sorted(left^right):
                note=annotations.get(t,{})
                changes.append({'sample':sample,'variant':variant,'baseline':baseline,'tier':tier,
                    'change':'gained' if t in left else 'lost','species_taxid':t,'name':note.get('name',''),
                    'human_annotation':note.get('infection_status_human','NA')})
def write(name,rows,empty_header):
    fields=list(dict.fromkeys(k for r in rows for k in r)) or empty_header
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
write('sample-metrics.tsv',summary,['sample','arm'])
write('species-support.tsv',species_rows,['sample','arm','species_taxid'])
write('species-changes.tsv',changes,['sample','variant','baseline','tier','change','species_taxid','name','human_annotation'])
write('pair-changes.tsv',pair_changes,['sample','variant','baseline'])
disagreement_file.close();evidence_file.close()
receipt={'analysis_source_commit':a.analysis_commit,'read_disagreement_rows':disagreement_count,'completed_samples':[x.name for x in samples],'metrics':summary,'changes':changes,'pair_contrasts':pair_changes,
    'calling_rule':'Descriptive threshold of at least three production alignment-and-similarity duplicate exemplars per species per sample.',
    'validated_rule':'Each of at least three exemplars has its own aligned validation result whose natural LCA resolves to the same species as its viral alignment.',
    'limits':['Validation uses production merged/joined queries, not the earlier NR both-mate criterion.',
              'Production separates artificial references when computing the natural LCA; combined and artificial fields remain in the full read-level output.',
              'Not-sampled and no-alignment rows are retained and are not reclassified as false positives.',
              'Human-infection annotations are copied from the frozen reference and are not evidence of infection.',
              'Absence of a validated call under cap20 is not absence of that species from the sample.'],
    'evidence':'Every full validation_hits.tsv.gz, filtered BLAST result, query preparation and task diagnostic is retained in the private sample checkpoint.'}
(out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'completed_samples':len(samples),'changed_calls':len(changes)}))
