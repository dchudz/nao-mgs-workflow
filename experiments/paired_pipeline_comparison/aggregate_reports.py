"""Summarize integrity-checked paired reports without adding a calling rule."""
import argparse
import collections
import csv
import json
from pathlib import Path
import statistics


def read(path):
    return json.loads(path.read_text())


def collect(results, selection, orders):
    rows=[]
    for sample in (r['run_accession'] for r in selection['selected']):
        for depth in ([5,25] if sample in selection['depth_runs'] else [5]):
            label=f'{sample}_{depth}m';d=results/label
            row={'sample':sample,'depth_million_pairs':depth,'label':label,
                 'feasibility_run':sample in selection['pilot_runs'],
                 'validation_first':orders['orders'][label][0],
                 'status':'not_collected'}
            report=d/'comparison-report'
            if not (report/'accounting.json').exists():
                if d.exists():row['status']='incomplete'
                rows.append(row);continue
            try:
                a=read(report/'accounting.json');v=read(report/'validation-status-counts.json')
                integrity=read(report/'integrity.json');complete=read(d/'complete.json')
                assert a['completed'] and integrity['passed']
                assert complete['both_arms_complete'] and complete['downstream_parameters_equal']
                assert a['sample']==label and a['pairs_per_arm']==depth*1_000_000
                assert read(d/'synthetic-control-gate.json')['passed']
                assert read(d/'bowtie-output-integrity.json')['passed']
                with (report/'standard/viral_read_status.tsv').open() as f:
                    statuses=[x for x in csv.DictReader(f,delimiter='\t') if x['scope']=='all']
                assert len(statuses)==1 and statuses[0]['group']==label
                s={k:int(statuses[0][k]) for k in ('n_reference','n_candidate','n_shared','n_same','n_reassigned','n_lost','n_gained')}
                assert s['n_reference']==s['n_shared']+s['n_lost']==v['production']['rows']
                assert s['n_candidate']==s['n_shared']+s['n_gained']==v['protein_variant']['rows']
                assert s['n_shared']==s['n_same']+s['n_reassigned']==a['shared_output_pair_ids']
                assert s['n_lost']==a['production_only_output_pair_ids']
                assert s['n_gained']==a['protein_only_output_pair_ids']
                for arm in ('production','protein_variant'):
                    assert sum(v[arm]['validation_status'].values())==v[arm]['rows']
                    assert 0<=v[arm]['joint_exemplars']<=v[arm]['rows']
                proof=read(d/'completed-evidence-backup.json')
                assert proof['readback_verified'] and len(proof['sha256'])==64
                row.update(status='complete',source_commit=complete['source_commit'],
                           unavailable_outputs=','.join(a['shared_unavailable_outputs']),
                           stage_minutes=a['sum_stage_wall_seconds']/60,
                           production_rows=s['n_reference'],protein_rows=s['n_candidate'],
                           shared=s['n_shared'],same_assignment=s['n_same'],reassigned=s['n_reassigned'],
                           production_only=s['n_lost'],protein_only=s['n_gained'],
                           net_rows=s['n_candidate']-s['n_reference'])
                for arm in ('production','protein_variant'):
                    row[arm+'_joint_exemplars']=v[arm]['joint_exemplars']
                    for status in sorted({'aligned','no_alignment','not_sampled'}|set(v[arm]['validation_status'])):
                        row[arm+'_'+status]=v[arm]['validation_status'].get(status,0)
            except (AssertionError,KeyError,ValueError,OSError) as error:
                row.update(status='blocked',error=type(error).__name__+': '+str(error))
            rows.append(row)
    return rows


def summarize(rows):
    strata=[]
    for depth in (5,25):
        all_rows=[r for r in rows if r['depth_million_pairs']==depth]
        done=[r for r in all_rows if r['status']=='complete']
        summary={'depth_million_pairs':depth,'planned_units':len(all_rows),
                 'status_counts':dict(collections.Counter(r['status'] for r in all_rows))}
        if done:
            summary['completed_run_count']=len(done)
            for metric in ('production_rows','protein_rows','production_only','protein_only','reassigned','net_rows','stage_minutes'):
                values=[r[metric] for r in done]
                summary[metric]={'sum':sum(values),'median':statistics.median(values),'min':min(values),'max':max(values)}
        strata.append(summary)
    return {'strata':strata,'interpretation':[
        'Units are sequencing runs, not proven independent specimens.',
        'Output counts and changed assignments are not sensitivity, specificity, or confirmed new detections.',
        'No additional taxonomic agreement or exemplar threshold is imposed.',
        'Validation statuses have different sampled sets; unaligned and unsampled records are kept separate.',
        'Stage minutes exclude acquisition, boot, preservation and idle time, and share profiling/QC work.',
        'Not-collected units may be running remotely; absence from this local snapshot is not a zero result.',
        'Depth strata are nested; their rows must not be pooled as independent observations.',
        'Unavailable Bracken output is not zero abundance.']}


def main():
    p=argparse.ArgumentParser()
    for name in ('results','selection','orders','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();rows=collect(a.results,read(a.selection),read(a.orders))
    a.out.mkdir(parents=True,exist_ok=True)
    columns=list(dict.fromkeys(key for row in rows for key in row))
    with (a.out/'runs.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
    (a.out/'summary.json').write_text(json.dumps(summarize(rows),indent=2)+'\n')
    print(json.dumps({'status_counts':dict(collections.Counter(r['status'] for r in rows))}))
    if any(r['status']=='blocked' for r in rows):raise SystemExit('Blocked result; inspect runs.tsv')


if __name__=='__main__':main()
