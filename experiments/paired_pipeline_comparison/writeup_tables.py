"""Build reproducible, stratified report tables from accepted frozen outputs."""
import argparse,collections,csv,datetime as dt,hashlib,json,re,statistics,subprocess
from pathlib import Path

def read(p):return json.loads(p.read_text())
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def duration(value):
    scale={'ms':.001,'s':1,'m':60,'h':3600,'d':86400}
    parts=re.findall(r'([\d.]+)\s*(ms|[smhd])',value)
    assert parts and re.sub(r'([\d.]+)\s*(ms|[smhd])','',value).strip()==''
    return sum(float(n)*scale[u] for n,u in parts)

def main():
    p=argparse.ArgumentParser()
    for key in ('study','audit','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--trace-root',type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    inventory=rows(a.study/'results/aggregate/runs.tsv');done=[r for r in inventory if r['status']=='complete']
    assert all((a.audit/r['label']/'complete.json').exists() for r in done)
    input_hashes={str(a.study/'results/aggregate/runs.tsv'):sha(a.study/'results/aggregate/runs.tsv')};tables=collections.defaultdict(list)
    for depth in (5,25):
        cases=[r for r in done if int(r['depth_million_pairs'])==depth]
        audits={r['label']:read(a.audit/r['label']/'complete.json') for r in cases}
        screens=collections.Counter();species=[]
        for r in cases:
            label=r['label'];d=a.study/'results'/label;proof=audits[label]
            input_hashes[str(a.audit/label/'complete.json')]=sha(a.audit/label/'complete.json')
            input_hashes[str(a.audit/label/'species.tsv')]=sha(a.audit/label/'species.tsv')
            assert sha(a.audit/label/'species.tsv')==proof['outputs']['species.tsv']
            species.extend(rows(a.audit/label/'species.tsv'))
            screen=read(d/'union-counts.json');screens.update(screen)
            input_hashes[str(d/'union-counts.json')]=sha(d/'union-counts.json')
            assert proof['fastp']['production']['before_pairs']==screen['nucleotide_pairs']
            assert proof['fastp']['protein_variant']['before_pairs']==screen['union_pairs']
            for arm,key in [('production','production_rows'),('protein_variant','protein_rows')]:
                assert proof['arms'][arm]['rows']==int(r[key])
        tables['screen-funnel'].append(dict(depth_million_pairs=depth,units=len(cases),**screens,
            nucleotide_only=screens['nucleotide_pairs']-screens['shared_pairs']))
        for stage in ('input_pairs','screened_pairs','after_fastp_pairs','output_rows','joint_exemplars','species_resolved_exemplars','threshold_3_calls'):
            line=dict(depth_million_pairs=depth,units=len(cases),stage=stage)
            for arm in ('production','protein_variant'):
                if stage=='input_pairs':value=len(cases)*depth*1000000
                elif stage=='screened_pairs':value=screens['nucleotide_pairs' if arm=='production' else 'union_pairs']
                elif stage=='after_fastp_pairs':value=sum(x['fastp'][arm]['after_pairs'] for x in audits.values())
                elif stage=='output_rows':value=sum(x['arms'][arm]['rows'] for x in audits.values())
                elif stage=='joint_exemplars':value=sum(x['arms'][arm]['joint_exemplars'] for x in audits.values())
                elif stage=='species_resolved_exemplars':value=sum(sum(x['arms'][arm]['species_exemplars'].values()) for x in audits.values())
                else:value=sum(int(x[arm+'_exemplars'])>=3 for x in species)
                line[arm]=value
            tables['funnel'].append(line)
        for arm in ('production','protein_variant'):
            resolved=sum(sum(x['arms'][arm]['species_exemplars'].values()) for x in audits.values())
            joint=sum(x['arms'][arm]['joint_exemplars'] for x in audits.values())
            ranks=collections.Counter()
            for x in audits.values():ranks.update(x['arms'][arm]['unresolved_exemplars_by_rank'])
            assert resolved+sum(ranks.values())==joint
            tables['taxonomic-resolution'].append(dict(depth_million_pairs=depth,units=len(cases),arm=arm,
                joint_exemplars=joint,species_resolved=resolved,above_species_or_unresolved=sum(ranks.values()),
                species_resolved_fraction=resolved/joint if joint else None))
            for rank,count in sorted(ranks.items()):
                tables['unresolved-ranks'].append(dict(depth_million_pairs=depth,units=len(cases),arm=arm,rank=rank,joint_exemplars=count))
        volume=dict(depth_million_pairs=depth,units=len(cases))
        for field in ('production_rows','protein_rows','shared','same_assignment','reassigned','production_only','protein_only','net_rows','production_joint_exemplars','protein_variant_joint_exemplars'):
            volume[field]=sum(int(r[field]) for r in cases)
        volume['reassigned_fraction_of_shared']=volume['reassigned']/volume['shared'] if volume['shared'] else None
        assert volume['shared']+volume['production_only']==volume['production_rows']
        assert volume['shared']+volume['protein_only']==volume['protein_rows']
        assert volume['same_assignment']+volume['reassigned']==volume['shared']
        tables['output-volume'].append(volume)
        for n in (2,3,5):
            states=collections.Counter(x['threshold_'+str(n)] for x in species)
            pc=sum(int(x['production_exemplars'])>=n for x in species)
            qc=sum(int(x['protein_variant_exemplars'])>=n for x in species)
            assert pc==states['retained']+states['lost']
            assert qc==states['retained']+states['gained_from_zero']+states['promoted']
            tables['threshold-calls'].append(dict(depth_million_pairs=depth,units=len(cases),threshold=n,
                production_calls=pc,protein_variant_calls=qc,net_calls=qc-pc,
                **{k:states[k] for k in ('retained','gained_from_zero','promoted','lost','below_threshold_both')}))
        for n in (2,3,5):
            per_run=[]
            for case in cases:
                label=case['label'];group=[x for x in species if x['group']==label]
                pc=sum(int(x['production_exemplars'])>=n for x in group)
                qc=sum(int(x['protein_variant_exemplars'])>=n for x in group)
                per_run.append(qc-pc)
                tables['calls-per-run'].append(dict(depth_million_pairs=depth,label=label,threshold=n,
                    production_calls=pc,protein_variant_calls=qc,net_calls=qc-pc))
            tables['call-distribution'].append(dict(depth_million_pairs=depth,units=len(cases),threshold=n,
                median_net_calls=statistics.median(per_run),minimum_net_calls=min(per_run),maximum_net_calls=max(per_run),
                runs_with_gain=sum(x>0 for x in per_run),runs_with_no_net_change=sum(x==0 for x in per_run),runs_with_loss=sum(x<0 for x in per_run)))
        for arm in ('production','protein_variant'):
            for category in ('all_called','retained','gained_from_zero','promoted','lost'):
                support=[int(r[arm+'_exemplars']) for r in species if int(r[arm+'_exemplars'])>=3 and (category=='all_called' or r['threshold_3']==category)]
                if not support:continue
                tables['support-per-call'].append(dict(depth_million_pairs=depth,units=len(cases),arm=arm,threshold=3,category=category,
                    calls=len(support),joint_exemplars_in_called_species=sum(support),median_exemplars=statistics.median(support),
                    mean_exemplars=statistics.mean(support),minimum_exemplars=min(support),maximum_exemplars=max(support)))
        for r in species:tables['species-transitions'].append(dict(depth_million_pairs=depth,**r))
        for arm,side in [('production','reference'),('protein_variant','candidate')]:
            stats=collections.Counter();n_dist=agree=0;sum_dist=0
            for r in cases:
                d=a.study/'results'/r['label']/'comparison-report'
                stats.update(read(d/'validation-status-counts.json')[arm]['validation_status'])
                for path in (d/'validation-status-counts.json',d/'standard/viral_validation_agreement.tsv'):input_hashes[str(path)]=sha(path)
                standard=rows(d/'standard/viral_validation_agreement.tsv');assert len(standard)==1
                x=standard[0];n=int(x['n_validated_'+side]);n_dist+=n
                if n:
                    raw=float(x['agreement_rate_'+side])*n;assert abs(raw-round(raw))<1e-6
                    agree+=round(raw);sum_dist+=float(x['mean_distance_'+side])*n
            attempted=stats['aligned']+stats['no_alignment']
            tables['validation-audit'].append(dict(depth_million_pairs=depth,units=len(cases),arm=arm,
                rows=sum(stats.values()),aligned=stats['aligned'],no_alignment=stats['no_alignment'],
                not_sampled=stats['not_sampled'],other_status_rows=sum(v for k,v in stats.items() if k not in ('aligned','no_alignment','not_sampled')),
                attempted=attempted,no_alignment_rate=stats['no_alignment']/attempted if attempted else None,
                rows_with_taxonomic_distance=n_dist,agreement_rows=agree,
                agreement_rate=agree/n_dist if n_dist else None,mean_taxonomic_distance=sum_dist/n_dist if n_dist else None))
        for stage in ('baseline-run','aa-screen','union','prepare','protein-run','production-downstream','protein_variant-downstream'):
            wall=[];alloc=0;task_wall=0;task_count=0;derived=0;traced_units=0;missing_traces=[]
            for r in cases:
                d=a.study/'results'/r['label'];x=read(d/(stage+'.json'))
                input_hashes[str(d/(stage+'.json'))]=sha(d/(stage+'.json'))
                seconds=(dt.datetime.fromisoformat(x['finished_utc'])-dt.datetime.fromisoformat(x['started_utc'])).total_seconds()
                assert abs(seconds-x['wall_seconds'])<2
                wall.append(seconds)
                trace=d/stage/'trace.tsv'
                if not trace.exists() and a.trace_root:trace=a.trace_root/r['label']/stage/'trace.tsv'
                if not trace.exists():
                    missing_traces.append(r['label']);continue
                traced_units+=1
                input_hashes[str(trace)]=sha(trace)
                for task in rows(trace):
                    assert task['status']=='COMPLETED' and task['exit'] in ('0','-')
                    realtime=duration(task['realtime'])
                    if 'cpus' in task:cpus=int(task['cpus'])
                    else:
                        name=task['name'].split(' (')[0]
                        cpus={'FASTP':8,'PROTEIN_QUERY':1,'PROTEIN_ALIGN':32,'PROTEIN_CANDIDATES':32,'HUMAN':8,'OTHER':8,'PROTEIN_ADAPTER':32}[name]
                        derived+=1
                    alloc+=realtime*cpus/3600;task_count+=1
                    if 'complete' in task and 'start' in task:
                        task_wall+=(dt.datetime.fromisoformat(task['complete'])-dt.datetime.fromisoformat(task['start'])).total_seconds()
            tables['runtime'].append(dict(depth_million_pairs=depth,units=len(cases),stage=stage,
                workflow_wall_minutes=sum(wall)/60,median_workflow_minutes=statistics.median(wall)/60,
                min_workflow_minutes=min(wall)/60,max_workflow_minutes=max(wall)/60,
                task_allocation_hours=alloc if task_count else None,tasks=task_count,units_with_task_trace=traced_units,units_missing_task_trace=';'.join(missing_traces),
                tasks_with_cpu_allocation_from_frozen_config=derived,
                task_wall_hours_complete_minus_start=task_wall/3600 if task_count and not derived else None))
    for name,data in tables.items():
        with (a.out/(name+'.tsv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    proof=dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),accepted_comparisons=len(done),
        missing=[r['label'] for r in inventory if r['status']!='complete'],
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.study,text=True).strip(),
        analysis_sha256=sha(Path(__file__)),input_hashes=input_hashes,
        outputs={name+'.tsv':sha(a.out/(name+'.tsv')) for name in tables},
        definitions={'species':'Reference taxid_species of the aligner LCA; unresolved above-species evidence remains separate.',
        'threshold':'Added reporting rule at 2, 3, 5 joint exemplars, applied equally to both arms.',
        'validation':'No-alignment denominator is aligned + no_alignment; agreement and mean distance use rows with numeric aligner distance, pooled with those row counts as weights.',
        'runtime':'Workflow elapsed is finished-start; allocation-hours=sum(realtime*allocated cpus/3600). Rounded Nextflow realtime values limit precision. Variant prepare/protein cpus are resolved from frozen process declarations and runtime.config because those trace schemas omit cpus.',
        'scope':'Nested depth strata remain separate. New NR audit does not modify these pipeline outputs.'})
    (a.out/'provenance.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({'accepted_comparisons':len(done),'tables':list(tables),'missing':proof['missing']}))

if __name__=='__main__':main()
