"""Call summaries require independent exemplars and each read's own validation."""
import csv,gzip,json,subprocess,sys,tempfile,unittest
from pathlib import Path

class SummaryTests(unittest.TestCase):
    def test_unsampled_and_duplicate_rows_do_not_add_validated_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'reference/results';ref.mkdir(parents=True)
            (ref/'taxonomy-nodes.dmp').write_text('1\t|\t1\t|\tno rank\t|\n100\t|\t1\t|\tgenus\t|\n111\t|\t100\t|\tspecies\t|\n112\t|\t111\t|\tstrain\t|\n222\t|\t100\t|\tspecies\t|\n')
            with gzip.open(ref/'total-virus-db-annotated.tsv.gz','wt') as f:
                f.write('taxid\tname\tinfection_status_human\n111\tSynthetic A\t0\n222\tSynthetic B\t1\n')
            d=root/'cohort/sample';d.mkdir(parents=True);(d/'complete.json').write_text('{"all_four_outputs_completed":true}')
            def row(sid,t,status='aligned',exemplar=None):
                return dict(seq_id=sid,prim_align_dup_exemplar=exemplar or sid,sim_dup_exemplar=exemplar or sid,
                            aligner_taxid_lca=str(t),validation_staxid_lca=str(t) if status=='aligned' else 'NA',validation_status=status)
            baseline=[row(f'a{i}',112) for i in range(2)]+[row(f'b{i}',222) for i in range(3)]
            variant=[row(f'a{i}',112) for i in range(2)]+[row('a2',112,'not_sampled'),row('a-duplicate',112,exemplar='a0')]+[row(f'b{i}',222) for i in range(2)]+[row('unresolved',100)]
            for arm in ['production','validation_control','screen_control','protein_variant']:
                out=d/(arm+'-downstream')/'output/results_downstream';out.mkdir(parents=True)
                with gzip.open(out/'sample_validation_hits.tsv.gz','wt') as f:
                    w=csv.DictWriter(f,fieldnames=list(baseline[0]),delimiter='\t');w.writeheader();w.writerows(variant if arm=='protein_variant' else baseline)
            script=Path(__file__).with_name('summarize_full_depth.py')
            subprocess.run([sys.executable,str(script),str(root/'cohort'),'--reference',str(root/'reference')],check=True,capture_output=True)
            result=json.loads((root/'cohort/analysis/summary.json').read_text())
            changes={(r['tier'],r['change'],r['species_taxid']) for r in result['changes'] if r['variant']=='protein_variant' and r['baseline']=='production'}
            self.assertEqual(changes,{('candidate','gained',111),('candidate','lost',222),('validated','lost',222)})
            metric=next(r for r in result['metrics'] if r['arm']=='protein_variant')
            self.assertEqual(metric['species_concordant_validated_exemplars'],4)
            self.assertEqual(metric['exemplars_without_species_resolution'],1)
            with gzip.open(root/'cohort/analysis/call-evidence.tsv.gz','rt') as f:
                evidence={r['seq_id']:r for r in csv.DictReader(f,delimiter='\t') if r['arm']=='protein_variant'}
            self.assertNotIn('a-duplicate',evidence)
            self.assertEqual(evidence['a2']['species_concordant_validated'],'False')
            with gzip.open(root/'cohort/analysis/read-disagreements.tsv.gz','rt') as f:
                differences={r['seq_id']:r for r in csv.DictReader(f,delimiter='\t') if r['variant']=='protein_variant' and r['baseline']=='production'}
            self.assertEqual(differences['a2']['change'],'candidate_gained')
            self.assertEqual(differences['b2']['change'],'candidate_lost')
            self.assertEqual(differences['b2']['baseline_aligner_taxid_lca'],'222')

if __name__=='__main__':unittest.main()
