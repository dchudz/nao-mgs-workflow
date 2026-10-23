"""Production duplicate fields retain original mates and measured HSP geometry."""
import dataclasses
import unittest
import tempfile
from pathlib import Path
from diamond_adapter import Hit, paired_row, group_from_hits, grouped_hits, near_best, lca
from diamond_geometry import Hsp

PARENTS = {1:1,10239:1,10:10239,11:10,12:10,20:10239,81077:1,99:81077}

def hit(taxon=11, subject='g~0~30', query='p/1', frame=1):
    h=Hsp(subject,150,20,21,frame,4 if frame>0 else 18,12 if frame>0 else 10,2,4,'ACD','ACD')
    return Hit(query,taxon,50.,1e-8,3,h)

def row(hs, sid='original'):
    return paired_row(sid,'sample',('A'*21,'C'*21),('I'*21,'?'*21),hs,PARENTS,{81077,99},set())

class AdapterTests(unittest.TestCase):
    def test_original_mates_and_no_synthetic_unaligned_coordinates(self):
        r=row([[],[hit(query='p/2',frame=-1)]])
        self.assertEqual(r['query_seq'],'A'*21)
        self.assertEqual(r['query_qual_rev'],'?'*21)
        self.assertEqual(r['prim_align_ref_start'],'NA')
        self.assertEqual(r['prim_align_query_rc'],'NA')
        self.assertEqual((r['prim_align_ref_start_unclipped_rev'],r['prim_align_ref_end_unclipped_rev']),(30,50))
        self.assertTrue(r['prim_align_query_rc_rev'])
        self.assertEqual(r['prim_align_best_alignment_score_rev'],'NA')
        self.assertEqual(r['prim_align_pair_status'],'NA')
        self.assertEqual(r['protein_bitscore_rev'],50.)

    def test_split_genomes_keep_mate_order_for_production_to_canonicalize(self):
        r=row([[hit(subject='z~0~30')],[hit(taxon=20,subject='a~3~30',query='p/2')]])
        self.assertEqual(r['prim_align_genome_id_all'],'z/a')
        self.assertEqual(r['aligner_taxid_lca'],'10239')
        self.assertFalse(r['prim_align_query_rc'])
        self.assertTrue(r['prim_align_query_rc_rev'])
        self.assertEqual(r['query_seq_rev'],'C'*21)

    def test_all_near_best_taxa_contribute_and_primary_is_order_independent(self):
        a=hit(taxon=11,subject='a~0~30');b=hit(taxon=12,subject='b~0~30')
        left=row([[a,b],[]]);right=row([[b,a],[]])
        self.assertEqual(left,right)
        self.assertEqual(left['aligner_taxid_lca'],'10')
        self.assertEqual(left['protein_subject'],'a~0~30')
        self.assertEqual(left['aligner_n_assignments'],'2')

    def test_batched_group_matches_production_updates(self):
        retained=[dataclasses.replace(hit(taxon=t),score=score) for t,score in
                  [(99,60.),(12,60.),(11,60.),(11,54.),(12,54.)]*17]
        for artificial,unclassified in [({99},{12}),({99},set()),(set(),{11,12})]:
            expected=lca.new_group('original',retained[0].taxid,retained[0].score,artificial,unclassified)
            for h in retained[1:]:
                expected=lca.update_group(expected,h.taxid,h.score,artificial,unclassified)
            actual=group_from_hits('original',retained,artificial,unclassified)
            self.assertEqual(actual,expected)
            self.assertEqual(lca.output_group(lca.summarize_group(actual,PARENTS,{})[0]),
                             lca.output_group(lca.summarize_group(expected,PARENTS,{})[0]))

    def test_streaming_keeps_boundary_hits_and_rejects_unsorted_pairs(self):
        def record(query,score):
            return '\t'.join(map(str,[query,'g~0~0',score,1e-8,25,100,75,1,75,1,25,1,30,'A'*25,'A'*25]))+'\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits.tsv'
            p.write_text(record('p/1',89)+record('p/1',90)+record('p/1',100)+record('p/2',50)+record('q/1',100))
            observed=list(grouped_hits(p,{'g':(11,90)},{'p'}))
            self.assertEqual(len(observed),1)
            self.assertEqual([[h.score for h in mate] for mate in observed[0][1]],[[100,90],[50]])
            p.write_text(record('q/1',100)+record('p/1',100))
            with self.assertRaisesRegex(ValueError,'sorted'):
                list(grouped_hits(p,{'g':(11,90)},{'p','q'}))

    def test_constructs_follow_production_separate_taxonomy(self):
        r=row([[hit()],[hit(taxon=99,subject='a~0~30',query='p/2')]])
        self.assertEqual(r['aligner_taxid_lca'],'11')
        self.assertEqual(r['aligner_taxid_lca_combined'],'1')
        self.assertEqual(r['aligner_taxid_lca_artificial'],'99')
        self.assertEqual(row([[hit()],[]],sid='another')['seq_id'],'another')
        with self.assertRaises(ValueError):row([[],[]])
        wrong=dataclasses.replace(hit(),hsp=dataclasses.replace(hit().hsp,query_length=22))
        with self.assertRaises(ValueError):row([[wrong],[]])

if __name__=='__main__':unittest.main()
