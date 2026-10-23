"""Protein HSP geometry obeys the pipeline's genomic endpoint convention."""
import dataclasses,unittest
from diamond_geometry import Hsp,geometry

class GeometryTests(unittest.TestCase):
    def test_six_reference_frames_and_both_query_strands(self):
        # Original read has 3 left and 9 right clipped nt in its translated orientation.
        for frame in range(6):
            offset=30+frame%3
            for qframe,start,end in [(1,4,12),(-1,18,10)]:
                h=Hsp(f'g~{frame}~{offset}',150,20,21,qframe,start,end,2,4,'ACD','ACD')
                x=geometry(h)
                left,right=(offset+3,offset+12) if frame<3 else (150-offset-12,150-offset-3)
                # Strand-relative clips reverse with the subject's orientation.
                clips=(3,9) if frame<3 else (9,3)
                self.assertEqual((x.ref_start,x.ref_end_exclusive),(left,right))
                self.assertEqual((x.ref_start_unclipped,x.ref_end_unclipped),(left-clips[0],right-1+clips[1]))
                self.assertEqual(x.query_rc,(frame>=3)!=(qframe<0))

    def test_gaps_and_nonzero_query_phase_preserve_reference_span(self):
        for qframe,start,end in [(2,5,13),(3,6,14),(-2,19,11),(-3,18,10)]:
            h=Hsp('g~0~30',150,20,23,qframe,start,end,2,4,'A-CD','AB-D')
            x=geometry(h)
            self.assertEqual((x.ref_start,x.ref_end_exclusive),(33,42))
            self.assertEqual(x.ref_end_unclipped-x.ref_start_unclipped+1,23)
        # A subject deletion consumes more reference than query, so the unclipped span grows.
        h=Hsp('g~0~0',60,20,12,1,1,9,1,4,'A-CD','ABCD')
        x=geometry(h)
        self.assertEqual((x.ref_start_unclipped,x.ref_end_unclipped),(0,14))

    def test_invalid_coordinates_are_not_converted_to_duplicate_keys(self):
        h=Hsp('g~0~30',150,20,21,1,4,12,2,4,'ACD','ACD')
        for changes in [dict(query_frame=0),dict(query_start=5),dict(query_end=11),dict(subject='g~1~30'),dict(subject_end=21),dict(query_aligned='A/CD'),dict(subject_aligned='AC'),dict(genome_length=35)]:
            with self.subTest(changes=changes),self.assertRaises(ValueError):geometry(dataclasses.replace(h,**changes))

if __name__=='__main__':unittest.main()
