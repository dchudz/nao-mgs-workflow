"""Map one real translated HSP into production duplicate-coordinate fields.

Aligned strings must come from DIAMOND qseq_gapped and sseq_gapped.
The ordinary qseq field is nucleotide sequence in blastx and is not suitable.

This only supplies geometry. It does not assign taxonomy, choose primary hits,
merge mates, invent nucleotide alignment scores, or manufacture SAM records.
"""
from dataclasses import dataclass
from translated_coordinates import coordinates

@dataclass(frozen=True)
class Hsp:
    subject: str
    genome_length: int
    peptide_length: int
    query_length: int
    query_frame: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_aligned: str
    subject_aligned: str

@dataclass(frozen=True)
class Geometry:
    genome: str
    ref_start: int
    ref_end_exclusive: int
    ref_start_unclipped: int
    ref_end_unclipped: int  # production convention: zero-based inclusive
    query_rc: bool

def geometry(h: Hsp) -> Geometry:
    """Require translated residue counts and frame phase to match reported bounds."""
    if h.query_frame not in (-3,-2,-1,1,2,3):
        raise ValueError('A six-frame query alignment requires frame ±1, ±2, or ±3')
    try:
        _, frame, offset = h.subject.rsplit('~',2)
        frame, offset = int(frame), int(offset)
    except (ValueError,TypeError) as e:
        raise ValueError('Expected the preserved genome~frame~offset peptide ID') from e
    if not (0<=frame<6 and offset>=0 and offset%3==frame%3):
        raise ValueError('Reference peptide frame and nucleotide offset disagree')
    if not (1<=h.subject_start<=h.subject_end<=h.peptide_length):
        raise ValueError('Subject coordinates are outside the preserved peptide')
    qlo,qhi=min(h.query_start,h.query_end)-1,max(h.query_start,h.query_end)
    if not (0<=qlo<qhi<=h.query_length):
        raise ValueError('Query alignment bounds exceed the original read')
    if (h.query_end>h.query_start)!=(h.query_frame>0):
        raise ValueError('Query coordinate direction disagrees with its frame')
    phase=h.query_start-1 if h.query_frame>0 else h.query_length-h.query_start
    if phase%3!=abs(h.query_frame)-1:
        raise ValueError('Query start does not match the declared codon phase')
    q,s=h.query_aligned,h.subject_aligned
    if not q or len(q)!=len(s) or any(a==b=='-' for a,b in zip(q,s)):
        raise ValueError('Aligned strings must have equal length and no double-gap columns')
    if set(q+s)-set('ABCDEFGHIJKLMNOPQRSTUVWXYZ*-'):
        raise ValueError('Frameshift or malformed alignment strings are unsupported')
    if 3*len(q.replace('-',''))!=qhi-qlo or len(s.replace('-',''))!=h.subject_end-h.subject_start+1:
        raise ValueError('Aligned residue counts do not reproduce the reported bounds')
    if offset+3*h.peptide_length>h.genome_length:
        raise ValueError('Preserved peptide extends beyond its original genome')
    genome,strand,lo,hi=coordinates(h.subject,h.subject_start,h.subject_end,h.genome_length)
    reverse=(h.query_frame<0)!=(strand=='-')
    left_clip,right_clip=(h.query_length-qhi,qlo) if reverse else (qlo,h.query_length-qhi)
    return Geometry(genome,lo,hi,lo-left_clip,hi-1+right_clip,reverse)
