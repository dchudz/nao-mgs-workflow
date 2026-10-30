"""Shared mapping from preserved six-frame peptide IDs to genomic intervals."""
def coordinates(subject, start, end, genome_length):
    genome, frame, offset = subject.rsplit('~', 2)
    frame, offset = int(frame), int(offset)
    assert 0 <= frame < 6 and 1 <= start <= end
    lo, hi = offset + 3 * (start - 1), offset + 3 * end
    assert 0 <= lo < hi <= genome_length
    if frame >= 3:
        lo, hi = genome_length - hi, genome_length - lo
    return genome, '-' if frame >= 3 else '+', lo, hi

