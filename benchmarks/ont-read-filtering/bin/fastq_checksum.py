#!/usr/bin/env python3
"""Summarize an uncompressed FASTQ stream from stdin without keeping any reads.

Prints a header and one tab-separated row:
  variant, label, reads, bases, masked_bases (N), md5_ordered,
  hash_unordered, hash_ids_unordered, hash_seqs_unordered,
  hash_quals_outside_n_unordered

md5_ordered covers the exact byte stream. hash_unordered is the sum, modulo
2^64, of a per-record hash, so it matches whenever two outputs hold the same
records in any order. The other hashes do the same for parts of each record:
  - hash_ids_unordered: read IDs only, which shows whether two filters kept
    the same reads.
  - hash_seqs_unordered: read ID and sequence, ignoring quality, which shows
    whether two maskers masked the same bases.
  - hash_quals_outside_n_unordered: read ID and quality string, ignoring the
    quality of N bases. If this and hash_seqs_unordered both match, two outputs
    differ at most in the quality scores of masked bases.
"""

import hashlib
import re
import sys

MASK64 = (1 << 64) - 1
N_RUN = re.compile(rb"N+")


def record_hash(data: bytes) -> int:
    """Return the first 8 bytes of the record's MD5 as an integer."""
    return int.from_bytes(hashlib.md5(data).digest()[:8], "big")


def blank_n_quals(seq: bytes, qual: bytes) -> bytes:
    """Return qual with the score at every N base of seq replaced by '!'."""
    if b"N" not in seq:
        return qual
    blanked = bytearray(qual)
    for run in N_RUN.finditer(seq):
        blanked[run.start() : run.end()] = b"!" * (run.end() - run.start())
    return bytes(blanked)


def main() -> None:
    """Read FASTQ records from stdin and print the summary row."""
    variant, label = sys.argv[1], sys.argv[2]
    stream_md5 = hashlib.md5()
    reads = bases = masked = 0
    h_records = h_ids = h_seqs = h_quals = 0
    lines = iter(sys.stdin.buffer)
    for header in lines:
        seq, plus, qual = next(lines), next(lines), next(lines)
        record = header + seq + plus + qual
        stream_md5.update(record)
        reads += 1
        read_id = header.split()[0]
        seq = seq.rstrip(b"\n")
        qual = qual.rstrip(b"\n")
        bases += len(seq)
        masked += seq.count(b"N")
        h_records = (h_records + record_hash(record)) & MASK64
        h_ids = (h_ids + record_hash(read_id)) & MASK64
        h_seqs = (h_seqs + record_hash(read_id + b"\n" + seq)) & MASK64
        h_quals = (
            h_quals + record_hash(read_id + b"\n" + blank_n_quals(seq, qual))
        ) & MASK64
    print(
        "variant\tlabel\treads\tbases\tmasked_bases\tmd5_ordered\thash_unordered\thash_ids_unordered"
        "\thash_seqs_unordered\thash_quals_outside_n_unordered"
    )
    print(
        f"{variant}\t{label}\t{reads}\t{bases}\t{masked}\t{stream_md5.hexdigest()}\t{h_records:016x}\t{h_ids:016x}"
        f"\t{h_seqs:016x}\t{h_quals:016x}"
    )


if __name__ == "__main__":
    main()
