#!/usr/bin/env python3
"""Summarize an uncompressed FASTQ stream from stdin without keeping any reads.

Prints a header and one tab-separated row:
  variant, label, reads, bases, masked_bases (N), md5_ordered,
  hash_unordered, hash_ids_unordered

md5_ordered covers the exact byte stream. hash_unordered is the sum, modulo
2^64, of a per-record hash, so it matches whenever two outputs hold the same
records in any order. hash_ids_unordered does the same for read IDs only,
which shows whether two filters kept the same reads.
"""

import hashlib
import sys

MASK64 = (1 << 64) - 1


def record_hash(data: bytes) -> int:
    """Return the first 8 bytes of the record's MD5 as an integer."""
    return int.from_bytes(hashlib.md5(data).digest()[:8], "big")


def main() -> None:
    """Read FASTQ records from stdin and print the summary row."""
    variant, label = sys.argv[1], sys.argv[2]
    stream_md5 = hashlib.md5()
    reads = bases = masked = 0
    h_records = h_ids = 0
    lines = iter(sys.stdin.buffer)
    for header in lines:
        seq, plus, qual = next(lines), next(lines), next(lines)
        record = header + seq + plus + qual
        stream_md5.update(record)
        reads += 1
        seq = seq.rstrip(b"\n")
        bases += len(seq)
        masked += seq.count(b"N")
        h_records = (h_records + record_hash(record)) & MASK64
        h_ids = (h_ids + record_hash(header.split()[0])) & MASK64
    print(
        "variant\tlabel\treads\tbases\tmasked_bases\tmd5_ordered\thash_unordered\thash_ids_unordered"
    )
    print(
        f"{variant}\t{label}\t{reads}\t{bases}\t{masked}\t{stream_md5.hexdigest()}\t{h_records:016x}\t{h_ids:016x}"
    )


if __name__ == "__main__":
    main()
