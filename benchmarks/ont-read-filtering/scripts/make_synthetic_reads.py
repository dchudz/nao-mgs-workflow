#!/usr/bin/env python3
"""Write random ONT-like reads for testing the benchmark plumbing.

The reads exercise each filter: lengths fall on both sides of the 50 bp and
15 kb limits, mean error rates straddle 10%, and some reads carry
low-complexity stretches for masking. Usage:
  make_synthetic_reads.py N_READS SEED > reads.fastq
"""

import math
import random
import sys


def phred(error: float) -> int:
    """Convert an error probability to a Phred score capped at 50."""
    return min(50, max(1, round(-10 * math.log10(error))))


def make_read(rng: random.Random) -> tuple[str, str]:
    """Return one random read sequence and its quality string."""
    length = int(min(40000, max(20, rng.lognormvariate(7.2, 1.0))))
    seq = [rng.choice("ACGT") for _ in range(length)]
    if rng.random() < 0.3 and length > 300:
        start = rng.randrange(0, length - 250)
        unit = rng.choice(["A", "T", "AT", "CA", "GGC"])
        run = (unit * 200)[: rng.randint(40, 240)]
        seq[start : start + len(run)] = run
    mean_error = rng.uniform(0.02, 0.16)
    quals = [
        chr(33 + phred(min(0.75, rng.expovariate(1 / mean_error))))
        for _ in range(length)
    ]
    return "".join(seq), "".join(quals)


def main() -> None:
    """Print N_READS synthetic FASTQ records to stdout."""
    n_reads, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    out = sys.stdout
    for i in range(n_reads):
        seq, qual = make_read(rng)
        out.write(f"@synthetic-{seed}-{i:07d} runid=test\n{seq}\n+\n{qual}\n")


if __name__ == "__main__":
    main()
