#!/usr/bin/env python3
"""Write small synthetic inputs for checking the benchmark's plumbing locally.

Not realistic data: just enough schema for each step's script to run. Writes
the inputs, a toy ref_dir and a params file (toy-params.json) into OUTDIR.
"""

import argparse
import gzip
import json
import os
import random

BASES = "ACGT"


def write_gz(path: str, lines: list[str]) -> None:
    """Write lines to a gzipped text file."""
    with gzip.open(path, "wt") as f:
        f.writelines(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir")
    parser.add_argument("--pairs", type=int, default=20000)
    args = parser.parse_args()
    out = os.path.abspath(args.outdir)
    os.makedirs(f"{out}/ref/results", exist_ok=True)
    rng = random.Random(1)

    genomes = [
        (f"NC_{i:06d}.1", str(10000 + i), str(20000 + i // 3)) for i in range(50)
    ]
    write_gz(
        f"{out}/ref/results/virus-genome-metadata-gid.tsv.gz",
        ["genome_id\ttaxid\tspecies_taxid\n"]
        + [f"{g}\t{t}\t{s}\n" for g, t, s in genomes],
    )
    write_gz(
        f"{out}/ref/results/total-virus-db-annotated.tsv.gz",
        ["taxid\ttaxid_species\n"] + [f"{t}\t{s}\n" for _, t, s in genomes],
    )

    # Paired Bowtie2-style SAM sorted by read ID, plus matching interleaved FASTQ
    sam, fastq = [], []
    for i in range(args.pairs):
        qname = f"READ{i:09d}"
        genome = rng.choice(genomes)[0]
        pos = rng.randint(1, 5000)
        seqs = ["".join(rng.choices(BASES, k=150)) for _ in range(2)]
        score = rng.randint(150, 300)
        for mate, flag, p, pn, tlen in (
            (0, 99, pos, pos + 200, 350),
            (1, 147, pos + 200, pos, -350),
        ):
            tags = f"AS:i:{score}\tXS:i:{score - 10}\tNM:i:{rng.randint(0, 5)}\tYT:Z:CP\tYS:i:{score}"
            sam.append(
                f"{qname}\t{flag}\t{genome}\t{p}\t42\t150M\t=\t{pn}\t{tlen}\t{seqs[mate]}\t{'F' * 150}\t{tags}\n"
            )
            fastq.append(f"@{qname}/{mate + 1}\n{seqs[mate]}\n+\n{'F' * 150}\n")
    write_gz(f"{out}/toy.sam.gz", ["@HD\tVN:1.0\tSO:queryname\n"] + sam)
    write_gz(f"{out}/toy_reads.fastq.gz", fastq)

    # Hit tables for the downstream steps
    hits = []
    for i in range(args.pairs):
        taxid = rng.choice(genomes)[1]
        hits.append((f"READ{i:09d}", taxid, "".join(rng.choices(BASES, k=150))))
    header = "seq_id\taligner_taxid_lca\tquery_seq\n"
    write_gz(
        f"{out}/hits_by_seq_id.tsv.gz",
        [header] + [f"{s}\t{t}\t{q}\n" for s, t, q in hits],
    )
    by_taxid = sorted(hits, key=lambda h: h[1])
    write_gz(
        f"{out}/hits_by_taxid.tsv.gz",
        [header] + [f"{s}\t{t}\t{q}\n" for s, t, q in by_taxid],
    )
    db = sorted({(t, s if int(t) % 5 else "NA") for _, t, s in genomes})
    write_gz(
        f"{out}/db_by_taxid.tsv.gz",
        ["aligner_taxid_lca\ttaxid_species\n"] + [f"{t}\t{s}\n" for t, s in db],
    )
    species = dict(db)
    write_gz(
        f"{out}/joined.tsv.gz",
        ["seq_id\taligner_taxid_lca\tquery_seq\ttaxid_species\n"]
        + [f"{s}\t{t}\t{q}\t{species[t]}\n" for s, t, q in by_taxid],
    )

    params = {
        "ref_dir": f"{out}/ref",
        "fvs_sam": f"{out}/toy.sam.gz",
        "fvs_reads": f"{out}/toy_reads.fastq.gz",
        "pvbs_sam": f"{out}/toy.sam.gz",
        "ctd_tsv": f"{out}/hits_by_seq_id.tsv.gz",
        "jt_tsv1": f"{out}/hits_by_taxid.tsv.gz",
        "jt_tsv2": f"{out}/db_by_taxid.tsv.gz",
        "act_tsv": f"{out}/joined.tsv.gz",
    }
    with open(f"{out}/toy-params.json", "w") as f:
        json.dump(params, f, indent=2)


if __name__ == "__main__":
    main()
