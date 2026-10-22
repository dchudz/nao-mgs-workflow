# Amino-acid screening

`screen_amino_acids` is an experimental companion to the nucleotide Nucleaze
screen. It builds an exact peptide k-mer index from genome FASTA and partitions
interleaved FASTQ pairs using that index. The intended first evaluation is on
Nucleaze-rejected pairs, preserving the existing nucleotide positives. INDEX,
RUN and DOWNSTREAM do not invoke this tool yet.

The commands use the standard genetic code (NCBI code 1) and all three reading
frames on both DNA strands, for both reference and query. They do not infer ORFs
or require protein annotations. A stop or any codon containing an ambiguous
base breaks a peptide run. No window crosses a FASTA record, a read boundary,
or the boundary between mates. Lowercase input is accepted; U is treated as T.
Reference masking with N therefore remains effective. Alternative genetic
codes and annotated translation exceptions are not supported.

This is exact amino-acid matching, not protein alignment or taxonomic
classification. Translating all frames also indexes noncoding translations.
Candidate volume, sensitivity, index memory and screening throughput need
measurement before workflow integration. No amino-acid k has been selected or
validated: `--k` is required and accepts 1–12 amino acids (five bits per residue).
The low end exists for small correctness fixtures, not as a recommended setting.

## Local use

Build just this executable:

```bash
cargo build --release --manifest-path rust-tools/Cargo.toml -p screen_amino_acids
```

Choose `K` explicitly for the experiment, then build and use its index:

```bash
rust-tools/target/release/screen_amino_acids index \
    --reference virus-genomes-human-masked.fasta.gz \
    --output virus-genomes.aa.bin --k "$K"

rust-tools/target/release/screen_amino_acids screen \
    --index virus-genomes.aa.bin \
    --reads nucleotide-rejected.interleaved.fastq.gz \
    --matched protein-rescued.fastq \
    --unmatched still-rejected.fastq
```

Use the actual masked genome source that produced the nucleotide index when
comparing screens. A `.nucleaze.bin` file cannot be used as an amino-acid index.
Record the reference and amino-acid index SHA-256 digests, pipeline revision,
tool version (`--version`), command arguments, input digests, and Nucleaze
version/settings with each experiment. Keep masking and taxonomic scope fixed;
do not substitute an unrelated public protein database in the comparison.

Inputs can be plain or gzip, detected by content; concatenated gzip members are
supported. `--reads2` accepts a separate mate file; otherwise input is interleaved.
`--threads` controls parallel matching while preserving pair order.
`--unmatched` is optional to avoid writing the rejected fraction.
Output FASTQs are uncompressed and preserve every original header,
sequence, separator and quality line. Existing outputs are never overwritten.
Normal failures remove incomplete outputs. Files are visible during execution:
consume them only after the command exits successfully. Forced termination may
leave incomplete outputs to remove before retrying. This requires ordinary
file creation and writing, without hard-link or atomic-rename requirements.

Pairs must be adjacent in R1/R2 order, with matching first header tokens after
removing optional `/1` and `/2` suffixes. When suffixes or recognized mate
comments (`1`, `2`, `1:...`, `2:...`) are present, they must agree with position.
Missing mate annotations are allowed. Four-line FASTQ is required. Truncation,
an odd record count, mismatched pair IDs, invalid sequence characters, and
invalid quality lengths fail the command. Duplicate pair names are retained
by occurrence, not collapsed. Any matching k-mer in either mate sends both
records to `--matched` exactly once; the other pairs go to `--unmatched`.

The index has a versioned header recording amino-acid k, genetic code and frame
policy, followed by sorted unique encoded peptide words. It is deterministic
for the same peptide set regardless of reference ordering. Incompatible headers,
invalid or duplicate codes, truncation, trailing bytes and empty indexes fail.
The format is separate from Nucleaze's index. It is not a checksum or a record
of the source genomes; retain the provenance described above.

The tool reports reference-record and unique-k-mer counts when indexing, and
matched/total **pair** counts when screening. Screening an already rejected
fraction gives a disjoint rescue set by construction. It does not determine
which rescues survive the existing cleaning, alignment and contamination gates,
or establish improved final detection. Those remain separate measurements.
Do not combine rescue candidates with the original FASTP input without testing
effects on its automatic adapter detection.

## Correctness tests

```bash
cargo test --manifest-path rust-tools/Cargo.toml -p screen_amino_acids
```

The tests use short artificial sequences that are not organism-derived.
An independent explicit codon-table oracle checks all 4,096 six-base strings,
plus stops, ambiguities, lower case, U, short inputs and the maximum k.
Executable tests cover synonymous peptide matches, either-mate selection,
all six frames, exact FASTQ partitioning, malformed input, gzip integrity,
deterministic index serialization and output protection.
