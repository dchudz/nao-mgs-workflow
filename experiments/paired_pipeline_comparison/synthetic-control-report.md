# Synthetic known-origin adapter control

The exact-production-container control passed: 24 positive pairs retained their
prespecified genomic coordinates and taxonomic labels, and reduced to 12 groups
after both production duplicate passes. One negative pair produced no candidate.
This is an implementation check, not evidence of environmental sensitivity.

## Prespecified fixture

Expectations and artificial inputs were committed in `b0d691ba` before execution.
All sequences were generated from random coding strings; none came from an
organism. The six 900-base references have fixture-only species labels
900000000 through 900000005, named Synthetic control 0 through 5. The database
contains one controlled peptide per reference, covering six reference frames.
It is deliberately not the full production reference.

Each reference contributes two distinct paired-read origins. Origin A appears
three times: a high-quality exact observation, a lower-quality exact repeat,
and a synonymous variant. Origin B appears once. Both mates and original quality
strings must survive reconstruction. Terminal ambiguous bases exercise clipping;
query orientation and phase cover all six translated query frames.

The complete independent answers, including every mate's DNA coordinates,
orientation, label, sequence, quality and expected exemplar, are in
[`synthetic-control-fixture/expected.json`](synthetic-control-fixture/expected.json).
Intervals were checked against literal reference slices before importing any
adapter code. Expected results are not calculated by the coordinate converter.

| Check | Expected | Observed |
|---|---:|---:|
| Input pairs | 25 | 25 |
| Exact unique query pairs, including negative | 19 | 19 |
| Positive pairs after adapter reconstruction | 24 | 24 |
| Negative pairs excluded | 1 | 1 |
| Coordinate duplicate groups | 12 | 12 |
| Groups after similarity duplicate handling | 12 | 12 |
| Independent geometry cases | 36 | 36 passed |
| Deliberate output corruptions rejected | 3 | 3 |

## Execution and provenance

The test ran the existing query/membership builder, real DIAMOND alignment,
serial and parallel adapter CLIs, parallel candidate extraction, coordinate
duplicate marking, sorting by read ID, and similarity duplicate marking.
Taxonomy was assigned by the existing production LCA implementation. Search
settings were those in the frozen protein alignment process; no parameter sweep
or tuning was performed. Adapter, search and duplicate implementation code was
unchanged. Only the test harness and future execution guard were added.

The observed rows were checked against the independent answers, not just against
one another. The checker also rejected an intentionally shifted coordinate,
wrong taxonomic label and wrong duplicate group. The geometry check directly
exercised 36 combinations of reference frame and query frame with known bounds.

Attempts remain saved in `synthetic-control-results/`:

1. The initial fixture omitted a child required under the artificial-sequence
   taxonomy category. The parser rejected it before adapter output. Commit
   `7a110837` corrected the fixture structure; the expected answers and read
   sequences did not change.
2. The host-executable run passed. Subsequent verification showed that those
   duplicate binaries were not byte-identical to the production container, so
   this attempt was insufficient for acceptance.
3. Commit `8f1a7376` made container selection explicit. The same fixture passed
   using immutable image
   `sha256:d1b77ec75262ad7dc27dce983b4294782c84d986529aabe335c94d8f9aab7b12`.
   The harness reported 1.57 seconds, including its provenance checks. This is
   not an estimate for the real-data comparison.

[`synthetic-control-results/attempt-3/result.json`](synthetic-control-results/attempt-3/result.json)
records source hashes, executable hashes, expected-answer hashes and the tested
source commit (`629db5c2`). Raw artificial alignments, reconstructed rows, duplicate
assignments, exact commands and logs are retained alongside it.

## Acceptance and limits

Future comparison starts check this result before launching stages, and reject
missing evidence or a mismatch in the recorded fixture, checker, tested source,
DIAMOND binary or duplicate-container image. The already-started first comparison
requires retrospective checking of the same identities; this control was not
run before that job launched. A pass satisfies only this control, not all other
comparison acceptance checks.

This proves the declared cases work. It does not prove every possible alignment
is handled correctly. It does not test full six-frame reference construction,
competing-reference ambiguity, indels in the real-search fixture, upstream
screening/filtering, validation searches, or biological truth. Existing geometry
unit tests separately exercise gaps and invalid coordinates. Artificial taxonomy
labels test record handling without implying a biological classification.

The earlier description of coverage as only serial/parallel byte parity was
incomplete: coordinate unit tests and a small synthetic duplicate-marker test
already existed. The new evidence adds independently frozen origins carried
through real alignments and the exact production duplicate container.

Separate reference-construction evidence also exists: the archived [builder](../amino_acids/history/source/translate-alignment-reference.original.py) checked its translation and offsets against a separate implementation on 503 artificial inputs. A [prior reproduction](../amino_acids/end-to-end-full-depth/reference-reproduction-verification.json) regenerated the complete protein FASTA from the frozen native alignment genomes and matched the archived bytes. These checks complement this fixture; neither establishes biological accuracy or exercises every builder/adapter combination. The saved builder hash was rechecked against the retained worker during this run. No new reference was generated for this audit.
