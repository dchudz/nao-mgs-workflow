# Searching for viruses in protein space

**A paired comparison of the production wastewater pipeline against a protein-search variant, over 60 sequencing-run comparisons.**

*Written by two AI agents — Codex, which ran the experiment, and Claude, which reviewed and wrote it up. No human has verified the claims in this document. Figures are read from the frozen tables listed in [Provenance](#provenance); the reasoning around them should be checked, not assumed. Experiment run 2026-10-24/25.*

---

## Summary

The production pipeline looks for viral reads by matching DNA against a collection of viral genomes. DNA sequence diverges fast; the protein it encodes diverges slowly. So a virus distant enough from anything in the reference collection should be invisible to a DNA search while still being recognisable as a protein. That is the premise we set out to test.

We built a variant of the pipeline that searches in protein space — a translated word screen unioned with the existing DNA screen, and DIAMOND `blastx` in place of Bowtie2 — and ran both arms over identical reads from 50 public wastewater sequencing runs at 5 million read pairs each, plus ten of those runs again at 25 million.

**The variant reports 55–65% more species calls.** The increase is spread across the cohort: 47 of 50 runs gain, one is unchanged, two lose ground. It holds at all three reporting thresholds we tried, and at both sequencing depths.

**We cannot tell you whether those calls are real.** These are environmental samples with unknown contents, so there are no truth labels anywhere in this data. What we do have is the pipeline's own validation step, which re-searches a sample of output reads against a comprehensive DNA database. Split by category, it says two different things. On the calls both arms make, the variant's reads look much like production's. On the calls the variant adds, they do not: 3 of 1,051 checked rows are compatible with their assignment, at a mean taxonomic distance of 7 steps rather than about 1.

That is consistent with having found genuinely divergent viruses — a divergent virus is exactly what a near-exact DNA search should fail to confirm. It is equally consistent with having matched protein regions that viruses share with other organisms. **The experiment as run cannot separate those two explanations**, and the check that would have — a search against a comprehensive protein database — was specified and costed but not performed. [Section 8](#8-what-we-learned-and-what-we-did-not) says what would settle it.

---

## Contents

1. [Why we tried this](#1-why-we-tried-this)
2. [What we changed](#2-what-we-changed)
3. [How the comparison was run](#3-how-the-comparison-was-run)
4. [What counts as a detection](#4-what-counts-as-a-detection)
5. [What we found: more output](#5-what-we-found-more-output)
6. [What we found: more species calls](#6-what-we-found-more-species-calls)
7. [Whether the new calls look real](#7-whether-the-new-calls-look-real)
8. [What we learned, and what we did not](#8-what-we-learned-and-what-we-did-not)
9. [What it cost](#9-what-it-cost)
10. [Provenance](#provenance)

---

## 1. Why we tried this

### The problem

The pipeline's job is to find viral sequence in wastewater. It does that by comparison against a reference collection of known viral genomes: a read is viral if it looks like something we already have on file.

That design has a known weakness, and it is the weakness that matters most for biosurveillance. A novel or highly divergent virus — the thing the system exists to catch — is by definition not well represented in the reference collection. The question is how much divergence the pipeline tolerates before a read stops looking like anything it knows.

### Why protein space should help

Protein-coding sequence diverges unevenly. The genetic code is redundant: most amino acids have several codons, and the third position of a codon is frequently free to change without changing the amino acid. Selection acts on the protein, so synonymous changes accumulate while the protein sequence stays put. Two homologous viral genes can be 70% identical at the nucleotide level and 95% identical at the amino-acid level.

Both of the pipeline's DNA-matching steps are sensitive to exactly that kind of change:

- The **initial screen** requires an *exact* 24-nucleotide word shared between the read and some reference genome (`nucleaze_k = 24`, `configs/index.config`). A single substitution anywhere in a 24-base window destroys that window. At even moderate divergence, finding 24 consecutive identical bases becomes unlikely.
- The **aligner**, Bowtie2, runs with `--local --very-sensitive-local`, whose preset seeds alignments on 20-nucleotide exact matches with no mismatches permitted inside the seed. An alignment that cannot be seeded is never found, however good the extension would have been.

Searching the same read as protein changes the arithmetic. A 12-residue exact amino-acid match covers 36 nucleotides of the original read and survives every synonymous substitution in that span. And a protein aligner scores with a substitution matrix that knows some amino-acid swaps are chemically conservative and should be cheap, rather than treating every difference as equally disqualifying.

### What we tried first, and why it was not enough

Before this comparison we ran a separate phase that changed **only the screen**, leaving Bowtie2 in place, across 1,000 samples and 5 billion input read pairs. The result is the reason this comparison exists.

| Screen | Final pairs | Added vs. nt24 | Lost | Added / baseline | Adds beyond nt20 | Adds beyond the aa12 union |
|---|---:|---:|---:|---:|---:|---:|
| **nt24** (production) | 1,478,260 | — | — | — | — | — |
| nt22 | 1,493,077 | +14,817 | 0 | +1.00% | 0 | 336 |
| nt20 | 1,498,327 | +20,067 | 0 | +1.36% | 0 | 474 |
| aa12 | 1,495,224 | +19,637 | 2,673 | +1.33% gross, +1.15% net | 44 | 0 |
| aa12 ∪ nt24 | 1,497,897 | +19,637 | 0 | +1.33% | 44 | 0 |

*Source: `experiments/amino_acids/history/results/first-step-1000/heldout/summary.json`, 1,000 completed samples; last changed in commit `864cd3af`.*

Three things follow.

**The amino-acid screen's gain was matched by simply shortening the DNA word.** It added 1.33% more final pairs gross, 1.15% net of its own 2,673 losses. Dropping the exact DNA word from 24 letters to 20 added 1.36%, with no losses at all.

**Almost none of its gain was unique.** Measured against the shorter DNA screen, the amino-acid screen contributed 44 additional pairs out of 1.48 million. The shorter DNA screen contributed 474 beyond the amino-acid union — an order of magnitude more.

**The binding constraint was downstream.** If a shorter nucleotide screen recovers nearly everything the protein screen finds, then the screen is not what is limiting the pipeline. The reads were getting through the screen and then failing to align. That points at Bowtie2, and it is consistent with Bowtie2's 20-nucleotide exact seed requirement, though it does not prove the seed is the specific cause.

So the present comparison changes the screen **and** the aligner together, deliberately. That buys a measurable effect at the cost of not being able to attribute it to either part — a trade we return to in [section 8](#8-what-we-learned-and-what-we-did-not).

### What we expected

The expectation going in was a low single-digit percentage increase in output, consistent with the earlier phase. A result of 55–65% is far outside that, which is itself a reason to look hard at where the extra output comes from rather than to treat it as a win.

---

## 2. What we changed

### The pipeline, stage by stage

The production path for viral reads runs roughly like this. Raw read pairs arrive. A fast k-mer screen throws away everything that shares no exact 24-base word with any reference viral genome; this is a cheap filter that reduces hundreds of millions of pairs to tens of thousands. FASTP then trims adapters and low-quality bases. Bowtie2 aligns what survives against the viral reference and assigns a taxon. The read is then checked against human and other contaminant references and discarded if it matches one. Duplicates are collapsed in two passes — first by alignment coordinates, then by sequence similarity. What remains is counted per taxon and annotated. Finally a validation step takes a sample of output reads and re-searches them against core_nt, a comprehensive DNA database, to see whether the assignment holds up.

The variant replaces two of those stages and leaves the rest alone.

### Every stage, and whether it changed

| # | Stage | production | protein variant | Changed? |
|---|---|---|---|---|
| 0 | Source viral reference | Viral nucleotide genome collection | The same source genomes, translated in six reading frames | **Derived representation** |
| 1 | Initial screen | `nucleaze` exact 24-nt word match against viral genomes | **Union** of that screen and an exact 12-residue amino-acid word match | **Yes** |
| 2 | Adapter and quality trimming | FASTP on screened pairs | *same code and settings* | No |
| 3 | Viral alignment | Bowtie2 `--local --very-sensitive-local --score-min G,0.1,19 -k 10 -X 850` | DIAMOND `blastx --very-sensitive --algo 0 --min-orf 1 --query-gencode 1 --evalue 1e-5 --query-cover 30 --top 10` | **Yes** |
| 4 | Alignment record conversion | Production alignment-to-hit-table path | **New adapter** supplying reference coordinates, paired-read fields and taxonomic assignments in the expected hit-table schema | **Yes** |
| 5 | Human and contaminant removal | Bowtie2 against human and other-contaminant references | *same code and settings* | No |
| 6 | Duplicate handling | Alignment-coordinate marking, then similarity marking | *same code and settings* | No |
| 7 | Downstream counting and annotation | Clade counts and annotated hit records | *same code and settings* | No |
| 8 | Validation | megablast against core_nt, capped sample | *same code and settings* | No |

*Sources: `subworkflows/local/extractViralReadsShort/main.nf`, `configs/index.config`, `experiments/amino_acids/protein-production.nf`, `experiments/amino_acids/protein-parallel.nf`, `experiments/paired_pipeline_comparison/run_pair.py`, `workflows/downstream.nf`, at frozen source commit `c77cdde0`.*

**"Same code and settings" does not mean the same output.** Every unchanged stage receives different candidate reads in the two arms, so its results differ. Stage 5 is worth being explicit about: the human and contaminant filters *do* run on protein-arm reads, identically to production. They are not bypassed.

### Stage 1: the new screen

The amino-acid screen translates each read in all six reading frames and looks for an exact 12-residue word shared with the translated reference.

Twelve is not a tuned parameter. It is the implementation's ceiling. The screen packs each residue into 5 bits and stores a window as a single `u64`:

```rust
// rust-tools/screen_amino_acids/src/lib.rs
let mask = (1u64 << (5 * k)) - 1;
```

Five bits times twelve residues is 60 of the available 64. Thirteen would need 65. The code enforces it directly:

```rust
if !(1..=12).contains(&k) {
    return Err(invalid("amino-acid k must be between 1 and 12"));
}
```

So k=12 is the largest word the current implementation can represent, not a value chosen because it performed best. Nothing here establishes that 12 is the right length — only that it is the length available.

**The variant takes the union of both screens rather than substituting.** This matters, and the earlier phase is why: used alone, the amino-acid screen *loses* 2,673 pairs that the production screen keeps. The two screens overlap only partly, and neither contains the other. Substituting one for the other would drop reads the pipeline currently finds.

### Stage 3: the new aligner

Bowtie2 is a short-read DNA aligner built for mapping reads to a reference they closely resemble. It finds a seed — an exact 20-base match under the configured preset — and extends it, scoring matches and mismatches at the nucleotide level, requiring a final score above a floor that scales with read length (`--score-min G,0.1,19`).

DIAMOND `blastx` translates each read in six frames and aligns the resulting peptides against a protein reference. It is also seed-based — it uses spaced seeds over a reduced amino-acid alphabet, so it is not the case that it does away with seeding — but its scoring uses a protein substitution matrix, so a chemically conservative amino-acid change costs little.

**The two aligners accept alignments on criteria that are not comparable.** Bowtie2 accepts on a nucleotide alignment score. DIAMOND accepts on an E-value of 1e-5 with the alignment covering at least 30% of the query, reporting all hits within 10% of the best bit score (`--top 10`). These are different kinds of threshold and they were not calibrated against each other. A share of the output difference is therefore attributable to the acceptance criteria rather than to protein space as such, and this experiment cannot say how large that share is.

One further setting is worth naming: `--min-orf 1` tells DIAMOND not to require a long open reading frame before searching a frame. Short reads frequently contain stop codons in some frames; requiring a long ORF would discard them. Setting it to 1 means essentially every frame of every read is searched.

### Stage 4: the adapter

DIAMOND's output is not in the format the production downstream workflow consumes, so the variant adds a conversion step. It reconstructs reference coordinates, paired-read fields and taxonomic assignments into the expected hit-table schema. This is new code on the critical path, which is why it was checked against prespecified answers before the comparison was allowed to proceed.

<details>
<summary><strong>Known-answer control for the adapter</strong></summary>

The adapter was checked against answers specified before execution, using artificial coding sequences rather than an organism's genome. The check exercised the real adapter and the production duplicate-marking container.

| Known-answer control | Expected | Observed |
|---|---:|---:|
| Positive read pairs reconstructed | 24 | 24 |
| Negative pair excluded | 1 | 1 |
| Distinct groups after both duplicate passes | 12 | 12 |
| Coordinate / orientation cases | 36 | 36 passed |
| Deliberately corrupted outputs rejected | 3 | 3 |

*Source: `experiments/paired_pipeline_comparison/synthetic-control-report.md`, prespecified fixture commit `b0d691ba`, container-selection correction `8f1a7376`. The first comparison had started before this check was requested and was checked retrospectively; all subsequent execution was gated on the verified control.*

This tests adapter and integration correctness. It does not test sensitivity, and it is not a biological control.

</details>

<details>
<summary><strong>Frozen settings</strong></summary>

| Frozen setting | Value |
|---|---|
| Nucleotide screen word length | 24 nucleotides, exact |
| Amino-acid screen word length | 12 residues, exact, six-frame |
| Protein search mode | `--very-sensitive --algo 0` |
| Protein search acceptance | `--evalue 1e-5 --query-cover 30 --top 10 --min-orf 1` |
| Nucleotide alignment | `--local --very-sensitive-local --score-min G,0.1,19 -k 10 -X 850` |
| Alignment duplicate coordinate tolerance | 1 nucleotide |
| Validation sample cap | 20 eligible records per taxonomic group within each sample group |
| Added reporting thresholds | 2, 3 and 5 joint exemplars per group/species |

The three reporting thresholds are analysis parameters added for this comparison. They are not pipeline settings — see [section 4](#4-what-counts-as-a-detection).

</details>

---

## 3. How the comparison was run

| Stratum | Sequencing runs | Input per run | Paired comparisons |
|---|---:|---:|---:|
| Main cohort | 50 | 5 million read pairs | 50 |
| Nested depth comparison | 10 of those same 50 | 25 million read pairs | 10 |

All 60 comparisons completed and were accepted.

Runs were drawn at random from a saved inventory of public wastewater sequencing runs from SecureBio's study. Reads are sampled uniformly across each source archive, and **both arms receive exactly the same sampled reads** — the comparison is paired at the read level, not merely at the run level. In the depth comparison the 5M input is contained within the 25M input, so the two strata are nested and must not be pooled as independent observations.

Both arms use the same source viral genome collection and the same downstream settings. Production's viral stage runs first in every pair, so the alignment-stage order is fixed; the downstream validation order is prespecified and counterbalanced across the cohort.

Two caveats on the cohort. These runs are not a random sample of wastewater sequencing in general — they are a random draw from one saved inventory. And whether any two runs derive from the same specimen was not checked against metadata, so their biological independence has not been established either way.

<details>
<summary><strong>Execution environment and acceptance</strong></summary>

Work ran on an AWS `r6i.12xlarge` (48 vCPUs, 384 GiB) in the user's us-east-1 account, with some work on a Spot instance. A comparison counts as accepted only when it has a source identifier, input checksums, execution receipts, the standard downstream reports, and a private archive verified by reading it back. The cohort inventory records completion separately for every run and depth. Worker interruption and resumption records are preserved; resumption excluded already-accepted comparisons and retained the frozen scientific settings.

</details>

---

## 4. What counts as a detection

Before any numbers, a definition problem that changes the headline by orders of magnitude.

**The pipeline does not emit detections.** `countReadsPerClade` emits per-taxon read counts and annotated hit records. There is no minimum-read rule, no binary calling rule, nothing that converts counts into "this species is present" — verified by searching across `modules/`, `subworkflows/`, `workflows/` and `bin/`. Deciding what counts as a detection is a judgement made *after* the pipeline runs.

So this report is careful about three different units:

| Unit | What one of them is | Why it moves |
|---|---|---|
| **Row** | One read pair with one taxonomic assignment in the viral hits table | Inflated by sequencing depth, by duplicate reads, and by abundant species owning thousands of rows each. A read that merely changes its label counts in both arms. |
| **Joint exemplar** | One row surviving both duplicate passes (`seq_id == prim_align_dup_exemplar == sim_dup_exemplar`) | Removes duplicate amplification within a sample. Not monotone: removing records can *raise* exemplar counts, because a different representative survives. |
| **Thresholded call** | One (group, species) pair with at least N joint exemplars | **Not a pipeline output.** N is chosen by the analysis. The earlier 200-sample analysis used N=3, so we report at N=2, 3 and 5 rather than letting one arbitrary N carry the headline. |

Throughout this report, "call" means the third of these, at the stated threshold. It is our reporting rule, not the pipeline's.

---

## 5. What we found: more output

### The funnel

Here is what happens to 250 million input read pairs in each arm, stage by stage. Note that the last three rows change units — read pairs, then duplicate exemplars, then group/species combinations.

**5 million read pairs per input, 50 comparisons**

| Stage / unit | production | protein variant | ratio |
|---|---:|---:|---:|
| Selected input pairs | 250,000,000 | 250,000,000 | — |
| Pairs passing the raw screen | 134,101 | 3,768,248 | **28.1×** |
| Pairs after FASTP, entering alignment | 127,784 | 571,748 | 4.47× |
| Final hit-table rows | 70,084 | 108,862 | 1.55× |
| Joint duplicate exemplars | 67,452 | 99,206 | 1.47× |
| Exemplars assigned at species rank or below | 63,455 | 55,404 | **0.87×** |
| Group/species with at least 3 exemplars | 661 | 1,058 | 1.60× |

**25 million read pairs per input, 10 nested comparisons**

| Stage / unit | production | protein variant | ratio |
|---|---:|---:|---:|
| Selected input pairs | 250,000,000 | 250,000,000 | — |
| Pairs passing the raw screen | 170,557 | 1,824,348 | **10.7×** |
| Pairs after FASTP, entering alignment | 167,211 | 635,155 | 3.80× |
| Final hit-table rows | 91,821 | 129,965 | 1.42× |
| Joint duplicate exemplars | 77,223 | 101,037 | 1.31× |
| Exemplars assigned at species rank or below | 76,203 | 45,760 | **0.60×** |
| Group/species with at least 3 exemplars | 231 | 380 | 1.65× |

*Source: `funnel.tsv`, generated by `writeup_tables.py` (commit `906b4fba`) from accepted result history through `3d847146`. FASTP counts come from the viral-screen task, not the shared community-profiling subset.*

Reading down the table: the protein screen admits one to two orders of magnitude more pairs, quality trimming removes most of them, and each later stage narrows the gap. By the time we reach output rows the difference is 1.5×, not 28×.

**One row runs the other way.** The variant produces more joint exemplars than production but *fewer* exemplars resolved to species rank — 55,404 against 63,455 at 5M, and 45,760 against 76,203 at 25M. It reports more calls while resolving fewer reads to a species.

### The two screens miss different things

| Depth | Nucleotide screen | Amino-acid screen | Caught by both | Nucleotide only | Union |
|---|---:|---:|---:|---:|---:|
| 5M | 134,101 | 3,721,860 (28×) | 87,713 | 46,388 (34.6%) | 3,768,248 |
| 25M | 170,557 | 1,767,401 (10×) | 113,610 | 56,947 (33.4%) | 1,824,348 |

*Source: `screen-funnel.tsv`, all 60 accepted comparisons.*

The amino-acid screen admits an order of magnitude more pairs and still independently misses about a third of what the nucleotide screen admits. This is the empirical justification for taking the union.

It also shows how coarse the amino-acid screen is as a filter. FASTP keeps 95.3% and 98.0% of production's screened pairs, against 15.2% and 34.8% of the variant's. Most of what the protein screen admits is material that fails basic quality filtering — the screen is passing a large amount of low-quality sequence that an exact 24-base DNA match would never have let through.

### Most shared reads are relabelled, not added

This is the result that is easiest to miss and probably the most important for interpreting everything else.

| Measurement | 5M × 50 | 25M × 10 |
|---|---:|---:|
| production rows | 70,084 | 91,821 |
| protein variant rows | 108,862 | 129,965 |
| Pairs present only in production | 412 | 491 |
| Pairs present only in the variant | 39,190 | 38,635 |
| Pairs present in both arms | 69,672 | 91,330 |
| — of those, unchanged assignment | 29,434 | 28,678 |
| — of those, **changed** assignment | 40,238 | 62,652 |
| Net additional rows in the variant | 38,778 | 38,144 |
| **Reassigned share of shared pairs** | **57.8%** | **68.6%** |

*Source: `output-volume.tsv`.*

The variant is not production plus extras. Of the read pairs both arms recover, **most get a different taxonomic assignment** in the protein arm — 58% at 5M and 69% at 25M. The variant also loses 412 and 491 pairs that production finds, so it is not a superset even at the read level, despite the screen being a union. Trimming, alignment and duplicate selection can all change which records survive.

This has a direct consequence for [section 7](#7-whether-the-new-calls-look-real): a species that appears "new" in the variant may be new reads, or it may be reads production already had, relabelled.

### Fewer reads resolve to species

| Depth | Arm | Joint exemplars | Resolved to species | Share | Above species rank or unresolved |
|---|---|---:|---:|---:|---:|
| 5M | production | 67,452 | 63,455 | 94.1% | 3,997 |
| 5M | variant | 99,206 | 55,404 | 55.8% | 43,802 |
| 25M | production | 77,223 | 76,203 | 98.7% | 1,020 |
| 25M | variant | 101,037 | 45,760 | 45.3% | 55,277 |

*Source: `taxonomic-resolution.tsv`.*

Production resolves 94–99% of its exemplars to species rank. The variant resolves 45–56%. In absolute terms it leaves 43,802 and 55,277 exemplars stuck at a broader taxonomic label, against production's 3,997 and 1,020.

This is what you would expect if protein-level matches are frequently to regions conserved across a family or order, so the evidence genuinely does not identify a species. It is also what you would expect if some of those matches are to protein regions that are not diagnostic of viruses at all. The table does not distinguish those, and it does not establish that the above-species exemplars are all newly admitted reads — shared records can move to broader labels, and [the reassignment table above](#most-shared-reads-are-relabelled-not-added) shows most of them move somewhere.

---

## 6. What we found: more species calls

### At every threshold, in both strata

| Depth | Rule | production | variant | Retained in both | Gained from zero | Promoted | Lost | Net |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 5M | N=2 | 835 | 1,325 | 726 | 522 | 77 | 109 | **+490** |
| 5M | N=3 | 661 | 1,058 | 582 | 391 | 85 | 79 | **+397** |
| 5M | N=5 | 511 | 806 | 442 | 278 | 86 | 69 | **+295** |
| 25M | N=2 | 275 | 472 | 248 | 197 | 27 | 27 | **+197** |
| 25M | N=3 | 231 | 380 | 202 | 152 | 26 | 29 | **+149** |
| 25M | N=5 | 176 | 306 | 160 | 113 | 33 | 16 | **+130** |

*Source: `threshold-calls.tsv`. "Gained from zero" means production had no species-resolved joint exemplar for that group/species at all. "Promoted" means production had support below N and the variant reached N. "Lost" means production met N and the variant did not. "Retained" does not imply unchanged support, nor the same read pairs. These count group/species combinations, not globally distinct species.*

The relative increase is 57.7–60.1% at 5M and 64.5–73.9% at 25M. It does not depend on which of the three thresholds you pick, though only three were tried and the magnitude does vary.

Note the losses. The variant loses 69–109 calls per stratum that production reports. Gains outnumber losses by roughly five to one, but the variant is not strictly finding more of everything.

### The increase is not driven by a few runs

Net change in calls per run at N=3, each `█` one sequencing run:

```
5M, 50 runs                          25M, 10 runs
 +15 │██                              +27 │█
 +14 │█                               +22 │█
 +13 │██                              +19 │█
 +12 │█████                           +18 │█
 +11 │██                              +16 │█
 +10 │█████                           +12 │█
  +9 │█████                           +11 │█
  +8 │██████                          +10 │█
  +7 │███████                          +8 │█
  +6 │█████                            +6 │█
  +5 │██                                  │
  +4 │██                              median +14
  +3 │██                              range +6 to +27
  +1 │█                               10 of 10 gain
   0 │█
  -1 │██

median +8 · range −1 to +15 · 47 gain, 1 unchanged, 2 lose
```

*Source: `calls-per-run.tsv`, `call-distribution.tsv`.*

47 of 50 runs at 5M show a net gain; the largest loss anywhere is a single call. All ten runs at 25M gain. So the aggregate difference is a property of the cohort, not of a few unusual runs. A run with no net change can still contain gained and lost species that cancel in the total. These are descriptive counts for this cohort, not a rate at which gains would occur in new data.

### Both arms gain with depth

| Arm | Rule | Calls at 5M | Calls at 25M | Added at 25M | Lost at 25M | Median net per run |
|---|---|---:|---:|---:|---:|---:|
| production | N=2 | 155 | 275 | 120 | 0 | +13.0 |
| production | N=3 | 126 | 231 | 105 | 0 | +10.0 |
| production | N=5 | 89 | 176 | 87 | 0 | +8.0 |
| variant | N=2 | 268 | 472 | 204 | 0 | +18.5 |
| variant | N=3 | 218 | 380 | 162 | 0 | +13.0 |
| variant | N=5 | 163 | 306 | 143 | 0 | +13.5 |

*Source: `depth-summary.tsv`, `tabulate_depth.py` commit `760c004f`. The ten runs measured at both depths only.*

Every one of those ten runs gains calls at greater depth in each arm and at each threshold, and no shallow call disappears at the larger depth. Duplicate selection does not guarantee that property in general, but it held throughout this cohort. The variant gains more with depth than production does.

---

## 7. Whether the new calls look real

### What the validation step actually does

After assignment, the pipeline samples output reads — up to 20 eligible records per taxonomic group within each sample group — and re-searches them with **megablast against core_nt**, a comprehensive nucleotide database far larger than the viral reference collection. Each checked row gets one of three outcomes:

- **No alignment.** megablast found nothing at all.
- **Aligned, with a taxonomic distance.** It found something, and the result is scored by how many taxonomic steps separate what the pipeline assigned from what megablast found. Distance 0 means compatible, which here allows a *directional* parent/child relationship — a parent assignment agreeing with a more specific child counts as compatible. It is not exact species agreement and it is not accuracy.
- Not sampled. Most rows; the cap samples abundant groups more sparsely than small ones.

This is the only evidence in the experiment that bears on whether an assignment is correct. It is a DNA search, which matters enormously for interpreting it.

### Split by call category

Pooling the audit over an entire arm hides the thing worth knowing. Splitting it by what the call *is* — a call both arms make, a call the variant promoted over the threshold, or a call the variant invented from nothing — separates populations that behave completely differently.

| Depth | Arm | N=3 category | Checked rows | Coverage | No alignment | | Compatible | Mean distance |
|---|---|---|---:|---:|---:|---|---:|---:|
| 5M | production | Shared calls | 7,713 | 11.9% | **0.74%** | `▏` | 49.26% | 0.998 |
| 5M | variant | Shared calls | 7,574 | 18.1% | **6.59%** | `█▌` | 53.22% | 0.978 |
| 5M | variant | Promoted | 876 | 21.1% | **41.67%** | `██████████` | 16.63% | 2.591 |
| 5M | variant | Gained from zero | 4,284 | 25.9% | **32.42%** | `████████` | **0.35%** | **6.474** |
| 25M | production | Shared calls | 2,752 | 3.1% | **0.18%** | `▏` | 49.04% | 0.993 |
| 25M | variant | Shared calls | 2,876 | 7.1% | **7.16%** | `█▋` | 51.11% | 1.107 |
| 25M | variant | Promoted | 281 | 9.2% | **40.93%** | `█████████▊` | 13.86% | 3.880 |
| 25M | variant | Gained from zero | 1,746 | 9.2% | **39.81%** | `█████████▌` | **0.29%** | **7.053** |

*Bar column scaled to 45%. Coverage is checked rows ÷ assigned rows. Compatibility and mean distance use only rows carrying a numeric taxonomic-distance annotation, so their denominators are smaller than the checked-row count. Categories include all rows assigned to the category's species, not just its joint exemplars. Source: `validation-by-category.tsv`; category method `76a198f9`, aggregation `bda00ad4`. All categories, including those omitted here, reconcile exactly to the pooled audit.*

Three readings.

**On the calls both arms make, the two metrics point in different directions.** Compatibility is similar — 51.1% against 49.0% at 25M, 53.2% against 49.3% at 5M — and mean distance is about 1.0 in every case. But the same retained calls return no alignment in 7.2% and 6.6% of checked rows, against production's 0.2% and 0.7%. We ran no equivalence test, so these are two observed metrics pointing different ways, not a finding of equivalence and not a finding of degradation. The 25M ratio rests on only 5 production events and its magnitude should be treated as unstable; the 5M stratum has 57 events and gives a ninefold ratio.

**On the calls gained from zero, compatibility essentially collapses.** Of the checked rows that got a numeric distance, 3 of 1,051 and 10 of 2,889 are compatible with their assignment. Mean distance is 6.5 to 7.1 taxonomic steps rather than about 1.0. This is the lowest compatibility of any category by a wide margin.

**But gained calls are not the worst on every metric, and they are not where most unalignable rows live.** On no-alignment rate alone, *promoted* calls are worse — 40.9% and 41.7% against 39.8% and 32.4%. And counting rows rather than rates: at 25M the variant's 2,005 no-alignment rows are 844 from rows assigned above species rank or unresolved, 695 gained-from-zero, 206 retained, 145 below threshold in both arms, and 115 promoted. At 5M the 3,939 rows are 1,389 gained-from-zero, 1,353 above-species, 499 retained, 365 promoted, 333 below-threshold. Gained and promoted together are 40.4% and 44.5%. In both strata a majority of unalignable rows sit **outside** the newly gained calls.

### Why this does not settle the question

The obvious reading is that the new calls are wrong. That reading is not available, for a specific reason.

Megablast against core_nt is a **near-exact DNA search**. Consider the two competing explanations for a new call:

- **It is a genuinely divergent virus.** Then it has no close DNA relative in core_nt either, and megablast returns no alignment. That is the expected result, not a refutation.
- **It is a non-viral read matching a protein region viruses share with other organisms.** Then megablast may well find its true source organism — some bacterium or host sequence — and report a large taxonomic distance from the viral assignment.

Both predictions follow from what the search does, not from anything measured here. For the 32–40% of rows returning no alignment, the evidence is genuinely silent: that bucket is exactly where both explanations predict the same observation.

Where megablast *does* return an alignment to a non-viral sequence, that is real evidence against the assignment. The audit is not uninformative. What it cannot do is resolve the no-alignment population, and in the protein arm that population is large.

Two further limits on this table. These are row-level statistics from sampled validation, not percentages of false calls — coverage is 9.2% and 25.9% for the gained-from-zero category. And a large taxonomic distance does not by itself establish a non-viral origin.

### The open question

About 3,900 of the variant's checked rows **did** align, six to seven taxonomic steps away from their assignment. Those rows are informative in a way the no-alignment rows are not, and nobody has looked at them.

A read with no DNA relative anywhere would land in the no-alignment bucket. These did not. So megablast found *something* for them, and core_nt is far larger than the viral reference collection. Two explanations fit:

1. A non-viral read matched a conserved protein region, and megablast then found its actual source organism.
2. A read production already had, reassigned by the protein arm — recall that 58–69% of shared pairs change assignment.

**Checking whether those alignments land under Viruses (taxid 10239) would distinguish them**, and it is cheap: the per-row data already exists in `*_validation_hits.tsv.gz`, read by `tabulate_validation_categories.py`. The files are in the preserved S3 archives rather than locally, so it needs the archives mounted and the tabulation re-run. This is the single highest-value follow-up from the experiment as it stands.

<details>
<summary><strong>The pooled audit, before the category split</strong></summary>

| Depth | Arm | Checked rows | No alignment | Rate | Rows with distance | Compatible | Rate | Mean distance |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 5M | production | 9,708 | 65 | 0.67% | 9,643 | 5,088 | 52.76% | 0.960 |
| 5M | variant | 20,796 | 3,939 | 18.94% | 16,846 | 7,812 | 46.37% | 2.205 |
| 25M | production | 3,604 | 10 | 0.28% | 3,594 | 1,972 | 54.87% | 0.920 |
| 25M | variant | 8,457 | 2,005 | 23.71% | 6,445 | 2,843 | 44.11% | 2.620 |

*Source: `validation-audit.tsv`. Unchecked rows do not inherit the sampled rows' verdicts.*

The pooled view shows the variant doing worse on both metrics, which is true but unhelpfully coarse — it averages a retained population that looks like production together with a gained population that does not.

</details>

<details>
<summary><strong>Support per call</strong></summary>

An earlier draft of this analysis claimed the variant's calls rest on thinner evidence. That was an inference from pooled means and it does not hold: the medians run in opposite directions by stratum — 15 against production's 12 at 25M, and 11 against 12 at 5M. The pooled means are lower in the variant (119.3 against 328.9 at 25M) and the maxima differ (7,661 against 21,794), but nothing here shows the maximum accounts for the gap. The defensible statement is about categories rather than arms.

| Depth | Arm, category | Calls | Median exemplars | Mean | Max |
|---|---|---:|---:|---:|---:|
| 25M | variant, retained | 202 | 20.5 | 166.6 | 7,661 |
| 25M | variant, gained from zero | 152 | 9.0 | 62.2 | 2,126 |
| 25M | production, lost | 29 | 3.0 | 4.2 | 10 |
| 5M | variant, retained | 582 | 15 | 67.3 | 4,341 |
| 5M | variant, gained from zero | 391 | 9 | 30.8 | 2,100 |
| 5M | production, lost | 79 | 4 | 4.6 | 14 |

*At N=3, exemplars restricted to called species. Source: `support-per-call.tsv`.*

At 25M, newly gained calls sit at a median of 9 exemplars against 20.5 for calls both arms make. The calls production loses have medians of 3 and 4 against a threshold of 3, so they cluster at the reporting boundary — though their maxima are 10 and 14, and exemplar counts are not a measure of biological support.

</details>

### The check that was not run

The adjudication this experiment needed is a search of the new candidates against a **comprehensive protein database** (NR). That changes the search modality, so unlike megablast it can distinguish a divergent virus from a conserved non-viral protein region.

It was specified and costed, and David authorized it within the existing budget. **Codex declined to run it**, because validating newly recovered potentially human-infecting viral candidates would extend sequence-level discovery. That is an execution limitation on the agent's side — not a result, not a budget decision, and not evidence for or against the additional calls. No recovery or search was launched.

---

## 8. What we learned, and what we did not

### What the experiment establishes

For this cohort, on identical inputs: how many more viral read records the protein pipeline emits, how many assignments change, how many thresholded calls are gained and lost at each threshold reported, and what the change cost in compute. These are measurements of these 60 runs. They are not estimates of what the change would do in production on new data.

Three findings seem solid:

1. **The output increase is large, consistent, and not an artifact of the reporting threshold.** 55–65% more calls, 47 of 50 runs gaining, stable across N=2, 3 and 5 and across both depths.
2. **The variant is a different pipeline, not a more sensitive one.** It relabels 58–69% of the reads both arms recover, loses 412–491 pairs production finds, and resolves roughly half as many reads to species rank. Describing it as "production plus extra detections" is wrong.
3. **The added calls carry evidence unlike the shared calls.** Whatever they are, they are not more of the same. That is the single most decision-relevant output of the experiment.

### What it does not establish

**Whether the pipeline is better.** There are no truth labels anywhere in this data. These are environmental samples whose real contents are unknown, so more output is not evidence of more sensitivity.

**Which change produced the effect.** The screen and the aligner were replaced together, deliberately, because [section 1](#what-we-tried-first-and-why-it-was-not-enough) showed that changing the screen alone achieved almost nothing with Bowtie2 still in place. That makes this a test of the combined change, so it cannot say whether a cheaper subset would capture most of the effect.

**How much of the effect is the acceptance criteria rather than protein space.** Bowtie2 and DIAMOND accept alignments on thresholds that were never calibrated against each other. Some share of the extra output follows from that rather than from searching in protein space, and this design cannot size that share.

**Whether the new calls are correct.** [Section 7](#why-this-does-not-settle-the-question) covers this: the available validation is a near-exact DNA search, which is expected to fail on a genuinely divergent virus and on a conserved non-viral protein region alike.

**How often viruses share protein regions with other organisms in this data.** That is the concern motivating the whole specificity question, and nothing here measures it. What can be said about the pipeline's own defences is narrow: the human and contaminant filters do run on protein-arm reads, identically to production, but they remove reads matching *specific contaminant genomes* rather than testing whether a read is viral at all, which is the question a protein-similarity hit raises.

**Anything about community abundance.** Bracken abundance is unavailable in every accepted unit, because the retained taxonomy and report combination lacks the requested rank rows. The assignment comparison completed; this is not a complete community-abundance result.

### What would settle it

In rough order of cost:

1. **Check whether the aligned-but-incompatible rows land under Viruses.** The data already exists. See [the open question](#the-open-question). Cheap, and it would convert the largest ambiguous population in the audit into a usable signal.
2. **Run the NR protein search on the gained candidates.** This is the adjudication the experiment was designed around, and it is the thing that can distinguish a divergent virus from a conserved protein domain. Not performed; see above.
3. **Build a known-divergence sensitivity panel.** Viral sequences at measured divergence from the reference at known abundance, with recall measured as a function of divergence. This is the only thing that answers the original question — does protein search find divergent viruses the DNA pipeline misses — rather than the proxy question of whether it emits more output. Nothing like it exists here and none is planned or committed.

A note on 3: it is the measurement that should probably have come first. The experiment we ran can tell you the pipeline's output changed a lot. Only a panel with known contents can tell you whether that change is the one you wanted.

---

## 9. What it cost

| Quantity | Value | Meaning |
|---|---|---|
| Compute resource-time ceiling | $131.62 | Disjoint worker intervals at conservative hourly ceilings |
| EBS accrual model | $27.57 | Six retained volumes, study window only |
| **Subtotal of those two models** | **$159.19** | Excludes the unpriced components below |
| Current private S3 inventory | 168.27 GB, 12,615 objects | Stored bytes, not a storage bill |
| **EBS if all six volumes remain** | **$967.28 per 30 days** | Charges continue now that compute has stopped |
| Experiment cap | $500 | Budget limit, not spend |
| NR secondary audit | $0 | Not performed |

*Observed 2026-10-25T04:37:24Z. Sources: `spend-supervised-0437.json`, `worker-intervals.json`, `experiments/paired_pipeline_comparison/budget.json`.*

Both workers are shut down — the primary instance stopped, the Spot instance terminated — so compute intervals are closed. **Storage is the live cost.** Six volumes remain retained and accrue about $967 per 30 days until someone decides to clean them up. Cost Explorer figures are delayed, account-wide and overlap these models, so they are not added to the subtotal; snapshot storage, S3 storage age and requests, transfer and taxes still require invoice reconciliation.

<details>
<summary><strong>Runtime by stage, and input acquisition</strong></summary>

| Stage | 5M total min. | 5M median per comparison | 25M total min. | 25M median per comparison |
|---|---:|---:|---:|---:|
| Production RUN, including shared QC and profiling | 263.1 | 4.90 | 43.1 | 4.27 |
| Additional amino-acid screen | 63.7 | 1.19 | 38.8 | 3.90 |
| Screen union | 0.4 | 0.00 | 0.2 | 0.02 |
| Variant preparation | 3.5 | 0.06 | 1.0 | 0.10 |
| Protein alignment and record handling | 174.0 | 3.21 | 55.6 | 4.73 |
| Production downstream | 320.7 | 3.70 | 36.4 | 3.68 |
| Variant downstream | 261.6 | 3.88 | 62.1 | 4.12 |

*Source: `runtime.tsv`. Workflow elapsed time is finish minus start for the whole stage. The variant reuses production's shared QC and community-profiling outputs. Production RUN precedes protein RUN, so alignment-stage order is fixed; downstream validation order is prespecified and counterbalanced. These measurements describe this execution and are not a deployment speed ratio.*

Input acquisition is shared work across both arms and both depths:

| Measurement | Runs with timing | Missing | Recorded hours | Median min. per run |
|---|---:|---:|---:|---:|
| Recorded downloads | 48 / 50 | 2 | 2.42 | 2.82 |
| Final acquisition, excluding download | 50 / 50 | 0 | 14.56 | 16.66 |
| — of which archive validation | 50 / 50 | 0 | 5.39 | 5.57 |
| — remainder, including sampling | 50 / 50 | 0 | 9.17 | 10.24 |

*Source: `acquisition-summary.tsv`, `tabulate_acquisition.py` commit `f5cc6823`. Each source run is counted once across arms and depths. Archive validation is inside acquisition time; the last two rows partition it and must not be added to it. Missing download timings are unknown, not zero. These are summed elapsed times, not worker-hours.*

The method acquired complete source archives before selecting the smaller comparison inputs, so acquisition cost is not proportional to the 5M or 25M pairs actually processed. Reusing the preserved sampled inputs avoids repeating this, but a comparison on new source runs incurs it again.

</details>

---

## Provenance

All 60 accepted comparisons have a source identifier, input checksums, execution receipts, standard downstream reports, and a private archive verified by reading it back.

<details>
<summary><strong>Records, commits and archive paths</strong></summary>

| Item | Record |
|---|---|
| Frozen selection | `selection.json`, introduced in `b3976772` |
| First comparison source | `629db5c2` |
| Subsequent scientific source | `c77cdde0` |
| Input identity and sampling | Per-comparison `input-acquisition.json`; private `inputs/<accession>/verified.json` and `backup-readback.json` |
| Viral reference collection | `reference-collection-identity.json`; same genome identifiers and lengths, with the recorded masking treatments |
| Worker hardware | AWS `r6i.12xlarge`, 48 vCPUs and 384 GiB, us-east-1 |
| Runtime and placement | Per-comparison stage receipts, Nextflow traces, `execution-worker.json`, disjoint worker-interval ledger |
| Full preserved comparison | `s3://…/comparisons/<run_depth>/completed-evidence.tar.gz` |
| Archive verification | Per-comparison `completed-evidence-backup.json`, including SHA-256 and successful readback |
| Completion inventory | `results/aggregate/runs.tsv` |
| Added species tabulation | `tabulate_species.py`, commits `32ba3d63`, `6cc6c24c`, counts-only export `22fe0d42` |
| Write-up table generation | `writeup_tables.py` through `906b4fba` |
| Matched-depth tabulation | `tabulate_depth.py`, commit `760c004f` |
| Validation by category | `tabulate_validation_categories.py` `76a198f9`; `summarize_validation_categories.py` `bda00ad4` |
| Shared acquisition timing | `tabulate_acquisition.py`, `f5cc6823` |
| Initial trace recovery | `trace-recovery/receipt.json`; recovered from the checksum-verified preserved archive |
| Final presentation tables | committed in `385c9ac8` |
| Closure records | committed in `5dd872a6` |

The count-table inputs and recovered runtime traces are preserved in S3 under SHA-256 `e95fac284c4382dba14c0762fe01419bb38de37dd5f1ab11e5153e162a19420b`, verified by readback. That archive uses an explicit whitelist and contains no candidate sequence exports.

The branch and this report remain private and unpushed. Historical experiment commits retain their own source identities; later reporting code does not replace them.

</details>

### How this document was produced

The experiment was run by Codex. The underlying lab record is `WRITEUP.md` in this directory, co-authored section by section with Claude at commit `b24c96f3`, with presentation tables at `385c9ac8` and closure receipts at `5dd872a6`. This document is a rewrite of that record for a reader who has not followed the experiment.

Each agent corrected the other's claims in writing during drafting, and `WRITEUP.md`'s changelog records every correction. Corrections that changed a stated conclusion included: a claim that a three-exemplar detection threshold existed in the pipeline (it does not); a claim that the variant's output is a superset of production's (it is not); a claim that the variant's calls rest on thinner evidence (the medians run the other way); a claim that contaminant filters do not run on protein-arm reads (they do); and a claim that the category split substitutes for the unperformed NR audit (it does not).

**Nothing in this document has been reviewed by a human.** Every number is read from the frozen tables named in each caption rather than carried forward from an earlier draft, but the reasoning around them has only been checked by the two agents that wrote it.
