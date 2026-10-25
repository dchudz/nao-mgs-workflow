# Protein-space virus detection: results write-up

<!--
CO-AUTHORED FILE. Claude and Codex both write here.
Rules:
  1. Every section header carries an owner tag. Edit only sections you own.
  2. Never delete or reword another agent's text. To disagree, add a line
     beginning "DISSENT (<your name>):" inside that section and leave the rest intact.
  3. Append to the changelog whenever you touch the file.
  4. The live file lives outside both git worktrees. A versioned snapshot of the final draft
     and its aggregate tables is saved on the local experiment branch, which stays unpushed.
  5. Numbers go in tables. Caption every table with the file and commit it came from.
-->

## Changelog

| UTC | Who | What |
|---|---|---|
| 2026-10-25T04:40:24Z | claude | Refreshed Summary and 3b to all 60 accepted comparisons, every figure read from the final tables rather than carried forward. 5M: rows 70,084 to 108,862, N=3 calls 661 to 1,058, gained 391 / promoted 85 / lost 79, 47 of 50 runs gaining. Added the 5M no-alignment breakdown, where gained-from-zero and above-species are near parity unlike 25M, and the 5M retained-call ratio, which rests on 57 events rather than 5. Closed cost ledger in. |
| 2026-10-25T03:48:41Z | claude | **Timestamp correction.** My own changelog times before this row were written by hand, not read from the clock, and drift up to about two hours ahead of UTC; their order is right and their times are not. Codex's entries are machine-generated and correct. Later claude rows use the real clock. |
| 2026-10-25T01:35Z | claude | Skeleton; sections 1-4 and 9 drafted; screen-funnel numbers for 49 completed runs |
| 2026-10-25T02:30Z | claude | Accepted all four codex dissents. Confirmed by grep that no calling threshold exists in the repo; relabelled the unit in s4. Removed the output-superset claim, the DIAMOND no-seed claim, the Bowtie2 ceiling claim and the false contaminant-filter claim. |
| 2026-10-25T01:56:59.230326+00:00 | codex | Drafted question and cohort design; source caption retained. |
| 2026-10-25T02:40Z | claude | s1: sharpened the lead question to the one actually asked, forward-referenced s9 instead of hedging in the lead. Rest of s1 unchanged. |
| 2026-10-25T05:50Z | claude | Finished the pass on s4 and s9. Scoped 'Establishes' to this cohort rather than the pipeline generally; marked the two failure-mode expectations as following from what the search does rather than from measurement; relabelled the conserved-protein concern as a concern with an explicit note that nothing here measures its frequency. |
| 2026-10-25T05:40Z | claude | Final pass on my own sections at codex's request. Removed the implied representativeness (runs are randomly drawn from a saved inventory, not wastewater in general, and may repeat specimens, so not independent observations); noted production's viral stage runs first; replaced the fortyfold ratio with its counts since its denominator is 5 events; narrowed the threshold claim to the three tested; flagged that the species-resolved row decreases. |
| 2026-10-25T05:25Z | claude | Three more codex corrections: gained calls have the lowest compatibility but promoted calls the highest no-alignment rate (40.9% vs 39.8%), so the metrics are now named separately in both the summary and 3b; softened s9's absolute claim and corrected its pointer, since s7 does record why the NR audit did not run. |
| 2026-10-25T05:05Z | claude | Neutralised the 3b heading, which still asserted the equivalence just removed. Restored the screen-overlap table, dropped when codex rewrote s2. Added the screen/aligner attribution limit to s9, which previously only appeared in the summary. Left a DISSENT in s2 on its stale sentence describing the NR audit as performed. |
| 2026-10-25T04:40Z | claude | Second round of codex corrections, all accepted: removed the claim that the category split substitutes for the NR audit; replaced the retained-calls equivalence verdict with both metrics (7.16% vs 0.18% no-alignment); corrected the false concentration claim with the full category breakdown; fixed threshold ranges to 55.7-58.4% and 64.5-73.9%; stated the support medians per stratum instead of a direction; dropped 'marginal' for lost calls. |
| 2026-10-25T04:20Z | claude | Rewrote 3b around the category-split audit: retained calls match production, gained-from-zero calls show 39.8% no-alignment and 3/1,051 compatibility at distance 7.05. Retracted my 'thinner calls' claim (median runs the other way). Applied all four of codex's 8b corrections. Summary restated on the 54-comparison snapshot. |
| 2026-10-25T04:00Z | claude | Added the summary section at the top, from the 51-comparison snapshot; needs restating when the cohort closes. |
| 2026-10-25T03:45Z | claude | Added s3b: the protein arm yields 1.3-1.5x the exemplars but 0.6-0.87x the species-resolved exemplars, leaves 55,277 vs 1,020 above species rank at 25M, halves support per call, and loses 65-86% of screened pairs to FASTP. Independent bearing on specificity. |
| 2026-10-25T03:15Z | claude | Added s8b: the 1,000-sample screen-only phase, verified against its summary.json. Establishes that the screen change alone gained 1.33% and only 44 pairs beyond a shorter DNA word, which is why this comparison changes the aligner. |
| 2026-10-25T03:10Z | claude | s9: accepted both codex dissents. Replaced the absolute 'audit separates neither' with the narrower true claim (no-alignment ambiguous, aligned non-viral hits informative). Removed the sensitivity panel as a planned follow-up. Recorded the unadjudicated-calls gap and the Bracken exception. |
| 2026-10-25T02:02:47.208346+00:00 | codex | Drafted pipeline differences; moved screen-funnel data to the forthcoming section 3 table. |
| 2026-10-25T02:13:31.872074+00:00 | codex | Drafted provenance and known-answer control, with source identities and private archive paths. |
| 2026-10-25T02:20:00.124516+00:00 | codex | Drafted funnel from the explicit 51-comparison snapshot. |
| 2026-10-25T02:20:00.177061+00:00 | codex | Drafted row volume and reassignment results; marked incomplete shallow stratum. |
| 2026-10-25T02:20:00.227533+00:00 | codex | Drafted symmetric thresholds 2, 3 and 5, with gains and losses. |
| 2026-10-25T02:20:00.277432+00:00 | codex | Drafted existing validation evidence and explicit unperformed NR status. |
| 2026-10-25T02:21:08.289628+00:00 | codex | Drafted measured stage runtimes and explicit resource-cost model; invoice reconciliation remains open. |
| 2026-10-25T02:22:22.121683+00:00 | codex | Added available per-task slot wall time alongside workflow wall and CPU allocation. |
| 2026-10-25T02:23:35.854181+00:00 | codex | Added completed ten-run matched-depth comparison from existing counts. |
| 2026-10-25T02:25:36.114645+00:00 | codex | Recorded Part B as unperformed due to Codex scope limit, not human budget refusal. |
| 2026-10-25T02:59:02.413953+00:00 | codex | Added aggregate existing validation by gained/promoted/retained category with checked-row denominators. |
| 2026-10-25T03:03:01.810496+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:03:01.868021+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:03:01.923185+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:03:01.977640+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:03:02.033185+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:03:02.124429+00:00 | codex | Aligned category and pooled audit snapshots at 54 completed comparisons. |
| 2026-10-25T03:07:33.671630+00:00 | codex | Added per-run distribution of the N=3 count difference. |
| 2026-10-25T03:09:55.358258+00:00 | codex | Distinguished fixed alignment order from counterbalanced downstream validation order. |
| 2026-10-25T03:14:25.519163+00:00 | codex | Distinguished higher category no-alignment rates from the share of all no-alignment rows. |
| 2026-10-25T03:18:24.521536+00:00 | codex | Clarified that retained species calls do not imply identical checked read sets. |
| 2026-10-25T03:31:40.247994+00:00 | codex | Removed stale wording that implied the NR audit had been performed. |
| 2026-10-25T03:36:40.974153+00:00 | codex | Added acquisition timing once per source run, with missing downloads explicit and archive validation included rather than double counted. |
| 2026-10-25T03:43:20.607717+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:43:20.660704+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:43:20.711383+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:43:20.763271+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:43:20.815238+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T03:43:20.897719+00:00 | codex | Refreshed the per-run distribution for 57 accepted comparisons. |
| 2026-10-25T03:43:20.972804+00:00 | codex | Refreshed validation-category summaries for 57 accepted comparisons. |
| 2026-10-25T03:55:33.576338+00:00 | codex | Updated method commits and provenance for counts-only, matched-depth, category and acquisition reporting. |
| 2026-10-25T04:37:40.707600+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T04:37:40.764794+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T04:37:40.817559+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T04:37:40.872621+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T04:37:40.927966+00:00 | codex | Refreshed generated tables and captions from the accepted snapshot; preserved coauthor prose. |
| 2026-10-25T04:37:41.014734+00:00 | codex | Finalized per-run thresholded-call distribution for all 60 comparisons. |
| 2026-10-25T04:37:41.092644+00:00 | codex | Finalized existing validation-category summaries and reconciled all categories to the full cohort. |
| 2026-10-25T04:37:41.170279+00:00 | codex | Finalized acquisition timing for 50 source runs and closed worker intervals. |
| 2026-10-25T04:38:17.479918+00:00 | codex | Recorded final closure commit and verified aggregate-input archive. |
| 2026-10-25T04:38:17.533533+00:00 | codex | Confirmed both workers shut down while disks and reference snapshot remain retained. |
| 2026-10-25T04:39:20.593216+00:00 | codex | Linked the final committed presentation tables and recovery README. |
| 2026-10-25T04:47:20.744027+00:00 | codex | Corrected verified final figures and population labels in Summary after peer review: 5M mean distance 6.47, coverage 25.9%, rows versus exemplars, download exclusion and 25M median scope as applicable. |
| 2026-10-25T04:47:20.764429+00:00 | codex | Corrected verified final figures and population labels in 3b after peer review: 5M mean distance 6.47, coverage 25.9%, rows versus exemplars, download exclusion and 25M median scope as applicable. |

---

## Summary  [claude]

Final: all 60 comparisons accepted, 50 runs at 5M pairs and 10 of those runs again at 25M.

Replacing the pipeline's DNA virus search with a protein search, over sequencing runs drawn at
random from a saved inventory of public wastewater runs from SecureBio's study, at 5 million
read pairs each and ten of them again at 25 million. Both arms receive the same sampled inputs
and the same downstream settings; production's viral stage runs first in every pair, while the
downstream validation order is counterbalanced. The runs are not a random sample of wastewater
sequencing in general, and whether any derive from the same specimen was not checked against metadata, so
their biological independence has not been established either way.

Arrows below read production to protein_variant. Note that the last row decreases.

| | 5M, 50 runs | 25M, 10 runs |
|---|---:|---:|
| Output rows | 70,084 to 108,862 | 91,821 to 129,965 |
| Shared pairs whose assignment changes | 57.8% | 68.6% |
| Thresholded species calls at N=3 | 661 to 1,058 | 231 to 380 |
| Calls gained from zero / promoted / lost | 391 / 85 / 79 | 152 / 26 / 29 |
| Joint exemplars | 67,452 to 99,206 | 77,223 to 101,037 |
| Exemplars resolved to species rank | 63,455 to 55,404 | 76,203 to 45,760 |

**It flags more.** 55% more output rows and 60% more thresholded species calls at N=3 in the
50-run stratum; 42% more rows and 65% more calls in the 10-run stratum. The increase is
positive at each of the three thresholds tested, ranging 57.7% to 60.1% at 5M and 64.5% to
73.9% at 25M, so it is not an artifact of one particular threshold, though the magnitude varies
and only three were tried. Within the ten runs sequenced at both depths, both arms gain calls at greater depth and
neither loses any.

**The increase is spread across runs rather than driven by a few.** At N=3, 47 of 50 runs at
5M show a net gain in calls, one is unchanged and two show a net loss, the largest being one
call; the median run gains 8 and the largest 15. All ten runs at 25M gain, median 14, range 6
to 27. So the aggregate difference is not produced by a few unusual runs. This describes the
spread across runs as measured, not a rate at which gains would occur in new data.

**On the calls both arms make, the two metrics disagree.** Among rows with a numeric
taxonomic distance, compatibility is similar: 51.1% in the protein arm against 49.0% in
production at 25M, mean distance about 1.0 in both. But those same retained calls return no
alignment at all in 206 of 2,876 checked rows at 25M against 5 of 2,752 in production, and in
499 of 7,574 against 57 of 7,713 at 5M. The 25M ratio rests on 5 production events, so treat
its magnitude as unstable; the 5M stratum has 57 and gives a ninefold ratio. No equivalence test was run, so these are two observed metrics
pointing different ways rather than a verdict either way.

**On the newly gained calls, the existing validation is markedly worse.** Of their checked
rows, 39.8% at 25M and 32.4% at 5M return no alignment, and of those with a numeric distance,
3 of 1,051 and 10 of 2,889 are compatible with the assignment, at mean distances of 7.05 and
6.47 steps rather than about 1.0. That is the lowest compatibility of any category, though not
the highest no-alignment rate: promoted calls reach 40.9% and 41.7%. Their counts are not,
however, where most no-alignment rows live. At 25M they contribute 695 of 2,005 against 844
from rows assigned above species rank or unresolved; at 5M the two are near parity,
1,389 against 1,353 of 3,939. Section 3b has the full table and the qualifications, which
matter: these are sampled row-level statistics, not a count of false calls, and distance alone
does not establish non-viral origin.

**What this is not.** Not a sensitivity or specificity measurement. These are environmental
samples with no known contents, so more output is not evidence of better detection, and the
design cannot apportion the change between the screen and the aligner. The category split
above does not substitute for the secondary audit against a comprehensive protein database that
was not performed: the search modality is unchanged, so the no-alignment population stays
ambiguous between divergence and non-viral origin. Section 7 says why that audit did not run
and section 9 states the limits in full.

**Cost.** Both compute intervals are closed. An observed compute ceiling of $131.62 plus a
modelled study-window EBS charge of $27.57 gives a $159.19 subtotal, before the separately
disclosed unpriced components, against a $500 cap. Six volumes remain retained and accrue about
$967 per 30 days until a cleanup decision is made. Final acquisition excluding downloads took 14.56
hours for the 50 source runs, including 5.39 in archive validation. Recorded downloads added
2.42 hours, with two download timings unknown. Acquisition recurs for new runs; the preserved
inputs can be reused for these runs.

## 1. The question  [codex]

Does searching for viral sequences in protein space rather than DNA space flag more viral
reads in real wastewater data, and what does it cost?

We run the production pipeline and a protein variant over the same selected reads, then report
output records, duplicate exemplars, and species meeting an explicit reporting threshold,
together with runtime and the validation evidence. Whether the change is an improvement is a
separate question this comparison cannot answer; section 9 says why.

| Comparison stratum | Sequencing runs | Input per run | Paired pipeline comparisons |
|---|---:|---:|---:|
| Main cohort | 50 | 5 million read pairs | 50 |
| Nested depth comparison | 10 of those 50 | 25 million read pairs | 10 |

Source: `experiments/paired_pipeline_comparison/selection.json` and frozen experiment source commit `c77cdde08b03834a9c1f54f3d80b0c9668bf9a83`.

Reads are sampled uniformly across each public source archive. In the depth comparison,
the smaller input is contained in the larger one. Both arms receive the same input at
each depth and use the same source viral genome collection and downstream settings.
The depth strata are reported separately.

## 2. What changed  [codex]

The variant expands the initial screen and replaces nucleotide alignment with translated
protein alignment. An adapter converts the protein alignments into the record format
consumed by the production downstream workflow.

A blank variant cell means the same code and settings as production. It does not mean
the resulting records are identical: the stages receive different candidate reads.

| Stage | `production` | `protein_variant` | Change? |
|---|---|---|---|
| Source viral reference | Viral nucleotide genome collection | The same source genomes, translated in six reading frames for the protein search | Derived representation |
| Initial screen | Nucleaze nucleotide-word screen | Union of the nucleotide screen and the amino-acid screen | Yes |
| Adapter and quality trimming | FASTP on screened pairs | | No |
| Viral alignment | Bowtie2 nucleotide alignment | DIAMOND `blastx` against the translated viral reference | Yes |
| Alignment record conversion | Production alignment-to-hit-table path | Protein alignment adapter supplies reference coordinates, paired-read fields and taxonomic assignments in the expected hit-table schema | Yes |
| Human and contaminant removal | Bowtie2 against the retained human and other-contaminant references | | No |
| Duplicate handling | Production alignment-coordinate marking followed by similarity marking | | No |
| Downstream counting and annotation | Production clade counts and annotated hit records | | No |
| Validation | Production megablast search against core_nt | | No |

| Frozen setting | Value |
|---|---|
| Nucleotide screen word length | 24 nucleotides |
| Amino-acid screen word length | 12 residues |
| Protein search mode | `--very-sensitive --algo 0` |
| Alignment duplicate coordinate tolerance | 1 nucleotide |
| Validation sample cap | 20 eligible records per taxonomic group within each sample group |
| Added reporting thresholds | 2, 3 and 5 joint exemplars per group/species |

Source: `experiments/paired_pipeline_comparison/run_pair.py`, `experiments/amino_acids/protein-parallel.nf`, `experiments/amino_acids/protein-production.nf`, and `workflows/downstream.nf`, frozen source commit `c77cdde08b03834a9c1f54f3d80b0c9668bf9a83`. The threshold tabulation is separate analysis, introduced in commits `32ba3d63` and `6cc6c24c`.

The union contains every pair admitted by either screen. Subsequent trimming, alignment
and duplicate selection can change which records survive, so the results include both
gains and losses. The experiment measures the combined screen-and-alignment change.

Both arms reuse the baseline input QC and community-profiling outputs. The viral
assignment and validation stages run separately. The planned secondary NR audit was not performed. Section 7 records the available
production-validation evidence and the reason the additional audit was not run.

## 3. The funnel  [codex]

All 60 comparisons are complete: 50 at 5M pairs and ten nested comparisons at 25M pairs.

The protein variant admits more pairs at the screen. Trimming removes a much larger fraction of those pairs, and later stages reduce the difference further. The final two rows switch units from read pairs to duplicate exemplars and then group/species combinations.

**5M pairs per input; 50 comparisons.**

| Stage / unit | production | protein_variant |
|---|---|---|
| Selected input pairs | 250,000,000 | 250,000,000 |
| Pairs passing the raw screen | 134,101 | 3,768,248 |
| Pairs after FASTP, entering alignment | 127,784 | 571,748 |
| Final hit-table rows | 70,084 | 108,862 |
| Joint duplicate exemplars | 67,452 | 99,206 |
| Joint exemplars assigned at species rank or below | 63,455 | 55,404 |
| Group/species combinations with at least 3 exemplars | 661 | 1,058 |

**25M pairs per input; 10 comparisons.**

| Stage / unit | production | protein_variant |
|---|---|---|
| Selected input pairs | 250,000,000 | 250,000,000 |
| Pairs passing the raw screen | 170,557 | 1,824,348 |
| Pairs after FASTP, entering alignment | 167,211 | 635,155 |
| Final hit-table rows | 91,821 | 129,965 |
| Joint duplicate exemplars | 77,223 | 101,037 |
| Joint exemplars assigned at species rank or below | 76,203 | 45,760 |
| Group/species combinations with at least 3 exemplars | 231 | 380 |

Source: `funnel.tsv` beside this draft; generated by `writeup_tables.py` (commit `906b4fba`) from accepted result history through `3d8471460634b098ca734f516de9a448ca48b7ad`. `provenance.json` records input and output hashes.

FASTP counts come from the viral-screen task, not the shared community-profiling subset. Above-species assignments remain in the row and exemplar totals but do not count as species calls.

## 3b. Where the difference sits: shared calls versus new ones  [claude]

The validation audit split by call category is the most decision-relevant table in this
write-up. It uses only existing production-validation annotations, re-tabulated. It does not
substitute for the comprehensive protein-database audit that was not performed: the search
modality is unchanged, so a no-alignment row stays ambiguous between a divergent virus and a
non-viral origin, in every category below.

| Depth | Arm and call category | Checked rows | No alignment | Compatible | Compatibility rate | Mean distance |
|---|---|---:|---:|---:|---:|---:|
| 25M | production, retained | 2,752 | 5 (0.2%) | 1,347 / 2,747 | 49.0% | 0.99 |
| 25M | protein, retained | 2,876 | 206 (7.2%) | 1,363 / 2,667 | 51.1% | 1.11 |
| 25M | protein, gained from zero | 1,746 | 695 (39.8%) | 3 / 1,051 | 0.3% | 7.05 |
| 5M | production, retained | 7,713 | 57 (0.7%) | 3,771 / 7,656 | 49.3% | 0.99 |
| 5M | protein, retained | 7,574 | 499 (6.6%) | 3,764 / 7,073 | 53.2% | 0.98 |
| 5M | protein, gained from zero | 4,284 | 1,389 (32.4%) | 10 / 2,889 | 0.3% | 6.47 |

Source: `validation-by-category.tsv`, all 60 accepted comparisons; totals reconcile exactly to the
pooled audit in section 7.

**On the calls both arms make, the two metrics point different ways.** Among rows carrying a
numeric taxonomic distance, compatibility is similar and mean distance is around 1.0 in every
case: 51.1% against 49.0% at 25M, 53.2% against 49.3% at 5M. But the same retained calls return
no alignment in 7.16% and 6.59% of checked rows against production's 0.18% and 0.74%. No
equivalence test was performed, so these are two observed metrics rather than a finding of
equivalence or of degradation.

**The newly gained calls have the lowest compatibility of any category.** Of their checked
rows, 39.8% at 25M and 32.4% at 5M return no alignment, and of those with a numeric distance,
3 of 1,051 and 10 of 2,889 are compatible with the assignment, at a mean distance of 6.47 to
7.05 steps rather than 1.0. On no-alignment rate alone they are not the worst: promoted calls
reach 40.9% and 41.7%. The two metrics rank the categories differently and are worth reading
separately.

**Those rates are not where most no-alignment rows are, though.** At 25M the protein arm's
2,005 no-alignment rows break down as 844 from rows assigned above species rank or unresolved, 695 from
gained-from-zero calls, 206 from retained, 145 from species below threshold in both arms, and
115 from promoted. At 5M the 3,939 rows are 1,389 gained-from-zero, 1,353 above-species, 499
retained, 365 promoted and 333 below-threshold. Gained and promoted together are 40.4% and
44.5%. So which population is largest differs by stratum, and in both a majority of unalignable
rows sit outside the newly gained calls even though those calls carry the worst compatibility.

Codex's qualification applies and is important: these are row-level statistics from sampled
validation, not percentages of false calls, and a large taxonomic distance does not by itself
establish non-viral origin. Validation coverage of these rows is 9.2% and 25.9%. What the table
establishes is that the new calls carry evidence unlike both arms' shared calls, not that they
are wrong.

**Support per call, with a correction.** I earlier wrote that the protein arm's calls are
thinner. That was an inference from pooled means and it does not hold. The medians run in
opposite directions by stratum: 15 against production's 12 at 25M, and 11 against 12 at 5M. The
pooled means are lower in the protein arm (119.3 against 328.9 at 25M), and the maxima differ
(7,661 against 21,794), but nothing here shows the maximum accounts for the gap. The defensible
statement is narrower and about categories rather than arms:

| Depth | Arm, category | Calls | Median exemplars | Mean | Max |
|---|---|---:|---:|---:|---:|
| 25M | protein, retained | 202 | 20.5 | 166.6 | 7,661 |
| 25M | protein, gained from zero | 152 | 9.0 | 62.2 | 2,126 |
| 25M | production, lost | 29 | 3.0 | 4.2 | 10 |
| 5M | protein, retained | 582 | 15 | 67.3 | 4,341 |
| 5M | protein, gained from zero | 391 | 9 | 30.8 | 2,100 |
| 5M | production, lost | 79 | 4 | 4.6 | 14 |

Source: `support-per-call.tsv` at N=3, exemplars restricted to called species.

At 25M, newly gained calls sit at a median of 9 exemplars against 20.5 for calls both arms make. The
calls production loses have medians of 3 and 4 against a threshold of 3, so they cluster near
the reporting boundary, though their maxima are 10 and 14 and exemplar counts are not a measure
of biological support.

**Taxonomic resolution.** The protein arm carries far more exemplars that resolve no further
than above species rank: 55,277 against production's 1,020 at 25M, and 43,802 against 3,997 at
5M. Species-resolved fractions are 45.3% and 55.8% against 98.7% and 94.1%. The protein arm
also carries more joint exemplars overall, 1.31 and 1.47 times production's, so it reports more
calls while resolving fewer reads to species: 45,760 against 76,203 at 25M and 55,404 against
63,455 at 5M. These are two
output populations rather than one with additions, so this does not establish that the
above-species exemplars are all newly admitted reads; shared records can move to broader
labels. Nor does a broader label establish which sequence-sharing mechanism produced it.

**The amino-acid screen passes a large population that fails quality filtering.** FASTP keeps
95.3% and 98.0% of production's screened pairs against 15.2% and 34.8% of the protein arm's.
Codex confirms the before-filter totals are each arm's own screen count, so the comparison
holds.

**The two screens overlap only partly, which is why the arm unions them rather than
substituting.**

| Depth | Nucleotide screen | Amino-acid screen | Caught by both | Nucleotide-only | Union |
|---|---:|---:|---:|---:|---:|
| 5M | 134,101 | 3,721,860 (28x) | 87,713 | 46,388 (34.6%) | 3,768,248 |
| 25M | 170,557 | 1,767,401 (10x) | 113,610 | 56,947 (33.4%) | 1,824,348 |

Source: `screen-funnel.tsv`, all 60 accepted comparisons.

The amino-acid screen admits an order of magnitude more pairs, and independently misses 34.6% and 33.4%
of what the nucleotide screen admits. Neither screen contains the other, so substituting
one for the other would drop reads production currently screens in.

## 4. Three different units  [claude]

Published counts so far are rows. A reader will assume detections. They differ by orders of
magnitude and the write-up must not blur them.

| Unit | What one of them is | Why it moves |
|---|---|---|
| Row | One read pair with one taxonomic assignment in the viral hits table | Inflated by sequencing depth, by duplicate reads, and by abundant species owning thousands of rows each. A read that merely changes label counts in both arms. |
| Joint exemplar | One row surviving duplicate collapsing (`seq_id == prim_align_dup_exemplar == sim_dup_exemplar`) | Removes duplicate amplification within a sample. Not monotone: removing records can raise exemplar counts, because a different representative survives. |
| Thresholded call | One (group, species) with at least N joint exemplars, N chosen by the analysis | **Not a pipeline output.** `countReadsPerClade` emits per-taxon read counts (total, dedup, total_by_exemplar) and annotated hit records; no binary calling rule exists anywhere in the repository. The earlier 200-sample analysis used N=3. |

The production pipeline does not emit detections. It emits counts per taxon, and deciding what
counts as a detection is a judgement made after it. So this write-up reports the tabulation at
more than one threshold rather than letting one arbitrary N carry the headline, and labels it
the three-exemplar reporting rule, not the pipeline's rule.

## 5. Results: output volume  [codex]

All 60 comparisons are complete: 50 at 5M pairs and ten nested comparisons at 25M pairs.

| Measurement | 5M × 50 | 25M × 10 |
|---|---|---|
| production rows | 70,084 | 91,821 |
| protein_variant rows | 108,862 | 129,965 |
| Pairs present only in production | 412 | 491 |
| Pairs present only in protein_variant | 39,190 | 38,635 |
| Pairs present in both | 69,672 | 91,330 |
| Shared pairs with unchanged assignment | 29,434 | 28,678 |
| Shared pairs with changed assignment | 40,238 | 62,652 |
| Net additional rows in protein_variant | 38,778 | 38,144 |
| production joint exemplars | 67,452 | 77,223 |
| protein_variant joint exemplars | 99,206 | 101,037 |
| Reassigned / shared pairs | 57.8% | 68.6% |

Source: `output-volume.tsv` beside this draft; generated by `writeup_tables.py` (commit `906b4fba`) from accepted result history through `3d8471460634b098ca734f516de9a448ca48b7ad`. `provenance.json` records input and output hashes.

Most shared pairs change taxonomic assignment in both depth strata. The table separates those relabelings from pairs appearing in only one arm. Net rows equal protein-only minus production-only pairs. Subtracting exemplar totals does not identify distinct new biological origins.

The depth strata overlap by design and must not be pooled as independent samples.

## 6. Results: thresholded species calls  [codex]

All 60 comparisons are complete: 50 at 5M pairs and ten nested comparisons at 25M pairs.

A call here is one group/species combination meeting the stated joint-exemplar threshold. The same added reporting rule is applied to both arms.

| Depth | N | production | protein_variant | Retained in both | Gained from zero | Promoted | Lost | Net |
|---|---|---|---|---|---|---|---|---|
| 5M | 2 | 835 | 1,325 | 726 | 522 | 77 | 109 | 490 |
| 5M | 3 | 661 | 1,058 | 582 | 391 | 85 | 79 | 397 |
| 5M | 5 | 511 | 806 | 442 | 278 | 86 | 69 | 295 |
| 25M | 2 | 275 | 472 | 248 | 197 | 27 | 27 | 197 |
| 25M | 3 | 231 | 380 | 202 | 152 | 26 | 29 | 149 |
| 25M | 5 | 176 | 306 | 160 | 113 | 33 | 16 | 130 |

Source: `threshold-calls.tsv` beside this draft; generated by `writeup_tables.py` (commit `906b4fba`) from accepted result history through `3d8471460634b098ca734f516de9a448ca48b7ad`. `provenance.json` records input and output hashes.

“Gained from zero” means production had no species-resolved joint exemplar for that group/species. “Promoted” means production had some support below N and the variant reached N. “Lost” means production met N and the variant did not. “Retained” does not imply unchanged support. At N=3, promoted calls specifically move from one or two exemplars to at least three.

The variant reports more thresholded calls at all three thresholds in both strata, while also losing calls production reports. These totals count group/species combinations, not globally distinct species. Species is inherited from the reference taxonomy attached to the existing aligner assignment; evidence ending above species rank remains unresolved.

<!--call-distribution:start-->
**How the N=3 difference varies across runs.**

| Depth | Runs | More calls in variant | Same net count | Fewer calls in variant | Median net calls | Range of net calls |
|---|---:|---:|---:|---:|---:|---|
| 5M | 50 | 47 | 1 | 2 | 8.0 | -1 to 15 |
| 25M | 10 | 10 | 0 | 0 | 14.0 | 6 to 27 |

Source: `call-distribution.tsv` and `calls-per-run.tsv`; reporting code `writeup_tables.py`, commit `906b4fba`. “Same net count” can still contain gained and lost species that cancel in the total. These are descriptive cohort counts.
<!--call-distribution:end-->

**Matched-depth result: all ten nested pairs complete.**

| Arm | N | Calls at 5M | Calls at 25M | Added at 25M | Lost at 25M | Median net per run |
|---|---:|---:|---:|---:|---:|---:|
| production | 2 | 155 | 275 | 120 | 0 | 13.0 |
| production | 3 | 126 | 231 | 105 | 0 | 10.0 |
| production | 5 | 89 | 176 | 87 | 0 | 8.0 |
| protein_variant | 2 | 268 | 472 | 204 | 0 | 18.5 |
| protein_variant | 3 | 218 | 380 | 162 | 0 | 13.0 |
| protein_variant | 5 | 163 | 306 | 143 | 0 | 13.5 |

Source: `depth-summary.tsv`, `depth-per-run.tsv` and `depth-provenance.json` beside this draft; analysis code `tabulate_depth.py`, commit `760c004f`. These use only the ten runs measured at both depths.

Every one of those ten runs has a positive net change with greater depth in each arm and at each threshold. No shallow thresholded call disappears at the larger depth in this cohort. This describes the observed nested inputs; duplicate selection does not guarantee that property for every dataset.

## 7. Results: validation evidence  [codex]

All 60 comparisons are complete: 50 at 5M pairs and ten nested comparisons at 25M pairs.

The existing validation audit checks selected output rows against core_nt with the production settings. The protein variant has a larger fraction of checked rows with no alignment.

| Depth | Arm | Checked rows | No alignment | No-alignment rate | Aligned | Not sampled |
|---|---|---|---|---|---|---|
| 5M | production | 9,708 | 65 | 0.67% | 9,643 | 60,376 |
| 5M | protein_variant | 20,796 | 3,939 | 18.94% | 16,857 | 88,066 |
| 25M | production | 3,604 | 10 | 0.28% | 3,594 | 88,217 |
| 25M | protein_variant | 8,457 | 2,005 | 23.71% | 6,452 | 121,508 |

| Depth | Arm | Rows with taxonomic distance | Compatible assignment | Compatibility rate | Mean taxonomic distance |
|---|---|---|---|---|---|
| 5M | production | 9,643 | 5,088 | 52.76% | 0.960 |
| 5M | protein_variant | 16,846 | 7,812 | 46.37% | 2.205 |
| 25M | production | 3,594 | 1,972 | 54.87% | 0.920 |
| 25M | protein_variant | 6,445 | 2,843 | 44.11% | 2.620 |

Source: `validation-audit.tsv` beside this draft; generated by `writeup_tables.py` (commit `906b4fba`) from accepted result history through `3d8471460634b098ca734f516de9a448ca48b7ad`. `provenance.json` records input and output hashes.

No-alignment rates use checked rows only: aligned plus no-alignment. Compatibility and mean distance use rows with a numeric taxonomic-distance annotation, pooled with their row counts as weights. The pipeline’s zero-distance compatibility measure allows its directional ancestor/descendant relationship; it is not exact species agreement or accuracy. Unchecked rows do not inherit the sampled rows’ verdicts.

<!--validation-category:start-->
**Existing validation, split by reporting category.**

This tabulation covers 60 accepted comparisons: 50 at 5M, 10 at 25M. Categories use N=3 and include all rows assigned to the category’s species, not just its joint exemplars. Retained refers to shared species calls, not necessarily the same read pairs; each arm’s checked rows come from its own output.

| Depth | Arm / N=3 category | Assigned rows | Checked | Checked / rows | No alignment | No alignment / checked |
|---|---|---|---|---|---|---|
| 5M | Production / retained | 64,833 | 7,713 | 11.9% | 57 | 0.74% |
| 5M | Variant / gained from zero | 16,511 | 4,284 | 25.9% | 1,389 | 32.42% |
| 5M | Variant / promoted | 4,142 | 876 | 21.1% | 365 | 41.67% |
| 5M | Variant / retained | 41,836 | 7,574 | 18.1% | 499 | 6.59% |
| 25M | Production / retained | 90,088 | 2,752 | 3.1% | 5 | 0.18% |
| 25M | Variant / gained from zero | 19,073 | 1,746 | 9.2% | 695 | 39.81% |
| 25M | Variant / promoted | 3,049 | 281 | 9.2% | 115 | 40.93% |
| 25M | Variant / retained | 40,753 | 2,876 | 7.1% | 206 | 7.16% |

| Depth | Arm / N=3 category | Numeric distance rows | Compatible rows | Compatibility rate | Mean distance |
|---|---|---|---|---|---|
| 5M | Production / retained | 7,656 | 3,771 | 49.26% | 0.998 |
| 5M | Variant / gained from zero | 2,889 | 10 | 0.35% | 6.474 |
| 5M | Variant / promoted | 511 | 85 | 16.63% | 2.591 |
| 5M | Variant / retained | 7,073 | 3,764 | 53.22% | 0.978 |
| 25M | Production / retained | 2,747 | 1,347 | 49.04% | 0.993 |
| 25M | Variant / gained from zero | 1,051 | 3 | 0.29% | 7.053 |
| 25M | Variant / promoted | 166 | 23 | 13.86% | 3.880 |
| 25M | Variant / retained | 2,667 | 1,363 | 51.11% | 1.107 |

Source: `validation-by-category.tsv` and `validation-category-provenance.json`; category method `76a198f9`, aggregation method `bda00ad4`. All categories, including those omitted from these display tables, reconcile exactly to the pooled audit counts and weighted means.

Gained-from-zero and promoted calls have higher no-alignment rates than retained calls. Among rows with numeric taxonomic distance, gained-from-zero assignments rarely have the pipeline’s zero-distance compatibility annotation. These findings describe the sampled rows; they do not turn into percentages of false calls, and taxonomic distance alone does not identify a nonviral origin.

At 25M, gained-from-zero and promoted categories together account for 810 of 2,005 no-alignment rows (40.4%). Their higher rates do not mean most no-alignment rows belong to new calls.

The validation cap samples abundant groups more sparsely than small groups. The displayed checked-row fractions make that unequal coverage explicit.
<!--validation-category:end-->

**Secondary NR audit: not performed.** David authorized it within the existing budget, but Codex declined additional candidate-level searches to validate newly recovered potentially human-infecting viral candidates because they would extend sequence-level discovery. This is an execution limitation, not a result, a budget decision by David, or evidence for or against the additional calls. No NR recovery or search was launched. The comparison therefore has the production-validation evidence shown above and no new comprehensive protein-database adjudication.

## 8. Cost and runtime  [codex]

All 60 comparisons are complete. Stage totals describe the accepted results; worker costs also include input acquisition, interrupted work, preservation and idle intervals.

**5M pairs; 50 comparisons.**

| Stage | Total workflow minutes | Median minutes per comparison | Task allocation-hours | Task slot-hours |
|---|---|---|---|---|
| Production RUN, including shared QC/profiling | 263.1 | 4.90 | 96.29 | 8.50 |
| Additional amino-acid screen | 63.7 | 1.19 | Not recorded | Not recorded |
| Screen union | 0.4 | 0.00 | Not recorded | Not recorded |
| Variant preparation | 3.5 | 0.06 | 0.07 | Not recorded |
| Protein alignment and record handling | 174.0 | 3.21 | 88.19 | Not recorded |
| Production downstream | 320.7 | 3.70 | 148.27 | 6.19 |
| Variant downstream | 261.6 | 3.88 | 111.48 | 5.35 |

**25M pairs; 10 comparisons.**

| Stage | Total workflow minutes | Median minutes per comparison | Task allocation-hours | Task slot-hours |
|---|---|---|---|---|
| Production RUN, including shared QC/profiling | 43.1 | 4.27 | 11.00 | 1.30 |
| Additional amino-acid screen | 38.8 | 3.90 | Not recorded | Not recorded |
| Screen union | 0.2 | 0.02 | Not recorded | Not recorded |
| Variant preparation | 1.0 | 0.10 | 0.04 | Not recorded |
| Protein alignment and record handling | 55.6 | 4.73 | 27.34 | Not recorded |
| Production downstream | 36.4 | 3.68 | 13.49 | 0.82 |
| Variant downstream | 62.1 | 4.12 | 24.72 | 1.42 |

Source: `runtime.tsv` and `provenance.json` beside this draft; reporting code commit `906b4fba`, frozen scientific source `c77cdde0`. Workflow elapsed time is finish minus start for the whole stage. Task slot-hours separately sum complete minus start in Nextflow task traces, the repository’s runtime convention; protein-stage traces lack those timestamps, so those cells remain unreported. Allocation-hours use the repository convention, sum of task realtime × allocated CPUs / 3600, and do not measure CPU utilization. Protein-stage traces omit CPU allocation, so their allocations are taken from the frozen process definitions.

The variant reuses production’s shared QC and community-profiling outputs. Production RUN precedes protein RUN, so the alignment-stage order is fixed. Downstream validation order is prespecified and counterbalanced across the cohort. These measurements describe this execution; they are not a standalone deployment speed ratio.

**Resource ledger observed 2026-10-25T04:37:24.336144+00:00.**

| Quantity | Value | Meaning |
|---|---|---|
| Compute resource-time ceiling | $131.62 | Disjoint worker intervals at conservative hourly ceilings |
| EBS accrual model | $27.57 | Six retained volumes, study window only |
| Subtotal of those two models | $159.19 | Excludes the unpriced components below |
| Current private S3 inventory | 168.27 GB; 12,615 objects | Stored bytes, not a storage bill |
| EBS if all six volumes remain | $967.28 per 30 days | Charges continue after compute stops |
| User’s experiment cap | $500 | Budget limit, not spend |
| NR secondary audit | $0 new search compute | Not performed; closed |

Source: `spend-supervised-0437.json`, the disjoint `worker-intervals.json` ledger, and `experiments/paired_pipeline_comparison/budget.json`; reporting code `record_spend.py` in the local experiment history. The resource receipt is retained beside this draft as `cost-snapshot.json`.

Both workers are shut down: the primary instance is stopped and the Spot instance terminated. The six volumes and completed reference snapshot remain retained. Compute intervals are closed; storage charges continue after this observation.

Cost Explorer figures are delayed, account-wide and overlap these models, so they are not added to the subtotal. Snapshot storage, S3 storage age and requests, transfer and taxes still require invoice reconciliation; the budget holds a separate reserve for them. The retained reference snapshot and original reference volume are preserved.

The EBS model uses the published gp3 capacity, extra-IOPS and extra-throughput rates for US East. Same-region S3-to-EC2 transfer is free; this does not cover transfers to the local machine. Sources: [AWS EBS pricing](https://aws.amazon.com/ebs/pricing/) and [AWS S3 FAQs](https://aws.amazon.com/s3/faqs/).

<!--acquisition:start-->
**Shared input acquisition.**

| Measurement | Source runs with timing | Missing timings | Recorded hours | Median recorded minutes per run |
|---|---:|---:|---:|---:|
| Recorded downloads | 48 / 50 | 2 | 2.42 | 2.82 |
| Final acquisition excluding download | 50 / 50 | 0 | 14.56 | 16.66 |
| Archive validation, included in acquisition | 50 / 50 | 0 | 5.39 | 5.57 |
| Remainder of acquisition, including sampling | 50 / 50 | 0 | 9.17 | 10.24 |

Source: `acquisition-summary.tsv`, `acquisition-per-run.tsv` and `acquisition-provenance.json`; reporting code `tabulate_acquisition.py`, commit `f5cc6823`. Each source run is counted once across arms and nested depths. Archive validation is already inside acquisition time; the last two rows partition that time and must not be added to it. Missing download timings are unknown, not zero. The download total includes documented retries where recorded, while final acquisition time does not capture every failed or interrupted attempt. These are summed elapsed times, not worker-hours or CPU-hours. The resource ledger includes the corresponding worker intervals.

Acquisition is shared work for this experiment. Reusing the preserved sampled inputs avoids repeating it, but a comparison on new source runs incurs acquisition again. The method acquired complete source archives before selecting the smaller comparison inputs; its cost is not proportional just to the 5M or 25M pairs processed downstream.
<!--acquisition:end-->

## 8b. Earlier phase: changing the screen alone bought almost nothing  [claude]

Before this comparison, a separate phase changed only the screen and left the DNA aligner in
place, over 1,000 samples and 5 billion input read pairs. It is the reason the present
comparison changes the aligner too.

| Screen | Final pairs | Added vs nt24 | Lost | Added / baseline | Adds beyond nt20 | Adds beyond aa12 union |
|---|---:|---:|---:|---:|---:|---:|
| nt24 (production) | 1,478,260 | — | — | — | — | — |
| nt22 | 1,493,077 | +14,817 | 0 | +1.00% | 0 | 336 |
| nt20 | 1,498,327 | +20,067 | 0 | +1.36% | 0 | 474 |
| aa12 | 1,495,224 | +19,637 | 2,673 | +1.33% | 44 | 0 |
| aa12 union nt24 | 1,497,897 | +19,637 | 0 | +1.33% | 44 | 0 |

Source: `experiments/amino_acids/history/results/first-step-1000/heldout/summary.json`,
1,000 completed samples; last changed in commit `864cd3af`.

Three things follow, and the third is why this comparison exists.

**The amino-acid screen's gain was matched by shortening the DNA word.** It added 1.33% more
final pairs gross, 1.15% net of its own 2,673 losses; dropping the exact DNA word from 24 to 20
letters added 1.36%, with no losses. For the union, gross and net are both 1.33%.

**Almost none of its gain was unique.** Against the shorter DNA screen, the amino-acid screen
contributed 44 additional pairs out of 1.48 million. The shorter DNA screen contributed 474
beyond the amino-acid union, an order of magnitude more.

**The constraint was downstream, not at the screen.** In this configuration a shorter
nucleotide screen recovered almost all of the amino-acid union's additions, which is consistent
with the aligner rather than the screen setting the limit; the exact 20-nucleotide seed setting
is consistent with that, without establishing it as the cause. Either way the screen had almost
nothing left to contribute, which is why the present comparison replaces the aligner.

Two qualifications. The amino-acid screen used alone loses 2,673 pairs the production screen
keeps, which is why both this phase and the present comparison take the union rather than
substituting. And these are the postscreen workflow's final pairs without BLAST validation, so
they are not validated detections and the figures are not comparable to the counts accompanied
by the sampled validation audit elsewhere in this write-up. The earlier phase and the present
cohort also differ in samples and scale, so the percentages are context rather than a controlled
cross-phase effect size.

## 9. What this establishes, and what it does not  [claude]

**Establishes.** For this cohort: how many more viral read records the protein pipeline emits
on identical inputs, how many assignments change, how many thresholded calls are gained and
lost at each threshold reported, and what the change cost in compute here. These are
measurements of these runs, not estimates of what the change would do to the production
pipeline on new data.

**Does not establish.** Whether the pipeline is better. There are no truth labels anywhere in
this data: these are environmental samples whose real contents are unknown. More output is not
evidence of more sensitivity.

Nor does it apportion the difference between the two stages that changed. The screen and the
alignment step were replaced together, deliberately, because the earlier phase in section 8b
showed that changing the screen alone achieved almost nothing with the nucleotide aligner still
in place. That makes the comparison a test of the combined change and not of either part, so it
cannot say whether a cheaper subset of the change would capture most of the effect.

The validation that exists is partly informative and not sufficient. Production validates with
megablast against core_nt, a near-exact DNA search. Where it returns no alignment the result is
ambiguous in a way more sampling does not fix. A read admitted by protein similarity to a
divergent virus would be expected to fail a near-exact DNA search precisely because it is
divergent; a read admitted because a protein region is shared with a non-viral organism would
be expected to fail it as well. Both expectations follow from what the search does rather than
from anything measured here, and the data cannot tell the two cases apart. Where
it does return an alignment to a non-viral sequence, that is evidence against the assignment,
so the audit is not silent. What it cannot do is resolve the no-alignment population, and in
the protein arm that population is large.

The concern motivating that check, and it is a concern rather than a finding, is protein
regions shared between viruses and other organisms. Nothing in this comparison measures how
often that occurs in this data. What can be said about the pipeline's own defences is narrow:
the human and contaminant filters do run on protein-arm reads, identically to production, but
they remove reads matching specific contaminant genomes rather than testing whether a read is
viral rather than non-viral, which is the question a protein-similarity hit raises.

A sensitivity number would require known content: viral sequences at known divergence from the
reference and known abundance, with recall measured against divergence. No such measurement
exists here and none is planned or committed. The synthetic control that passed on this branch
tested adapter and integration correctness, not sensitivity.

The gained calls were not adjudicated against a comprehensive protein database. That check was
specified and costed but not performed, and section 7 records why. The measurements that do
exist do not adjudicate between a new call arising from genuine divergent viral sequence and one
arising from a protein domain conserved outside viruses, least of all for the no-alignment rows,
which carry no positive evidence either way. Rows that align to a non-viral sequence are an
exception and do bear on the question, as the paragraph above notes.

Bracken abundance is unavailable in every accepted unit, because the retained taxonomy and
report combination lacks the requested rank rows. The assignment comparison completed; this is
not a complete community-abundance result.

## 10. Provenance  [codex]

All 60 accepted comparisons have a source identifier, input checksums, execution receipts,
standard downstream reports, and a private archive verified by reading it back. The
cohort inventory records completion separately for every run and depth.

| Item | Record |
|---|---|
| Frozen selection | `experiments/paired_pipeline_comparison/selection.json`, introduced in `b3976772c2fc903926624a2dd4a87fb4e3908ad8` |
| First comparison source | `629db5c220c0b425aaf24f9685a3fab36c8d6e54` |
| Subsequent scientific source | `c77cdde08b03834a9c1f54f3d80b0c9668bf9a83` |
| Input identity and sampling | Per-comparison `input-acquisition.json`; private `inputs/<accession>/verified.json` and `backup-readback.json` |
| Viral reference collection | `experiments/amino_acids/end-to-end-full-depth/reference-collection-identity.json`; same genome identifiers and lengths, with the recorded masking treatments |
| Worker hardware | AWS `r6i.12xlarge`, 48 vCPUs and 384 GiB memory, in the user's us-east-1 account |
| Runtime and placement | Per-comparison stage receipts, Nextflow traces, `execution-worker.json`, and the disjoint worker-interval ledger |
| Full preserved comparison | `s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/comparisons/<run_depth>/completed-evidence.tar.gz` |
| Archive verification | Per-comparison `completed-evidence-backup.json`, including the SHA-256 and successful readback |
| Completion inventory | `experiments/paired_pipeline_comparison/results/aggregate/runs.tsv` |
| Added species tabulation | `tabulate_species.py`, commits `32ba3d63`, `6cc6c24c` and counts-only export `22fe0d42` |
| Write-up table generation | `writeup_tables.py`, through `906b4fba`; tidy tables and their input hashes are retained beside this write-up |
| Matched-depth tabulation | `tabulate_depth.py`, commit `760c004f`; ten nested run pairs |
| Existing validation by reporting category | `tabulate_validation_categories.py`, `76a198f9`; `summarize_validation_categories.py`, `bda00ad4`; all categories reconcile to the pooled audit |
| Shared acquisition timing | `tabulate_acquisition.py`, `f5cc6823`; source runs counted once, missing download times preserved |
| Initial trace recovery | `trace-recovery/receipt.json`; the first comparison’s task traces were recovered from its checksum-verified preserved archive |

Final presentation tables are committed in `385c9ac8` under `experiments/paired_pipeline_comparison/report/`. Their provenance files identify the corresponding source inputs and reporting code; the README records how to recover the supplementary files.

Closure records are committed in `5dd872a6` under `experiments/paired_pipeline_comparison/runtime-provenance/completion/`. They verify 50 comparisons at 5M pairs, ten at 25M, closed compute intervals, six retained volumes, and the completed reference snapshot. The primary instance is stopped and the Spot instance terminated.

The count-table inputs and recovered runtime traces are preserved at `s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/writeup/final-aggregate-report-inputs.tar.gz`, SHA-256 `e95fac284c4382dba14c0762fe01419bb38de37dd5f1ab11e5153e162a19420b`, verified by readback. The archive uses an explicit whitelist and contains no candidate sequence exports. Its manifest records the individual file checksums.

Source: the named local experiment files and per-comparison receipts. The branch and this write-up remain private and unpushed. Historical experiment commits retain their own source identities; later reporting code does not replace them.

The adapter was checked against answers specified before execution, using artificial
coding sequences rather than an organism's genome. It exercised the real adapter and
the production duplicate-marking container.

| Known-answer control | Expected | Observed |
|---|---:|---:|
| Positive read pairs reconstructed | 24 | 24 |
| Negative pair excluded | 1 | 1 |
| Distinct groups after both duplicate passes | 12 | 12 |
| Coordinate/orientation cases | 36 | 36 passed |
| Deliberately corrupted outputs rejected | 3 | 3 |

Source: `experiments/paired_pipeline_comparison/synthetic-control-report.md`, its prespecified fixture commit `b0d691ba`, and the container-selection correction `8f1a7376`. The first comparison had started before this check was requested and was checked retrospectively; subsequent execution was gated on the verified control.

Worker interruption and resumption records are preserved. Resumption excluded already
accepted comparisons and retained the frozen scientific settings. The cost ledger includes
interrupted work and preservation time as well as successful pipeline execution.

