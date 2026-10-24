# What the paired comparison measures

The experimental unit is one selected sequencing run at one fixed processing
depth. Both arms receive exactly the same selected paired reads and qualities.
The 5-million-pair subset is nested inside the 25-million-pair subset when both
are run. These are subsets of complete archives, not full-depth executions.
The frozen design has 50 runs at 5 million pairs and 10 of those at 25 million;
the initial feasibility stage has 10 runs and two deeper cases. Execution
beyond that stage remains conditional on measured cost and complete outputs.

The 50 selected runs have 50 distinct BioSample identifiers. Collection date
and country were absent for 26 runs in the frozen ENA metadata. A subsequent
NCBI BioSample audit recovered dates and location labels for all 50, without
changing selection. There are 23 location labels and 47 collection dates,
spanning 2024-10-06 through 2026-06-17; no location/date combination repeats.
Location labels vary in granularity and do not establish 23 independent sites.
Repeated locations and temporal relationships remain relevant. The results
describe the selected run collection and do not estimate prevalence. The raw
public metadata, its checksum and the accession mapping are retained in
`runtime-provenance/biosample-metadata-audit.json` and `biosample-metadata.xml`.

| Measurement | Supports | Does not establish |
| --- | --- | --- |
| Extra or missing output pairs | A difference in pipeline output under the frozen policies | True positives, false negatives, or improved sensitivity |
| Changed assignments on shared pairs | Different taxonomic decisions for the same observations | Which decision is correct; a change can also be a change in resolution |
| Production validation status | What the same production validation policy reports for each arm | A paired accuracy comparison: the sampled records and candidate sets can differ |
| Joint exemplar counts | The production duplicate rules applied to each complete arm | Independent biological events, or distinct exemplars attributable to additions by subtracting totals |
| Nested depth differences | What changes when more pairs from the same run are processed | A monotone increase in calls, or independent observations at the two depths |
| Synthetic known-origin control | The declared artificial coordinate, label, and duplicate cases pass | Complete biological accuracy or coverage of every possible adapter input |
| Stage wall time | Elapsed time for the recorded stage on the recorded worker | Standalone deployment cost, actual CPU use, or a speed advantage from the first order-confounded case |
| Worker-time accounting | A compute spending bound from launch time and rate ceilings | An invoice, or inclusion of persistent storage and transfer charges |

The variant changes both its upstream screen and alignment representation.
It uses the same source genome collection, downstream taxonomy, duplicate
procedures, validation cap and validation mode. This is a comparison of two
frozen pipelines; it cannot isolate alignment as the sole cause of a change.
There is no NR agreement requirement or additional exemplar calling threshold.

Profiling and raw-input QC are reused between arms. Consequently, their
identity is expected and is not independent evidence of reproducibility. Their
wall time is paid once in this experiment. It must not be omitted from a future
standalone deployment estimate or charged twice in the experiment ledger.

The first run has an explicitly unavailable shared Bracken abundance profile.
The raw profiling task completed, but the frozen taxonomy/report rank mismatch
caused the production wrapper to emit an empty output. Viral assignment and
validation consume separate inputs. Assignment comparisons can therefore be
reported with that limitation, but full-schema success and abundance conclusions
cannot be claimed. The same exception is verified separately for every run.

Completed units enter the aggregate only after both arms, output checks, and
checksum-verified private preservation are present. Missing or blocked units
stay in the inventory rather than becoming zero results. Shallow and deep
strata are reported separately.
