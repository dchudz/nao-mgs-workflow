"""Build a private, self-contained report from saved scientific result tables."""
import base64, datetime, html, json
from pathlib import Path

r = Path(__file__).resolve().parent
results = r / 'results'
sections = []

def load(path):
    p = results / path
    return json.loads(p.read_text()) if p.exists() else None

def fmt(x):
    if x is None:
        return 'Pending'
    if isinstance(x, int):
        return f'{x:,}'
    if isinstance(x, float):
        return f'{x:,.2f}'
    return html.escape(str(x))

def table(headers, rows):
    return '<div class="table"><table><thead><tr>' + ''.join(f'<th>{html.escape(h)}</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(f'<td>{fmt(x)}</td>' for x in row) + '</tr>' for row in rows) + '</tbody></table></div>'

def section(title, body):
    sections.append(f'<section><h2>{html.escape(title)}</h2>{body}</section>')

def image(name, caption):
    p = r / 'figures' / f'{name}.png'
    if not p.exists():
        return ''
    return f'<figure><img alt="{html.escape(caption)}" src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}"><figcaption>{html.escape(caption)}</figcaption></figure>'

native = load('heldout/summary.json')
protein = load('protein200/summary.json')
full = load('protein-full/native-comparison.tsv.json')
repeat = load('repeat-confirmation/summary.json')
cost = load('cost-summary.json')
n = native['completed_samples'] if native else 0
nr_done = bool(protein and protein['nr_complete'])
required_results = ['native-validation/summary.json', 'validation-scaling/summary.json', 'validation-methods/summary.json', 'protein200/original-read-overlap.json', 'protein200-nt20/summary.json', 'protein-nt-panel/comparison.json', 'protein-resolution/summary.json', 'protein-cases/summary.json', 'protein-repeat-gate/summary.json', 'protein-native-losses/summary.json', 'preparation-discrepancies/nt-validation/comparison.json', 'blast-threading/summary.json', 'seed-nr-standalone/context-comparison.json', 'validation-methods-standalone/summary.json', 'protein-native-loss-alignments/summary.json', 'protein-screen-boundaries/summary.json']
done = n == 1000 and nr_done and bool(full and full['viral_only_NR_exact_pairs_absent_from_native'] is not None) and bool(repeat and repeat['completed_paired_samples'] == 500) and bool(cost and cost.get('all_owned_instances_terminated') and cost.get('private_backups_verified')) and all(load(name) is not None for name in required_results) and bool(load('protein200-nt20/summary.json') and load('protein200-nt20/summary.json')['nr_complete']) and bool(load('protein-nt-panel/comparison.json') and load('protein-nt-panel/comparison.json')['nr_complete']) and bool(load('protein200/original-read-overlap.json') and load('protein200/original-read-overlap.json')['repeat_gate_included'])
done = done and all(bool(load(name) and load(name).get('initial_screen_included')) for name in ['protein200/original-read-overlap.json','protein200-nt20/original-read-overlap.json'])
loss_readiness = load('protein-native-losses/summary.json') or {}
validation_readiness = load('native-validation/complete.json') or {}
screen_readiness = load('screen-burden/summary.json') or {}
done = done and loss_readiness.get('nr_controls_complete', False) and loss_readiness.get('raw_screen_membership_measured', False) and validation_readiness.get('completed_samples') == 1000 and screen_readiness.get('samples') == 1000
done = done and all(bool(load(name) and load(name).get('repeat_gate_included')) for name in ['protein200/original-read-overlap.json','protein200-nt20/original-read-overlap.json'])
done = done and all(bool(load(name) and 'pending' not in load(name).get('groups', {})) for name in ['protein200/original-read-overlap.json','protein200-nt20/original-read-overlap.json'])
state = 'Completed experiment' if done else ('Experiment ended with incomplete measurements' if cost and cost.get('all_owned_instances_terminated') else 'Interim report — incomplete results are explicitly marked')

findings = []
if native:
    a = native['arms']; u = a['union']['totals']; nt = a['nt']['totals']; aa = a['aa12']['totals']; short = a['nt20']['totals']
    findings.append(f'<li><strong>Changing only the first screen:</strong> NT24 OR AA12 adds {u["added_pairs"]:,} final read pairs to {nt["final_pairs"]:,} NT24 assignments ({100*u["added_pairs"]/nt["final_pairs"]:.2f}%) across {n:,} samples, losing {u["lost_pairs"]:,}. NT20 adds {short["added_pairs"]:,}. Only {u["added_beyond_nt20"]:,} AA-union additions are absent from NT20.</li>')
    findings.append(f'<li><strong>A replacement screen loses existing assignments:</strong> AA12 alone loses {aa["lost_pairs"]:,} NT24 pairs. The union retains every NT24 candidate at the screening boundary; final assignment preservation is measured separately.</li>')
if protein:
    v = protein['tiers']['viral_reference_candidate']
    if nr_done:
        vnr = protein['tiers']['nr_viral_only']
        findings.append(f'<li><strong>Protein processing:</strong> {v["added_exact_pairs"]:,} additional exact read pairs match the viral protein reference; {vnr["added_exact_pairs"]:,} have viral-only near-best NR evidence. These are environmental sequence assignments, not confirmed infections or novel viruses.</li>')
    else:
        findings.append(f'<li><strong>Protein processing:</strong> {v["added_exact_pairs"]:,} additional exact read pairs match the viral protein reference. Broad NR validation is pending; no validated gain is claimed.</li>')
if repeat and repeat['completed_paired_samples']:
    t = repeat['totals']
    findings.append(f'<li><strong>Repeat-filter confirmation:</strong> across {repeat["completed_paired_samples"]:,} prospectively selected samples, candidate volume falls {100*repeat["aggregate_candidate_reduction_fraction"]:.2f}%; {t["lost_original_union_pairs"]:,} original-union assignments and {t["lost_native_pairs"]:,} NT24 assignments are lost.</li>')
section('Measured findings', '<ul>' + ''.join(findings) + '</ul><p class="note">Read-pair counts, deduplicated support, and sample-species assignments answer different questions. None supplies known-truth sensitivity or specificity.</p>')

section('Question and comparison', '<p>Jeff suggested first changing the initial Nucleaze exact nucleotide k-mer screen to amino-acid space, with unknown strand and reading frame. This experiment separates that upstream change from replacing the later viral alignment and classification path with protein evidence.</p>' + table(['Comparison', 'Planned scope', 'Completed scope'], [
    ['Initial screen, unchanged nucleotide downstream', '1,000 BioSamples × 5 million pairs', f'{n:,} BioSamples; {n*5_000_000:,} pairs'],
    ['Protein viral processing and broad NR competition', 'First 200 fixed BioSamples', '200 with NR validation' if nr_done else ('200 viral-reference comparisons; NR pending' if protein else 'In progress')],
    ['Full-depth anchor', '247,633,776 pairs from SRR39694979', 'Native and viral-protein paths complete; see NR status below'],
    ['Secondary repeat-filter confirmation', 'Last 500 fixed BioSamples', repeat['completed_paired_samples'] if repeat else 0],
]) + '<p>All held-out samples come from public CASPER project <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1247874">PRJNA1247874</a>. The 1,000 BioSamples were selected by fixed SHA-256 rank from 2,149 eligible paired RNA-seq BioSamples, excluding the development and engineering pilots. Each uses four predetermined 1.25-million-pair blocks distributed across the run. This is five billion pairs at planned completion, not 1,000 full-depth runs. The samples span 31 location labels and collection dates from 2024-01-16 to 2026-06-30. The first 200 cover all 31 labels. The full-depth anchor is also a development BioSample and does not add an independent sample.</p>')

if native:
    labels = {'nt':'NT24 baseline','nt22':'NT22','nt20':'NT20','aa12':'AA12 alone','union':'NT24 OR AA12'}
    rows = []
    for key, label in labels.items():
        v = native['arms'][key]; t = v['totals']
        rows.append([label, t['final_pairs'], t['added_pairs'], t['lost_pairs'], t['joint_exemplars'], t['added_beyond_nt20'], v['samples_with_gain']])
    section('Initial-screen comparison', table(['Arm', 'Final pairs', 'Added to NT24', 'Lost from NT24', 'Joint duplicate exemplars', 'Added beyond NT20', 'Samples with gains'], rows) + '<p>NT20 and NT22 are simpler nucleotide controls. All arms use the same original FASTP, Bowtie2 viral alignment, contaminant filters, score threshold, LCA and duplicate processing after screening. AA12 requires an exact 12-residue match in any of six translations of a masked reference genome; a matching mate retains the whole pair. A longer exact peptide is not automatically a more selective biological criterion because translation collapses synonymous nucleotide differences and has strong composition effects.</p>' + image('nucleotide-comparison', 'Paired changes and computational burden. Screen CPU is measured; downstream CPU is estimated from task traces. Medians of stage CPU are not the median of total CPU.') + image('clustered-effects', 'Secondary uncertainty intervals resample complete location-label clusters. Repeated sites and dates are correlated; these intervals describe this cohort, not population prevalence.'))

species = load('heldout/species-summary.json')
if species:
    rows=[]
    for key,label in labels.items():
        t=species['counts'][f'{key}/minimum_3']
        rows.append([label,t['sample_species_calls'],t['gained_sample_species_calls'],t['lost_sample_species_calls'],t['human_species_calls'],t['gained_human_species_calls'],t['lost_human_species_calls']])
    section('Species-level support', table(['Arm', 'Sample-species calls', 'Gained', 'Lost', 'Human-annotated calls', 'Gained human-annotated', 'Lost human-annotated'], rows) + '<p>Shown at three or more original pipeline deduplicated clade exemplars. Support thresholds of one, two and three were analyzed descriptively. Human-infection annotation is a property of the fixed reference taxon; it does not diagnose a person or establish transmission.</p>')

support=load('heldout/support-transitions.json')
if support:
    rows=[]
    for key,label in labels.items():
        z=support['counts'].get(key+'/human_annotated',{})
        rows.append([label,z.get('gained_at_three',0),z.get('baseline_support_0',0),z.get('baseline_support_1',0),z.get('baseline_support_2',0)])
    section('Prior support for gained human-annotated species',table(['Arm','Gained at three exemplars','No baseline support','One baseline exemplar','Two baseline exemplars'],rows)+'<p>This descriptive breakdown separates additional species assignments from existing assignments that cross the three-exemplar threshold. Neither category establishes a biological detection.</p>')

if protein:
    rows=[]
    for tier,t in protein['tiers'].items():
        rows.append([tier.replace('_',' '),t['protein_exact_pairs'],t['added_exact_pairs'],t['lost_exact_pairs'],t['samples_with_gain']])
    body=table(['Evidence tier', 'Exact pairs', 'Added to native', 'Native exact pairs absent', 'Samples with gains'],rows)
    body+=table(['NR status', 'Exact pairs'],list(protein['nr_status'].items()))
else:
    body='<p>The first-200 protein comparison is still running. No completed broad-reference result is reported.</p>'
body+='<p>The protein path starts from AA12 candidates, retains original FASTP and human/other contaminant filtering, deduplicates complete trimmed pairs within each sample, aligns translated reads to six-frame viral reference segments, and assigns an LCA from eligible near-best hits. Each exact-sequence pair is searched once; original-read multiplicity describes membership rather than independent NR searches. Human/other contaminant filters precede protein alignment, whereas the native path performs viral nucleotide alignment first. Broad NR competition distinguishes viral-only, mixed, nonviral-only, unknown-taxonomy and unaligned evidence. Unaligned evidence is unvalidated. Near-best hits meet a 25-amino-acid minimum and 90% of the best eligible bitscore. Both-mate and taxonomically related subsets are secondary evidence summaries.</p><p>The native/protein comparison uses complete original FASTP pairs on both sides. Native TSVs can omit an unaligned mate, so comparing their stored sequence fields directly would inflate apparent gains. Exact-pair counts do not implement the native production coordinate/similarity duplicate policy. Keeping native assignments and adding NR-supported protein-only pairs is a distinct combined rule; it has not been validated as a diagnostic.</p>'
section('Full protein detection path',body + image('protein-evidence','NR competition and per-sample protein gains on common exact-pair units.'))

protein_nt20=load('protein200-nt20/summary.json')
if protein_nt20:
    rows=[]
    for tier,t in protein_nt20['tiers'].items():rows.append([tier.replace('_',' '),t['native_exact_pairs'],t['added_exact_pairs'],t['added_original_pairs'],t['lost_exact_pairs']])
    section('Protein comparison against NT20',table(['Evidence tier','NT20 exact pairs','Protein exact pairs absent from NT20','Original members absent from NT20','NT20 exact pairs absent'],rows)+'<p>The same first-200 protein outputs are also compared with the shorter nucleotide-screen control. Every count uses the same biological samples and complete original FASTP pairs. Additional candidate pairs remain unvalidated until broad-reference evidence is available.</p>')
boundary_rows=[]
for label,name in [('NT24','protein200/original-read-overlap.json'),('NT20','protein200-nt20/original-read-overlap.json')]:
    z=load(name)
    if not z or not z.get('initial_screen_included'):continue
    category='new_exact_pair_all_original_members_absent_from_native'
    for group in ['all','viral_only']:
        if group not in z['groups']:continue
        c=z['groups'][group]
        boundary_rows.append([label,group,c.get('initial_screen/'+category+'/none_retained',0),c.get('initial_screen/'+category+'/some_retained',0),c.get('initial_screen/'+category+'/all_retained',0),c.get('initial_screen_original_members/'+category+'/retained',0),c.get('initial_screen_original_members/'+category+'/excluded',0)])
if boundary_rows:
    section('Initial-screen versus later losses for protein additions',table(['Native comparator','Protein evidence','New exact queries with no member passing initial screen','Some members passed','All members passed','Original members passed screen but absent final calls','Original members absent screen'],boundary_rows)+'<p>Every query in this table has a new complete paired sequence and all of its original members are absent from the corresponding final native assignments. Raw FASTQ membership identifies whether those original IDs passed the first nucleotide screen; it does not identify which later operation removed them. Original members share an exact-sequence protein query and are not independent NR searches. <a href="results/protein-screen-boundaries/original-member-screen-membership.tsv">Original-ID membership evidence</a> is retained.</p>')

protein_species=load('protein200/species-summary.json')
if protein_species:
    rows=[]
    for tier,z in protein_species['counts'].items():
        if tier.endswith('/minimum_3'):rows.append([tier.split('/')[0].replace('_',' '),z['native_species'],z['protein_species'],z['gained_species'],z['lost_species'],z['gained_human_species'],z['lost_human_species']])
    section('Protein species-level support',table(['Evidence tier','Native sample-species calls','Protein calls','Gained','Lost','Gained human-annotated','Lost human-annotated'],rows)+'<p>Three or more exact complete read-pair keys support each sample-species assignment on both sides. Species labels and human-infection annotations come from the fixed viral reference. These counts differ from the original pipeline coordinate/similarity duplicate policy and do not establish infection.</p>')

resolution=load('protein-resolution/summary.json')
if resolution:
    rows=[]
    for group,z in resolution['groups'].items():
        rows.append([group,z['pairs'],z.get('species_relation/same_species',0),z.get('species_relation/different_species',0),z.get('species_relation/NR_above_species',0),z.get('species_relation/unresolved',0)])
    section('Taxonomic resolution of viral-only NR evidence',table(['Native overlap','Exact pairs','Same species as viral reference','Different species','NR LCA above species','Unresolved'],rows)+'<p>Viral-only NR evidence does not necessarily corroborate the species assigned by the smaller viral reference. These counts use the NR snapshot taxonomy to distinguish matching species from broader viral LCAs. The protein species table above describes viral-reference labels; broad NR support alone is not an independent species confirmation. <a href="results/protein-resolution/resolution.tsv">Pair-level taxonomic resolution</a> is retained.</p>')

cases=load('protein-cases/summary.json')
if cases:
    section('Inspectable protein disagreements',table(['Descriptive stratum','Candidate pairs'],list(cases['strata'].items()))+f'<p>{cases["cases"]:,} fixed-hash illustrative cases retain original read IDs, sequence lengths, homopolymer measurements and eligible near-best NR alignments from each competing clade. These examples are not an accuracy sample. <a href="results/protein-cases/cases.tsv">Selected read pairs</a> · <a href="results/protein-cases/case-alignments.tsv">Alignment evidence</a>.</p>')

losses=load('protein-native-losses/summary.json')
if losses:
    rows=[]
    for boundary,z in losses['loss_boundaries'].items():
        rows.append([boundary.replace('_',' '),z['original_pairs'],z['native_AA12_assignment_preserved'],z.get('present_in_original_AA12_screen'),z.get('prepared_sequences_identical_to_native'),z['nr_control_exact_pairs']])
    section('Native assignments absent from protein candidates',table(['Measured boundary','Original pairs','Retained by native AA12 path','Present in raw AA12 screen','Prepared sequences identical to native','Selected NR control sequences'],rows)+'<p>These first-200 losses are measured by original read identity. Preparation loss and absence of an eligible viral-protein alignment are separate observations. Raw-screen membership and selected control NR results are added when their jobs complete. A control result belongs to its explicitly selected query; no unsearched read receives an extrapolated verdict. <a href="results/protein-native-losses/lost-original-pairs.tsv">All lost original IDs and measured boundaries</a>.</p>')
loss_alignments=load('protein-native-loss-alignments/summary.json')
if loss_alignments:
    section('Reported alignments behind protein losses',table(['Observed boundary','Prepared query pairs','Original read pairs'],[[name.replace('_',' '),z['prepared_queries'],z['original_pairs']] for name,z in loss_alignments['boundaries'].items()])+'<p>The saved viral DIAMOND output was inspected for prepared native reads absent from protein candidates. Matches below the fixed 25-amino-acid summary cutoff had already passed the search reporting criteria. No reported alignment means no row under the frozen search settings, not demonstrated absence of homology. Primary classifications were not changed. <a href="results/protein-native-loss-alignments/query-loss-boundaries.tsv">Every affected query and original ID</a> · <a href="results/protein-native-loss-alignments/best-subthreshold-alignments.tsv">Best subthreshold alignments</a>.</p>')
header=load('preparation-discrepancies/header-control-summary.json')
if header:
    section('Read-header effect in eight preparation losses',table(['Header format','Repeat','Human mapped pairs','Retained pairs'],[[x['format'],x['repeat'],x['mapped_to_human_pairs'],x['retained_pairs']] for x in header['replicates']])+'<p>Eight native AA12 assignments absent from protein preparation were replayed with identical sequences, qualities, order and explicit Bowtie2 seed 0. Changing only terminal /1 and /2 to space-separated mate comments changed human-filter retention in both replicates; alignment multisets were identical within each header format. Whole source FASTQ files were byte-identical between workers, and the reconstructed header-only input matched the viral-alignment SAM roundtrip output byte-for-byte. Bowtie2 documents <a href="https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">read-name-dependent pseudo-random initialization</a>; the observed result is consistent with that mechanism. These eight outcome-selected examples do not estimate frequency or establish which assignment is biologically correct. Primary outputs were not changed. <a href="results/preparation-discrepancies/replay-summary.json">Replay evidence</a>.</p>')
header_nt=load('preparation-discrepancies/nt-validation/comparison.json')
if header_nt:
    rows=[]
    for group,z in header_nt['groups'].items():
        t=z['status'];rows.append([group,z['pairs'],t.get('viral_only',0),t.get('mixed_viral_nonviral',0),t.get('nonviral_only',0),t.get('unknown_taxonomy',0),t.get('no_alignment',0)])
    section('Comprehensive nucleotide evidence for the eight header cases',table(['Method / stratum','Pairs','Viral only','Mixed','Nonviral only','Unknown','No alignment'],rows)+'<p>The same eight complete pairs were searched against full core_nt in both nucleotide modes. Megablast reported 15 of 16 mates, with best alignments spanning 105–151 nucleotides at 96.644–100% identity. Sensitive blastn additionally found a 62-nucleotide, 85.484%-identity match on the remaining mate to a reference labeled Timema monikensis (E-value 1.10e-7); that pair retains viral evidence on its other mate. The human-filter replay had reported local 30–34-nucleotide matches. Reference labels and these outcome-selected examples do not establish biological truth or population frequency. <a href="results/preparation-discrepancies/nt-validation/comparison.tsv">Per-read competing evidence</a>.</p>')

oom=load('failed-attempts/nr24-oom-20261022T034211/attempt.json')
if oom:
    section('Incomplete NR attempt and memory recovery',table(['Reference block reached','Wall seconds','CPU seconds','Peak resident memory GiB','Completed alignment output bytes'],[[104,oom['search_wall_seconds'],oom['user_CPU_seconds']+oom['system_CPU_seconds'],oom['maximum_RSS_KiB']/1024**2,oom['completed_output_bytes']]])+'<p>The kernel killed DIAMOND for exceeding available memory on the96GiB worker. This attempt produced no completed validation result. Its runtime remains part of overall experiment spending. Recovery uses256GiB memory and a64GiB swap fallback; restarting the unchanged query pool and search parameters requires exact hardware/thread fixture alignment parity first. <a href="results/failed-attempts/nr24-oom-20261022T034211/kernel-oom.log">Kernel evidence</a> and the failed command/runtime log are retained.</p>')

controls=load('protein-nr/native-control-validation.json')
if controls:
    section('Native controls in the same NR search', table(['NR evidence status','Control exact pairs'],list(controls['status'].items())) + f'<p>{controls["control_pairs"]:,} native exact-pair controls were selected by fixed SHA rank, up to ten per first-200 sample. They received the same full-reference protein search. These environmental assignments are comparison controls, not known-truth positives; the selection is balanced by sample rather than viral species or read abundance.</p>')

native_validation=load('native-validation/summary.json')
if native_validation:
    rows=[]
    for arm,groups in native_validation['arms'].items():
        t=groups['all'];g=groups['added_to_baseline']
        rows.append([arm,t.get('pairs',0),t.get('status/aligned',0),t.get('status/no_alignment',0),t.get('status/not_sampled',0),g.get('status/aligned',0),g.get('evidence/aligned_natural_LCA_viral',0)])
    section('Original nucleotide validation',table(['Arm','Final pairs','Aligned','Selected, no alignment','Not sampled','Added pairs aligned','Added pairs with viral natural LCA'],rows)+'<p>The original production validation selects up to 20 duplicate exemplars per species per pooled arm and searches joined/merged mates with megablast. Only each read’s own alignment is summarized. Unsampled reads remain unvalidated. A natural LCA outside the viral clade can reflect mixed viral and nonviral evidence; it does not necessarily mean nonviral-only.</p>')

scaling=load('validation-scaling/summary.json')
if scaling:
    section('Validation sampling cost',table(['Cap per species','Queries','Run','Wall seconds','BLAST CPU seconds','Allocated worker cost USD'],[[x['cap'],x['queries'],x['run'],x['wall_seconds'],x['cpu_seconds'],x['allocated_worker_cost_usd']] for x in scaling['runs']])+'<p>The same first-200 NT24 input is processed by original species partitioning, hash sampling and mate merging. Search uses unchanged production megablast parameters and the complete existing core_nt. Caps select nested query sets. Timings exclude separately recorded preparation and reference staging. Other native tasks share the worker; the repeated cap-20 search measures temporal variation. Allocated worker cost charges the entire instance and its storage for that wall time, so it is not an incremental invoice charge or an isolated service estimate.</p>')
threading=load('blast-threading/summary.json')
if threading:
    section('BLAST query-threading comparison',table(['Run','Mode','Wall seconds','CPU seconds','Peak RSS GiB','Rows differing from original cap20','Queries with differences'],[[x['run'],x['mode'],x['wall_seconds'],x['cpu_seconds'],x['peak_rss_kib']/1024**2,x['baseline_parity']['rows_only_left']+x['baseline_parity']['rows_only_right'],x['baseline_parity']['differing_query_count']] for x in threading['runs']])+f'<p>Mode 0 is the pinned BLAST automatic choice; mode 1 explicitly splits queries. Every trial uses the same 2,820 cap-20 queries, full core_nt and eight threads. The pinned BLAST executable warns that this database is larger than its suggested range for query threading. Exact output comparison includes duplicate alignment multiplicity. Mode-1 repeats have identical alignment multisets: {threading["mode1_repeat_parity"]["alignment_multiset_equal"]}. The worker and warm cache are shared; all new trials have a 400 GiB container memory ceiling. <a href="results/blast-threading/summary.json">Commands, resource measurements and differing query IDs</a>.</p>')
three=load('validation-methods/summary.json')
if three:
    rows=[]
    for group,z in three['groups'].items():
        rows.append([group,z['pairs'],z['status_disagreements'],z['taxid_disagreements'],z['blastn_viral_only_where_megablast_unaligned'],z['protein_viral_only_where_megablast_unaligned']])
    section('Three methods on identical reads',table(['Stratum','Pairs','Status disagreements','Taxid disagreements','blastn viral-only / megablast unaligned','Protein viral-only / megablast unaligned'],rows)+'<p>Every comparison uses the exact same available individual-mate sequences. Nucleotide searches use core_nt and protein search uses NR, so reference coverage and representation also differ. Per-read assignments and all disagreements are retained in <a href="results/validation-methods/disagreements.tsv">the disagreement table</a>. The protein queries share one broad-reference job; no separate protein runtime is claimed for this stratum.</p>')

standalone=load('seed-nr-standalone/context-comparison.json')
standalone_methods=load('validation-methods-standalone/summary.json')
if standalone:
    section('Standalone protein search and query-context comparison',table(['Pairs','Query mates','Pooled alignment rows','Standalone alignment rows','Mates with different raw alignments','Pairs with different status','Pairs with different taxid'],[[standalone[k] for k in ['pairs','queries','pooled_alignment_rows','standalone_alignment_rows','raw_alignment_differing_query_mates','status_disagreements','taxid_disagreements']]])+'<p>The existing 746-pair, 1,436-mate diagnostic was searched alone against the same complete NR snapshot after the primary search and its NR-host analyses finished. Query bytes, executable and search parameters were preserved; only the surrounding query pool changed. Pooled primary results are retained separately. <a href="results/seed-nr-standalone/context-disagreements.tsv">Per-read context differences</a>.</p>')
if standalone_methods:
    section('Three methods with standalone protein validation',table(['Stratum','Pairs','Status disagreements','Taxid disagreements','blastn viral-only / megablast unaligned','Protein viral-only / megablast unaligned'],[[group,z['pairs'],z['status_disagreements'],z['taxid_disagreements'],z['blastn_viral_only_where_megablast_unaligned'],z['protein_viral_only_where_megablast_unaligned']] for group,z in standalone_methods['groups'].items()])+'<p>This secondary comparison uses the standalone protein output on the same available mates. Hardware, reference content and scoring differ between nucleotide and protein searches. <a href="results/validation-methods-standalone/disagreements.tsv">Every disagreement</a> is retained.</p>')

panel=load('protein-nt-panel/comparison.json')
if panel:
    rows=[]
    for group,z in panel['groups'].items():
        t=z['status'];rows.append([group,z['pairs'],t.get('viral_only',0),t.get('mixed_viral_nonviral',0),t.get('nonviral_only',0),t.get('unknown_taxonomy',0),t.get('no_alignment',0)])
    section('Protein-only candidates searched in nucleotide space',table(['Method / stratum','Pairs','Viral only','Mixed','Nonviral only','Unknown','No alignment'],rows)+'<p>One SHA-ranked protein candidate absent from both NT24 and NT20 was selected per each of the first 200 samples, plus one shared control from each of 32 SHA-ranked samples. Selection was frozen before any full-NR outcome was available. Every nucleotide query is copied unchanged from the full-NR query pool. This sample-balanced panel tests whether nucleotide search recovers the protein-only candidates; it is not a read-weighted accuracy estimate. Protein preparation collapses identical complete paired sequences within each sample before search. Original-read counts describe membership of those exact-sequence queries, not separate NR searches for every original member. The <a href="results/protein-nt-panel/original-members.tsv">original read-ID mapping</a> is retained. <a href="results/protein-nt-panel/comparison.tsv">Per-read assignments</a> and <a href="results/protein-nt-panel/eligible-nucleotide-alignments.tsv">eligible nucleotide alignments</a> are retained.</p>')

panel_rescues=load('protein-nt-panel/rescued-protein-panel-pairs.json')
panel_resources=load('protein-nt-panel/resource-summary.json')
if panel_rescues:
    z=panel_rescues
    section('Nucleotide evidence recovered for protein-only candidates',table(['Pairs unaligned by megablast, viral-only by blastn','Both mates reported','Median strongest alignment nt','Median identity %','Strongest hit E-value above 1e-5'],[[z['queries'],z['both_mates_aligned'],z['best_alignment_nt_median'],z['best_identity_pct_median'],z['best_Evalue_above_1e_minus5']]])+'<p>The weakest of these strongest-per-pair alignments is 54 nucleotides with E-value 0.000745. Nucleotide searches preserve the pinned BLAST default reporting E-value of 10; protein searches use 1e-5. The strength comparison does not reclassify the nucleotide outputs. Taxonomic resolution and agreement with the viral-reference label vary; NR results remain a separate measurement. <a href="results/protein-nt-panel/rescued-protein-panel-pairs.tsv">Every recovered query, original read ID and strongest alignment</a>.</p>')
if panel_resources:
    section('Matched protein-candidate panel nucleotide search work',table(['Method','Query mates','Wall seconds','CPU seconds','Whole-worker allocation USD'],[[x['task'],x['queries'],x['wall_seconds'],x['cpu_seconds'],x['allocated_whole_worker_cost_usd']] for x in panel_resources['runs']])+'<p>Both searches use the same 464 query mates, complete core_nt, pinned BLAST executable and eight threads. The primary worker and warm cache are shared; these allocations are not additive incremental invoice charges.</p>')

seed=load('seed-diagnostic200/comparison.json')
if seed:
    rows=[]
    for group,t in seed['groups'].items():
        z=t['status'];rows.append([group,t['pairs'],t['samples'],z.get('viral_only',0),z.get('mixed_viral_nonviral',0),z.get('nonviral_only',0),z.get('unknown_taxonomy',0),z.get('no_alignment',0)])
    section('Paired first-200 validation diagnostic',table(['Search / stratum','Pairs','Samples','Viral only','Mixed','Nonviral only','Unknown','No alignment'],rows)+'<p>Identical available mate sequences were searched separately with megablast and blastn. Each stratum selects up to two SHA-ranked pairs per sample; AA-added and NT20-added strata overlap. Changing the BLAST task changes seed and scoring defaults, so this is not an isolated word-size ablation. The separate-mate diagnostic is not a replay of joined-mate production validation or an accuracy estimate.</p>')

rescued=load('seed-diagnostic200/rescued-AA-pairs.json')
if rescued:
    section('Alignment strength of the recovered AA-added pairs',table(['Measure','Observed value'],[
        ['Pairs unaligned by megablast and viral-only with blastn',rescued['pairs']],
        ['Pairs with both mates aligned',rescued['aligned_mates'].get('2',0)],
        ['Best-alignment length range (nt)',f"{rescued['best_length_nt']['min']}–{rescued['best_length_nt']['max']}"],
        ['Best-alignment identity range (%)',f"{rescued['best_identity_pct']['min']:.2f}–{rescued['best_identity_pct']['max']:.2f}"],
        ['Largest best-hit E-value',f"{rescued['best_evalue']['max']:.2g}"],
        ['Exact original taxonomic assignment',rescued['relation'].get('exact',0)],
        ['Ancestor-related assignment',rescued['relation'].get('ancestor_related',0)],
        ['Different lineage',rescued['relation'].get('divergent',0)],
    ])+'<p>These long, strong nucleotide matches support viral sequence homology that megablast did not report. Their taxonomic resolution varies; several LCAs are broad viral groups. This is not evidence of a new infection or an unknown virus. <a href="results/seed-diagnostic200/rescued-AA-pairs.tsv">Every recovered read ID, taxon label and alignment-strength measurement</a> is retained.</p>')
if seed and any('viral_and_artificial_only' in x for x in seed['groups'].values()):
    section('What the competing nucleotide taxa represent',table(['Method / stratum','Viral plus artificial constructs only','Any cellular competitor'],[[k,v.get('viral_and_artificial_only',0),v.get('any_cellular_competitor',0)] for k,v in seed['groups'].items()])+'<p>Artificial constructs are a separate taxonomic clade and may contain viral sequence. The primary mixed/nonviral labels above are unchanged; these annotations prevent treating every nonviral competitor as cellular evidence. Counts can overlap other primary evidence categories.</p>')

if full:
    section('Full-depth anchor',table(['Measure','Pairs'],[
        ['Input read pairs',247633776],['Native final original pairs',full['native_pairs']],['Native exact complete trimmed pairs',full['native_exact_complete_trimmed_pairs']],['Viral-protein candidate exact pairs',full['protein_viral_reference_candidates']],['Protein candidate exact pairs absent from native',full['candidate_exact_trimmed_pairs_absent_from_native']],['Native exact pairs absent from protein candidates',full['native_exact_trimmed_pairs_absent_from_protein_candidates']],['NR viral-only exact pairs absent from native',full['viral_only_NR_exact_pairs_absent_from_native']],
    ])+'<p>For the first-screen-only change, NT24 produced 28,155 final original pairs. NT20, AA12 and their union each produced the same 28,178 final pairs, adding 23. This single deep sample is an engineering and depth check, not a population estimate.</p>')

section('Repetitive candidates and a prospective correction', '<p>In the triggering Miami sample, AA12 admitted 2,299,193 of five million raw pairs. In a fixed random reservoir of 1,024 exact-unique candidate pairs, 1,021 matched the 12-residue poly-proline/poly-glycine words; 1,020 contained a nucleotide G/C homopolymer run of at least 36 bases. These observations identify the matching mechanism; they do not establish the sequencing origin of those reads.</p><p>The secondary index omits AA12 words containing ten or more copies of one residue: 6,900 of 93,380,138 words. On five development samples and the triggering outlier, every original-union final pair ID and taxid was preserved. In the outlier, union candidates fell from 2,299,601 to 84,995. The rule and last-500 confirmation manifest were frozen before any of those 500 samples started. The original 1,000-sample and 200-sample protein comparisons were left unchanged. Preservation with nucleotide downstream processing does not establish preservation of protein-only assignments; the separate raw-gate check measures that attrition without claiming full downstream parity.</p>' + image('candidate-burden','Raw candidate burden varies substantially by sample; averages alone conceal high-volume outliers.') + image('repeat-filter-development','Exploratory development checks. The triggering outlier helped select the rule and is not independent confirmation.') + (table(['Prospective confirmation measure','Value'],[
    ['Completed paired samples',repeat['completed_paired_samples']],['Planned samples',500],['Aggregate candidate reduction',f'{100*repeat["aggregate_candidate_reduction_fraction"]:.2f}%' if repeat.get('aggregate_candidate_reduction_fraction') is not None else None],['Lost original-union final pairs',repeat.get('totals',{}).get('lost_original_union_pairs')],['Lost NT24 final pairs',repeat.get('totals',{}).get('lost_native_pairs')],['Shared-pair taxid changes',repeat.get('totals',{}).get('taxid_changes')],['Shared-pair sequence changes',repeat.get('totals',{}).get('shared_query_sequence_changes')],['Primary duplicate-representative changes',repeat.get('totals',{}).get('primary_duplicate_exemplar_changes')],['Similarity duplicate-representative changes',repeat.get('totals',{}).get('similarity_duplicate_exemplar_changes')],['Joint-exemplar membership changes',repeat.get('totals',{}).get('joint_exemplar_membership_changes')],
]) if repeat else '<p>Prospective confirmation remains pending.</p>'))

overlap=load('protein200/original-read-overlap.json')
if overlap and overlap.get('repeat_gate_included'):
    rows=[]
    for status,z in overlap['groups'].items():
        rows.append([status,z['candidate_exact_pairs'],z.get('repeat_gate/all_retained',0),z.get('repeat_gate/some_retained',0),z.get('repeat_gate/none_retained',0),z.get('repeat_gate/new_exact_pair_all_original_members_absent_from_native/none_retained',0)])
    section('Repeat filtering and protein-candidate support',table(['NR evidence','Exact pairs','All original members retained','Some retained','None retained','New-to-native pairs with none retained'],rows)+'<p>This secondary diagnostic applies the repeat-filtered raw screen to every original member of the frozen first-200 protein candidates. It measures loss of observed raw-screen support, without rerunning downstream processing or changing any original NR outcome. Preserving native nucleotide assignments does not by itself establish preservation of protein-only evidence.</p>')

section('Why upstream gains can remain small', '<p>Protein screening can pass reads that still fail the unchanged nucleotide aligner’s exact seeding requirement. A non-organism synthetic coding-sequence control with synonymous substitutions retained its amino-acid sequence but lacked an exact 20-base match: AA12 and translated alignment retained it, while the original Bowtie2 seed setting did not. Shortening the diagnostic seed recovered it. This illustrates a mechanism, not measured biological sensitivity.</p><p>Among the first five AA-union additions absent from raw NT20, two obtained an exact 20-base match after FASTP correction. The other three reflected the difference between the human-masked screen reference and the DUST-masked alignment reference. Those five examined reads do not establish the explanation for every later addition.</p><p>The original nucleotide validation uses megablast and can fail to align added reads that a more sensitive blastn diagnostic aligns. The controlled pilot compared 60 AA-added pairs and 60 fixed baseline controls using separate available mates. Megablast aligned 34/60 added pairs, versus 60/60 with blastn; the blastn results were 56 viral-only, one nonviral-only and three mixed. Separate-mate diagnostics are not exact replays of the production joined-mate validation.</p>')

section('Methods, limitations and reproducibility', '<ul><li>Five development BioSamples and a separate engineering pilot selected AA12 before the held-out comparison. AA8 passed 93.35% of the engineering pilot and was unsuitable. The synthetic IID-DNA null is a computational control, not environmental specificity.</li><li>References include 73,815 fixed viral genomes. Six-frame viral protein segments use genetic code 1 and a 25-residue minimum. The translated reference contains 44,892,414 segments. NR uses the 2026-10-17 NCBI BLAST snapshot; the actual alias reports 1,172,457,802 sequences and 441,882,733,350 residues.</li><li>DIAMOND viral alignment uses very-sensitive mode; NR uses sensitive mode. Both retain the frozen E-value, coverage, genetic-code and index-chunk parameters. Changing index chunks from four to one was rejected after 107/4,875 development-pair NR classifications changed. Further collapsing duplicate queries across the pooled NR input was also rejected: it changed alignment output for two of 54,468 original reads and saved only about 3% in a concurrent-workload benchmark.</li><li>Disabling cyclic garbage collection accelerated full-depth viral aggregation with byte-identical TSV, JSON and FASTA output. Bounded NR score storage passed an independent original-implementation parity oracle. Resource optimizations were not assumed equivalent from documentation alone.</li><li>There is no known-truth pathogen panel and no reference holdout. Extra assignments cannot be interpreted as sensitivity gains, infection, or novel-virus discovery. Repeated geography, collection dates, site composition and reference contamination limit generalization.</li><li>The worktree contains a standalone amino-acid screen and private experimental workflows. Production INDEX/RUN/DOWNSTREAM interfaces were not changed. Original binaries, input/reference hashes, frozen protocols, later addenda, pair membership tables, raw evidence, task traces and logs are retained privately.</li></ul><p>Primary sample-level bootstrap intervals use 2,000 resamples. Secondary whole-location intervals use 5,000 resamples and a sample-weighted ratio estimator. Both use fixed seeds. The last-500 repeat-filter experiment is a separate prospective secondary comparison, not a retuning of primary results.</p>')

runtime=load('search-runtimes.json')
if runtime:
    rows=[[x['method'],x['scope'],x['queries'],x.get('wall_seconds'),x.get('cpu_seconds') if x.get('cpu_seconds') is not None else 'Not measured',x.get('allocated_worker_cost_usd')] for x in runtime['runs']]
    section('Search commands and measured time',table(['Method','Query scope','Sequences','Wall seconds','CPU seconds','Allocated worker cost USD'],rows)+'<p>Commands below are the recorded invocations. Allocated cost charges the full worker and its configured EBS for the search wall time. Concurrent work shares that allocation, so these amounts must not be added together as invoice costs. Reference staging, interrupted attempts and idle time are included only in the experiment cost ledger.</p>'+''.join('<details><summary>'+html.escape(x['method']+' · '+x['scope'])+'</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">'+html.escape(x['command'])+'</pre></details>' for x in runtime['runs']))

stages=load('protein-nr/stage-timings.json')
if stages:
    section('NR search-time reference indexing',table(['Recorded stage group','Reported seconds'],[[k.replace('_',' '),v['logged_seconds']] for k,v in stages['groups'].items()])+f'<p>The completed search records {stages["logged_chunks"]:,} query/reference/shape/index-chunk combinations. These stage durations include rebuilding reference histograms and seed arrays while searching the existing native BLAST database. They are not CPU seconds or a separately measured persistent-index build; total search wall time also includes remaining setup, masking, cleanup and output work.</p>')

reference=load('reference-build/summary.json')
if reference:
    rows=[[x['reference'],x['build_wall_seconds'],x['build_CPU_seconds'],x['index_bytes']/2**30] for x in reference['measured_builds']]
    rows.extend([['Comprehensive core_nt', 'Prebuilt; not measured', 'Not measured',reference['storage']['core_nt_database_bytes']/2**30],['Comprehensive NR','Prebuilt; search-time indexing included in search','Not measured separately',reference['NR_snapshot_bytes']/2**30]])
    section('Reference construction and storage',table(['Reference','Build seconds','Build CPU seconds','Database GiB'],rows)+f'<p>Six-frame viral translation took {reference["viral_translation_wall_seconds"]:.2f} seconds before DIAMOND index construction. NR was used in its native BLAST format; per-search reference indexing is included in the search runtime. Full core_nt and NR snapshots are retained during computation. No reduced competitive reference was substituted.</p>')

if cost:
    section('Cost and cleanup',table(['Measure','Value'],list(cost.items())))
else:
    section('Cost and cleanup','<p>The experiment is still active under a $500 cap. The user removed the 24-hour deadline. The spending safeguard on the active workers is 2026-10-23 10:00 UTC. Applying the same guard when the stopped NRworker restarts gives a conservative account-wide bound of $463.52, including billing-lag and other-charge reserves. Cost Explorer reported $68.23 month to date at the latest check; that is delayed account billing, not live experiment spending. Completed workers will stop earlier. The failed96GiB NRworker is being replaced by a256GiB configuration with a slightly lower hourly compute rate; its frozen search parameters are retained. Existing budget shutdown safeguards remain in force. Cost Explorer is delayed; the planning reserve includes prior unbilled account usage. This report does not yet assert final cost or compute cleanup.</p>')

section('Source and result files','<p><a href="results/heldout/summary.json">Native comparison summary</a> · <a href="results/heldout/sample-metrics.tsv">Paired sample metrics</a> · <a href="results/screen-burden/sample-metrics.tsv">Candidate volumes</a> · <a href="frozen-protocol.json">Frozen protocol</a> · <a href="repeat-confirmation-protocol.json">Prospective repeat-filter protocol</a> · <a href="official-regional-prices.json">Official AWS rate evidence</a></p><p>Public data context: <a href="https://securebio.org/mgs-data/">SecureBio metagenomic data</a> · <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1247874">NCBI BioProject</a>. Tool documentation: <a href="https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">Bowtie2</a> · <a href="https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options">DIAMOND</a>.</p>')

doc='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Amino-acid detection experiment</title><style>body{font:16px/1.6 system-ui,sans-serif;color:#172635;background:#f4f6f8;margin:0}main{max-width:1100px;margin:auto;background:white;padding:48px 5vw}h1{font-size:36px;line-height:1.15;letter-spacing:-.8px}h2{font-size:23px;margin-top:45px;border-top:1px solid #dfe5e8;padding-top:22px}p,li{max-width:95ch}li{margin:9px 0}.status{color:#526779;font-size:14px}.note{background:#eef4f7;border-left:4px solid #3977aa;padding:12px 18px}.table{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}th,td{text-align:right;padding:9px 10px;border-bottom:1px solid #dde4e9}th:first-child,td:first-child{text-align:left}th{background:#eef3f6}figure{margin:28px 0}img{max-width:100%;height:auto}figcaption{font-size:13px;color:#526779}a{color:#176c92}@media print{body{background:white}main{padding:0}section{break-inside:avoid}h2{break-after:avoid}a{color:inherit}}</style><main><h1>Amino-acid detection experiment</h1><p class="status">Private experiment · '+state+' · Generated '+datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')+'</p>'+''.join(sections)+'</main></html>'
(r/'report.html').write_text(doc)
print(json.dumps({'report':str(r/'report.html'),'native_samples':n,'nr_complete':nr_done,'fully_complete':done}))
