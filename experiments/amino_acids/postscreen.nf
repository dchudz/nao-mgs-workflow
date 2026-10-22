include { EXTRACT_POSTSCREEN } from './postscreen-workflow'
include { MARK_VIRAL_DUPLICATES } from '../../subworkflows/local/markViralDuplicates'
include { VALIDATE_VIRAL_ASSIGNMENTS } from '../../subworkflows/local/validateViralAssignments'
workflow {
    reads = channel.of(tuple(params.sample, file(params.reads, checkIfExists: true)))
    results = EXTRACT_POSTSCREEN(reads, params.ref_dir, [adapters: file(params.adapters), minhits: 1, kmer_suffix: 'viral', db_download_timeout: 1200, aln_score_threshold: 20, taxid_artificial: '81077'])
    results.hits_final.view { "FINAL_HITS $it" }
    if (params.validate) {
        duplicates = MARK_VIRAL_DUPLICATES(results.hits_final, 1)
        validation = VALIDATE_VIRAL_ASSIGNMENTS(duplicates.hits, channel.value("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"), params.ref_dir, [platform:'illumina', validation_n_sample:20, blast_perc_id:60, blast_qcov_hsp_perc:30, blast_min_frac:0.9, blast_max_rank:10, taxid_artificial:81077, db_download_timeout:1200])
        validation.annotated_hits.view { "VALIDATED_HITS $it" }
    }
}
