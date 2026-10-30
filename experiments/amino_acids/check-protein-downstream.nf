include { MARK_VIRAL_DUPLICATES } from '../../subworkflows/local/markViralDuplicates'
include { VALIDATE_VIRAL_ASSIGNMENTS } from '../../subworkflows/local/validateViralAssignments'
workflow {
    hits = channel.of(tuple(params.sample,file(params.hits,checkIfExists:true)))
    marked = MARK_VIRAL_DUPLICATES(hits,1)
    settings=[platform:'illumina',validation_n_sample:20,blast_task:'blastn',blast_perc_id:60,blast_qcov_hsp_perc:30,blast_min_frac:0.9,blast_max_rank:10,taxid_artificial:81077,db_download_timeout:1200]
    checked = VALIDATE_VIRAL_ASSIGNMENTS(marked.hits,channel.value("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"),params.ref_dir,settings)
    checked.annotated_hits.view { "PROTEIN_DOWNSTREAM_CHECK $it" }
}
