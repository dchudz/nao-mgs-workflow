include { VALIDATE_VIRAL_ASSIGNMENTS } from '../../subworkflows/local/validateViralAssignments'
workflow {
    // Input was duplicate-marked separately within each biological sample and arm.
    hits = channel.fromPath(params.manifest).splitCsv(header:true,sep:'\t').map { row -> tuple(row.id,file(row.path,checkIfExists:true)) }
    validation = VALIDATE_VIRAL_ASSIGNMENTS(hits,channel.value("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"),params.ref_dir,[platform:'illumina',validation_n_sample:20,blast_perc_id:60,blast_qcov_hsp_perc:30,blast_min_frac:0.9,blast_max_rank:10,taxid_artificial:81077,db_download_timeout:1200])
    validation.annotated_hits.view { "VALIDATED_HITS $it" }
}
