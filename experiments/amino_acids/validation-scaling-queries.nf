include { SPLIT_VIRAL_TSV_BY_SELECTED_TAXID } from '../../subworkflows/local/splitViralTsvBySelectedTaxid'
include { DOWNSAMPLE_VIRAL_ASSIGNMENTS } from '../../subworkflows/local/downsampleViralAssignments'
include { CONCATENATE_FILES_BY_EXTENSION } from '../../modules/local/concatenateFilesByExtension'
include { COPY_FILE as COPY_QUERIES } from '../../modules/local/copyFile'
workflow {
    hits = channel.fromPath(params.manifest).splitCsv(header:true,sep:'\t').map { row -> tuple(row.id,file(row.path,checkIfExists:true)) }
    split = SPLIT_VIRAL_TSV_BY_SELECTED_TAXID(hits,channel.value("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"))
    sampled = DOWNSAMPLE_VIRAL_ASSIGNMENTS(split.tsv,params.validation_n_sample as Integer,channel.of(false),'seq_id,sim_dup_exemplar')
    joined = CONCATENATE_FILES_BY_EXTENSION(sampled.fasta,'sampled_reads')
    COPY_QUERIES(joined.output,'queries.fasta.gz')
}
