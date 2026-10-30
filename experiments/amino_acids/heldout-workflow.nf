include { EXTRACT_POSTSCREEN } from './postscreen-workflow'
include { MARK_VIRAL_DUPLICATES } from '../../subworkflows/local/markViralDuplicates'
include { ADD_SAMPLE_COLUMN as ADD_GROUP_COLUMN } from '../../modules/local/addSampleColumn'
include { COUNT_READS_PER_CLADE } from '../../modules/local/countReadsPerClade'
workflow {
    reads = channel.fromPath(params.manifest).splitCsv(header:true,sep:'\t').map { row -> tuple(row.id,file(row.path,checkIfExists:true)) }
    results = EXTRACT_POSTSCREEN(reads,params.ref_dir,[adapters:file(params.adapters),minhits:1,kmer_suffix:'viral',db_download_timeout:1200,aln_score_threshold:20,taxid_artificial:'81077'])
    grouped = ADD_GROUP_COLUMN(results.hits_final,'group','grouped')
    marked = MARK_VIRAL_DUPLICATES(grouped.output,1)
    counts = COUNT_READS_PER_CLADE(marked.hits,channel.value("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"))
    marked.hits.view { "MARKED_HITS $it" }
}
