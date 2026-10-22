include { EXTRACT_POSTSCREEN } from './postscreen-workflow'
workflow {
    reads = channel.fromPath(params.manifest).splitCsv(header:true, sep:'\t').map { row -> tuple(row.id, file(row.path, checkIfExists:true)) }
    results = EXTRACT_POSTSCREEN(reads, params.ref_dir, [adapters:file(params.adapters),minhits:1,kmer_suffix:'viral',db_download_timeout:1200,aln_score_threshold:20,taxid_artificial:'81077'])
    results.hits_final.view { "FINAL_HITS $it" }
}
