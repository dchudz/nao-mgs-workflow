include { BOWTIE2 as HUMAN } from '../../modules/local/bowtie2'
include { BOWTIE2 as OTHER } from '../../modules/local/bowtie2'
workflow {
    reads = channel.fromPath(params.manifest).splitCsv(header:true, sep:'\t').map { row -> tuple(row.id, file(row.path, checkIfExists:true)) }
    settings = [par_string:'--local --very-sensitive-local -X 850',remove_sq:true,debug:false,interleaved:true,db_download_timeout:1200]
    human = HUMAN(reads,"${params.ref_dir}/results/bt2-human-index",settings+[suffix:'human'])
    other = OTHER(human.reads_unmapped,"${params.ref_dir}/results/bt2-other-index",settings+[suffix:'other'])
    other.reads_unmapped.view { "PROTEIN_INPUT $it" }
}
