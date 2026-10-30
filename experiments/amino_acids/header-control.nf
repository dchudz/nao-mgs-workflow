include { BOWTIE2 as HUMAN } from '../../modules/local/bowtie2'
workflow {
    reads = channel.of(tuple('slash1',file(params.original,checkIfExists:true)),tuple('slash2',file(params.original,checkIfExists:true)),tuple('space1',file(params.normalized,checkIfExists:true)),tuple('space2',file(params.normalized,checkIfExists:true)))
    settings = [par_string:'--local --very-sensitive-local -X 850 --seed 0',remove_sq:true,debug:false,interleaved:true,db_download_timeout:1200,suffix:'human']
    HUMAN(reads,"${params.ref_dir}/results/bt2-human-index",settings)
}
