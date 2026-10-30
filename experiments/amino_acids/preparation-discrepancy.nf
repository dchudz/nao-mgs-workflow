include { BOWTIE2 as VIRUS } from '../../modules/local/bowtie2'
include { BOWTIE2 as HUMAN_AFTER_VIRUS } from '../../modules/local/bowtie2'
include { BOWTIE2 as OTHER_AFTER_VIRUS } from '../../modules/local/bowtie2'
include { BOWTIE2 as HUMAN_DIRECT } from '../../modules/local/bowtie2'
include { BOWTIE2 as OTHER_DIRECT } from '../../modules/local/bowtie2'
workflow {
    reads = channel.of(tuple('diagnostic',file(params.reads,checkIfExists:true)))
    base = [remove_sq:true,debug:false,interleaved:true,db_download_timeout:1200]
    host = base + [par_string:'--local --very-sensitive-local -X 850']
    viral = VIRUS(reads,"${params.ref_dir}/results/bt2-virus-index",base+[par_string:'--local --very-sensitive-local --score-min G,0.1,19 -k 10 -X 850',suffix:'viral'])
    human_after = HUMAN_AFTER_VIRUS(viral.reads_mapped,"${params.ref_dir}/results/bt2-human-index",host+[suffix:'human_after_virus'])
    OTHER_AFTER_VIRUS(human_after.reads_unmapped,"${params.ref_dir}/results/bt2-other-index",host+[suffix:'other_after_virus'])
    human_direct = HUMAN_DIRECT(reads,"${params.ref_dir}/results/bt2-human-index",host+[suffix:'human_direct'])
    OTHER_DIRECT(human_direct.reads_unmapped,"${params.ref_dir}/results/bt2-other-index",host+[suffix:'other_direct'])
}
