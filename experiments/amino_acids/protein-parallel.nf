// Execution-only parallel adapter; search processes are imported unchanged.
include { PROTEIN_QUERY; PROTEIN_ALIGN } from './protein-production'
include { BOWTIE2 as HUMAN } from '../../modules/local/bowtie2'
include { BOWTIE2 as OTHER } from '../../modules/local/bowtie2'

process PROTEIN_CANDIDATES {
    label 'python'
    label 'max'
    memory '240 GB'
    input:
        tuple val(sample), path(hits), path(metadata), path(members), path(reads)
    output:
        tuple val(sample), path('candidates.fastq.gz'), emit: reads
        tuple val(sample), path(hits), path(metadata), path(members), emit: evidence
    script:
        """
        python ${params.pipeline_root}/analysis-tools/parallel_diamond_adapter.py --workers ${task.cpus} --extract-only \\
            --hits ${hits} --query-metadata ${metadata} --members ${members} --reads ${reads} \\
            --sample ${sample} --reference ${params.ref_dir} --genome-lengths ${params.genome_lengths} --output candidates.fastq.gz
        """
}

process PROTEIN_ADAPTER {
    label 'python'
    label 'max'
    memory '240 GB'
    input:
        tuple val(sample), path(hits), path(metadata), path(members), path(reads)
    output:
        tuple val(sample), path("${sample}_virus_hits.tsv.gz")
    script:
        """
        python ${params.pipeline_root}/analysis-tools/parallel_diamond_adapter.py --workers ${task.cpus} \\
            --hits ${hits} --query-metadata ${metadata} --members ${members} --reads ${reads} \\
            --sample ${sample} --reference ${params.ref_dir} --genome-lengths ${params.genome_lengths} --output unsorted.tsv.gz
        zcat unsorted.tsv.gz | { IFS= read -r header; printf '%s\\n' "\${header}"; LC_ALL=C sort -k1,1; } | gzip > ${sample}_virus_hits.tsv.gz
        """
}

workflow {
    reads = channel.of(tuple(params.sample, file(params.reads,checkIfExists:true)))
    queries = PROTEIN_QUERY(reads)
    alignments = PROTEIN_ALIGN(queries)
    candidates = PROTEIN_CANDIDATES(alignments)
    settings = [par_string:'--local --very-sensitive-local -X 850',remove_sq:true,debug:false,interleaved:true,db_download_timeout:1200]
    human = HUMAN(candidates.reads,"${params.ref_dir}/results/bt2-human-index",settings+[suffix:'human'])
    other = OTHER(human.reads_unmapped,"${params.ref_dir}/results/bt2-other-index",settings+[suffix:'other'])
    assembled = candidates.evidence.join(other.reads_unmapped,by:0)
    hits = PROTEIN_ADAPTER(assembled)
    hits.view { "PROTEIN_FINAL $it" }
}
