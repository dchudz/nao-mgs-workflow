// Protein alignment on the same FASTP input as the nucleotide screening control.
include { BOWTIE2 as HUMAN } from '../../modules/local/bowtie2'
include { BOWTIE2 as OTHER } from '../../modules/local/bowtie2'

process PROTEIN_QUERY {
    label 'python'
    label 'single_cpu_16GB_memory'
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path('queries.fna'), path('queries.tsv.gz'), path('queries.fna.members.tsv.gz'), path(reads)
    script:
        """
        python ${params.pipeline_root}/experiments/amino_acids/protein_queries.py ${reads} queries.fna queries.tsv.gz --prefix ${sample}_
        """
}

process PROTEIN_ALIGN {
    label 'python'
    label 'max'
    input:
        tuple val(sample), path(fasta), path(metadata), path(members), path(reads)
    output:
        tuple val(sample), path('protein-hits.tsv'), path(metadata), path(members), path(reads)
    script:
        """
        if [[ -s ${fasta} ]]; then
        ${params.diamond} blastx --db ${params.protein_db} --query ${fasta} --out protein-hits.tsv \\
            --very-sensitive --algo 0 --threads ${task.cpus} --block-size 2 --index-chunks 4 \\
            --min-orf 1 --query-gencode 1 --evalue 1e-5 --query-cover 30 --top 10 \\
            --outfmt 6 qseqid sseqid bitscore evalue length pident qlen qstart qend sstart send qframe slen qseq_gapped sseq_gapped
        else
            : > protein-hits.tsv
        fi
        """
}

process PROTEIN_CANDIDATES {
    label 'python'
    label 'single_cpu_16GB_memory'
    input:
        tuple val(sample), path(hits), path(metadata), path(members), path(reads)
    output:
        tuple val(sample), path('candidates.fastq'), emit: reads
        tuple val(sample), path(hits), path(metadata), path(members), emit: evidence
    script:
        """
        python ${params.pipeline_root}/experiments/amino_acids/diamond_adapter.py --extract-only \\
            --hits ${hits} --query-metadata ${metadata} --members ${members} --reads ${reads} \\
            --sample ${sample} --reference ${params.ref_dir} --genome-lengths ${params.genome_lengths} --output candidates.fastq
        """
}

process PROTEIN_ADAPTER {
    label 'python'
    label 'single_cpu_32GB_memory'
    input:
        tuple val(sample), path(hits), path(metadata), path(members), path(reads)
    output:
        tuple val(sample), path("${sample}_virus_hits.tsv.gz")
    script:
        """
        python ${params.pipeline_root}/experiments/amino_acids/diamond_adapter.py \\
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
