/*
Benchmark output compression in the slowest Illumina Python steps.

Every module script writes its .gz outputs through Python's gzip.open, which
defaults to level 9. This runs five production scripts unchanged on real
inputs, while hook/sitecustomize.py switches the writer per variant (gzip at
another level, Zstandard, or none). Inputs stay as the production gzip files,
so only output compression changes.

Steps, named after the production processes they mirror:
  - FILTER_VIRAL_SAM, PROCESS_VIRAL_BOWTIE2_SAM (EXTRACT_VIRAL_READS_SHORT)
  - CHECK_TSV_DUPLICATES, JOIN_TSVS, ADD_CONDITIONAL_TSV_COLUMN
    (SPLIT_VIRAL_TSV_BY_SELECTED_TAXID, downstream)
A step runs only when its input params are set.

Each step runs once per variant per replicate, with production's container and
resource requests. Outputs stay in the work directory; only summaries are
published:
  - trace.tsv:      per-task timing, CPU, memory and I/O (Nextflow trace)
  - hosts.tsv:      per-task wall time of the core command and host size
  - checksums.tsv:  per-output format, sizes, line count and decompressed MD5
*/

// A production module script, run from this checkout exactly as committed
def moduleScript(module, name) {
    file("${projectDir}/../../modules/local/${module}/resources/usr/bin/${name}", checkIfExists: true)
}

// Shell lines that put the hook on PYTHONPATH and select the variant
def hookEnv(variant) {
    """
    export PYTHONPATH="\$PWD/hook"
    export BENCH_COMPRESSION="${variant}"
    """
}

// Shell snippet that records the core command's wall time and the host size
// (Batch shares CPUs rather than pinning them, so /proc shows the instance).
def hostRecord(step, variant, rep) {
    """
    end_s=\$(date +%s.%N)
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${step}" "${variant}" "${rep}" \\
        "\$(awk -v a="\$start_s" -v b="\$end_s" 'BEGIN {printf "%.1f", b - a}')" \\
        "\$(grep -c ^processor /proc/cpuinfo)" \\
        "\$(awk '/MemTotal/ {printf "%.0f", \$2 / 1048576}' /proc/meminfo)" \\
        "\$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//')" > host.tsv
    """
}

process FILTER_VIRAL_SAM {
    label "biopython"
    label "single"
    tag "${variant}|rep${rep}"
    input:
        tuple path(sam), path(reads), val(variant), val(rep)
        path(script)
        path(hook, stageAs: "hook/sitecustomize.py")
    output:
        tuple val("FILTER_VIRAL_SAM"), val(variant), val(rep), path("out.sam.gz"), emit: output
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        ${hookEnv(variant)}
        start_s=\$(date +%s.%N)
        python3 ${script} ${sam} ${reads} out.sam.gz ${params.bt2_score_threshold}
        ${hostRecord("FILTER_VIRAL_SAM", variant, rep)}
        """
}

process PROCESS_VIRAL_BOWTIE2_SAM {
    label "pysam_biopython"
    label "single"
    tag "${variant}|rep${rep}"
    input:
        tuple path(sam), val(variant), val(rep)
        path(script)
        path(genome_meta)
        path(virus_db)
        path(hook, stageAs: "hook/sitecustomize.py")
    output:
        tuple val("PROCESS_VIRAL_BOWTIE2_SAM"), val(variant), val(rep), path("out.tsv.gz"), emit: output
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        ${hookEnv(variant)}
        start_s=\$(date +%s.%N)
        zcat ${sam} | python3 ${script} -m ${genome_meta} -v ${virus_db} -o out.tsv.gz --paired
        ${hostRecord("PROCESS_VIRAL_BOWTIE2_SAM", variant, rep)}
        """
}

process CHECK_TSV_DUPLICATES {
    label "python"
    label "single"
    tag "${variant}|rep${rep}"
    input:
        tuple path(tsv), val(variant), val(rep)
        path(script)
        path(hook, stageAs: "hook/sitecustomize.py")
    output:
        tuple val("CHECK_TSV_DUPLICATES"), val(variant), val(rep), path("out.tsv.gz"), emit: output
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        ${hookEnv(variant)}
        start_s=\$(date +%s.%N)
        python3 ${script} -i ${tsv} -f "seq_id" -o out.tsv.gz
        ${hostRecord("CHECK_TSV_DUPLICATES", variant, rep)}
        """
}

process JOIN_TSVS {
    label "python"
    label "single"
    tag "${variant}|rep${rep}"
    input:
        tuple path(tsv1), path(tsv2), val(variant), val(rep)
        path(script)
        path(hook, stageAs: "hook/sitecustomize.py")
    output:
        tuple val("JOIN_TSVS"), val(variant), val(rep), path("out.tsv.gz"), emit: output
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        ${hookEnv(variant)}
        start_s=\$(date +%s.%N)
        python3 ${script} ${tsv1} ${tsv2} aligner_taxid_lca left out.tsv.gz
        ${hostRecord("JOIN_TSVS", variant, rep)}
        """
}

process ADD_CONDITIONAL_TSV_COLUMN {
    label "python"
    label "single"
    tag "${variant}|rep${rep}"
    input:
        tuple path(tsv), val(variant), val(rep)
        path(script)
        path(hook, stageAs: "hook/sitecustomize.py")
    output:
        tuple val("ADD_CONDITIONAL_TSV_COLUMN"), val(variant), val(rep), path("out.tsv.gz"), emit: output
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        ${hookEnv(variant)}
        start_s=\$(date +%s.%N)
        python3 ${script} --input ${tsv} --chk-col taxid_species --match-val NA \\
            --if-col aligner_taxid_lca --else-col taxid_species --new-hdr selected_taxid \\
            --output out.tsv.gz
        ${hostRecord("ADD_CONDITIONAL_TSV_COLUMN", variant, rep)}
        """
}

process CHECKSUM {
    label "python"
    label "checksum_resources"
    tag "${step}|${variant}"
    input:
        tuple val(step), val(variant), val(rep), path(output)
    output:
        path("checksum.tsv")
    script:
        """
        set -euo pipefail
        table_checksum.py ${step} ${variant} ${output} > checksum.tsv
        """
}

workflow {
    // CLI params arrive as strings; "1" in a Groovy range becomes the character code 49
    reps = channel.of(1..(params.replicates as int))
    // (variant, rep)
    runs = channel.fromList(params.variants.tokenize(",")).combine(reps)
    // The pysam_biopython container has Python 3.13, which has no zstd module
    runs_no_zstd = runs.filter { variant, _rep -> !variant.startsWith("zstd") }
    hook = file("${projectDir}/hook/sitecustomize.py", checkIfExists: true)

    fvs_in = params.fvs_sam ? channel.of(tuple(file(params.fvs_sam), file(params.fvs_reads))) : channel.empty()
    fvs = FILTER_VIRAL_SAM(fvs_in.combine(runs), moduleScript("filterViralSam", "filter_viral_sam.py"), hook)

    pvbs_in = params.pvbs_sam ? channel.of(file(params.pvbs_sam)) : channel.empty()
    pvbs = PROCESS_VIRAL_BOWTIE2_SAM(
        pvbs_in.combine(runs_no_zstd),
        moduleScript("processViralBowtie2Sam", "process_viral_bowtie2_sam.py"),
        file("${params.ref_dir}/results/virus-genome-metadata-gid.tsv.gz"),
        file("${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"),
        hook,
    )

    ctd_in = params.ctd_tsv ? channel.of(file(params.ctd_tsv)) : channel.empty()
    ctd = CHECK_TSV_DUPLICATES(ctd_in.combine(runs), moduleScript("checkTsvDuplicates", "check_tsv_duplicates.py"), hook)

    jt_in = params.jt_tsv1 ? channel.of(tuple(file(params.jt_tsv1), file(params.jt_tsv2))) : channel.empty()
    jt = JOIN_TSVS(jt_in.combine(runs), moduleScript("joinTsvs", "join_tsvs.py"), hook)

    act_in = params.act_tsv ? channel.of(file(params.act_tsv)) : channel.empty()
    act = ADD_CONDITIONAL_TSV_COLUMN(
        act_in.combine(runs),
        moduleScript("addConditionalTsvColumn", "add_conditional_tsv_column.py"),
        hook,
    )

    // Outputs should be deterministic, so checksum one replicate only
    outputs = fvs.output.mix(pvbs.output, ctd.output, jt.output, act.output)
    sums = CHECKSUM(outputs.filter { _s, _v, rep, _f -> rep == (params.checksum_replicate as int) })

    sums.collectFile(name: "checksums.tsv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: true)
    fvs.host.mix(pvbs.host, ctd.host, jt.host, act.host)
        .collectFile(name: "hosts.tsv", storeDir: params.outdir, seed: "step\tvariant\trep\tcore_wall_s\thost_vcpus\thost_mem_gib\tcpu_model\n", sort: true)
}
