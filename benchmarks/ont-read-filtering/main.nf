/*
Benchmark faster variants of two ONT read-filtering steps in
EXTRACT_VIRAL_READS_ONT: FILTLONG and MASK_FASTQ_READS.

Each variant runs on every sample in the samplesheet, params.replicates times,
with the same resource requests as production. Outputs stay in the work
directory; only summaries are published:
  - trace.tsv:      per-task timing, CPU, memory and I/O (Nextflow trace)
  - hosts.tsv:      per-task wall time of the core command and host size
  - checksums.tsv:  per-variant read/base counts and output checksums
*/

// Production thresholds, as passed in subworkflows/local/extractViralReadsONT
def filtlongArgs() { "--min_length 50 --max_length 15000 --min_mean_q 90" }
def maskArgs() { "window=25 entropy=0.55" }

// Shell snippet that records the core command's wall time and the host size
// (Batch shares CPUs rather than pinning them, so /proc shows the instance).
def hostRecord(variant, label, rep) {
    """
    end_s=\$(date +%s.%N)
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${variant}" "${label}" "${rep}" \\
        "\$(awk -v a="\$start_s" -v b="\$end_s" 'BEGIN {printf "%.1f", b - a}')" \\
        "\$(grep -c ^processor /proc/cpuinfo)" \\
        "\$(awk '/MemTotal/ {printf "%.0f", \$2 / 1048576}' /proc/meminfo)" \\
        "\$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//')" > host.tsv
    """
}

// Variants that run in the production Filtlong container (gzip only, no pigz)
process FILTLONG_PROD {
    label "filtlong_resources"
    container params.filtlong_container
    tag "${variant}|${label}|rep${rep}"
    input:
        tuple val(label), path(reads), val(variant), val(rep)
    output:
        tuple val("filtlong_${variant}"), val(label), val(rep), path("out.fastq.gz"), emit: reads
        path("host.tsv"), emit: host
    script:
        def verbose = variant == "current" ? "--verbose" : ""
        """
        set -euo pipefail
        start_s=\$(date +%s.%N)
        filtlong ${filtlongArgs()} ${verbose} ${reads} | gzip > out.fastq.gz
        ${hostRecord("filtlong_" + variant, label, rep)}
        """
}

// Variants that need pigz alongside the same Filtlong version
process FILTLONG_PIGZ {
    label "filtlong_resources"
    conda "bioconda::filtlong=0.2.1 conda-forge::pigz=2.8 conda-forge::procps-ng"
    tag "${variant}|${label}|rep${rep}"
    input:
        tuple val(label), path(reads), val(variant), val(rep)
    output:
        tuple val("filtlong_${variant}"), val(label), val(rep), path("out.fastq.gz"), emit: reads
        path("host.tsv"), emit: host
    script:
        def n = (params.filtlong_chunks ?: task.cpus) as int
        def cmd
        if (variant == "quiet_pigz") {
            cmd = "filtlong ${filtlongArgs()} ${reads} | pigz -p ${task.cpus} > out.fastq.gz"
        } else if (variant == "chunked") {
            // Deal reads round-robin into n chunks, filter the chunks in
            // parallel, then concatenate (concatenated gzip members are valid gzip).
            // Filtlong reads its input twice, so chunks must be files.
            cmd = """
            mkdir -p chunks filtered
            pigz -dc ${reads} | awk -v n=${n} '
                BEGIN { for (i = 0; i < n; i++) cmd[i] = "pigz -1 -p 1 > chunks/c" i ".fastq.gz" }
                { print | cmd[int((NR - 1) / 4) % n] }
                END { for (i = 0; i < n; i++) close(cmd[i]) }'
            ls chunks/*.fastq.gz | xargs -P ${n} -I{} sh -c \\
                'filtlong ${filtlongArgs()} "\$1" | pigz -p 1 > "filtered/\$(basename "\$1")"' _ {}
            cat filtered/*.fastq.gz > out.fastq.gz
            rm -rf chunks filtered
            """
        } else {
            error "Unknown FILTLONG_PIGZ variant: ${variant}"
        }
        """
        set -euo pipefail
        start_s=\$(date +%s.%N)
        ${cmd}
        ${hostRecord("filtlong_" + variant, label, rep)}
        """
}

// Multithreaded single-pass alternative. Filtlong's --min_mean_q 90 means a
// mean per-base error of at most 10%, which chopper expresses as -q 10.
process CHOPPER {
    label "filtlong_resources"
    conda "bioconda::chopper=0.14.1 conda-forge::pigz=2.8 conda-forge::procps-ng"
    tag "${variant}|${label}|rep${rep}"
    input:
        tuple val(label), path(reads), val(variant), val(rep)
    output:
        tuple val("filtlong_${variant}"), val(label), val(rep), path("out.fastq.gz"), emit: reads
        path("host.tsv"), emit: host
    script:
        """
        set -euo pipefail
        start_s=\$(date +%s.%N)
        pigz -dc ${reads} | chopper -q 10 -l 50 --maxlength 15000 -t ${task.cpus} | pigz -p ${task.cpus} > out.fastq.gz
        ${hostRecord("filtlong_" + variant, label, rep)}
        """
}

// All masking variants use the production BBTools container (BBMap 39.01, pigz)
process MASK {
    label "mask_resources"
    container params.bbtools_container
    tag "${variant}|${label}|rep${rep}"
    input:
        tuple val(label), path(reads), val(variant), val(rep)
    output:
        tuple val("mask_${variant}"), val(label), val(rep), path("out.fastq.gz"), emit: reads
        path("host.tsv"), emit: host
    script:
        def mem = task.memory.toGiga()
        def n = params.mask_chunks as int
        def cmd
        if (variant == "current") {
            // Same command as modules/local/maskRead
            cmd = "zcat -f ${reads} | bbmask.sh in=stdin.fastq out=stdout.fastq ${maskArgs()} -Xmx${mem}g | gzip > out.fastq.gz"
        } else if (variant == "pigz") {
            cmd = "pigz -dc ${reads} | bbmask.sh in=stdin.fastq out=stdout.fastq ${maskArgs()} -Xmx${mem}g | pigz -p ${task.cpus} > out.fastq.gz"
        } else if (variant == "native_gz") {
            // Let BBTools read and write gzip itself (it uses pigz when available)
            cmd = "bbmask.sh in=${reads} out=out.fastq.gz ${maskArgs()} threads=${task.cpus} -Xmx${mem}g"
        } else if (variant == "chunked") {
            // Deal reads round-robin into n bbmask processes, each holding 1/n
            // of the input in memory; no intermediate files.
            def xmx = Math.max(1, (mem * 0.85 / n) as int)
            def t = Math.max(1, (task.cpus / n) as int)
            cmd = """
            pigz -dc ${reads} | awk -v n=${n} '
                BEGIN { for (i = 0; i < n; i++) cmd[i] = "bbmask.sh in=stdin.fastq out=stdout.fastq ${maskArgs()} threads=${t} -Xmx${xmx}g 2> bbmask_" i ".log | pigz -p 2 > masked_" i ".fastq.gz" }
                { print | cmd[int((NR - 1) / 4) % n] }
                END { for (i = 0; i < n; i++) close(cmd[i]) }'
            if grep -l -E 'Exception|OutOfMemory' bbmask_*.log; then exit 1; fi
            cat masked_*.fastq.gz > out.fastq.gz
            rm -f masked_*.fastq.gz
            """
        } else if (variant == "bbduk") {
            // Streaming entropy masking with BBDuk instead of BBMask
            cmd = "pigz -dc ${reads} | bbduk.sh in=stdin.fastq int=f out=stdout.fastq entropy=0.55 entropywindow=25 entropyk=5 entropymask=t threads=${task.cpus} -Xmx8g | pigz -p ${task.cpus} > out.fastq.gz"
        } else {
            error "Unknown MASK variant: ${variant}"
        }
        """
        set -euo pipefail
        start_s=\$(date +%s.%N)
        ${cmd}
        ${hostRecord("mask_" + variant, label, rep)}
        """
}

// Counts and checksums of decompressed output; see bin/fastq_checksum.py
process CHECKSUM {
    label "checksum_resources"
    conda "conda-forge::python=3.12 conda-forge::pigz=2.8 conda-forge::procps-ng"
    tag "${variant}|${label}"
    input:
        tuple val(variant), val(label), val(rep), path(reads)
    output:
        path("checksum.tsv")
    script:
        """
        set -euo pipefail
        pigz -dc ${reads} | fastq_checksum.py ${variant} ${label} > checksum.tsv
        """
}

workflow {
    samples = channel.fromPath(params.samplesheet).splitCsv(header: true)
    // CLI params arrive as strings; "1" in a Groovy range becomes the character code 49
    reps = channel.of(1..(params.replicates as int))
    fl_variants = channel.fromList(params.filtlong_variants.tokenize(","))
    mask_variants = channel.fromList(params.mask_variants.tokenize(","))

    // (label, reads, variant, rep)
    fl_in = samples.map { r -> tuple(r.label, file(r.filtlong_input)) }.combine(fl_variants).combine(reps)
    mask_in = samples.map { r -> tuple(r.label, file(r.mask_input)) }.combine(mask_variants).combine(reps)

    fl_branch = fl_in.branch { _label, _reads, variant, _rep ->
        prod: variant in ["current", "quiet"]
        pigz: variant in ["quiet_pigz", "chunked"]
        chopper: variant == "chopper"
    }
    prod_out = FILTLONG_PROD(fl_branch.prod)
    pigz_out = FILTLONG_PIGZ(fl_branch.pigz)
    chopper_out = CHOPPER(fl_branch.chopper)
    mask_out = MASK(mask_in)

    // Outputs should be deterministic, so checksum one replicate only
    all_reads = prod_out.reads.mix(pigz_out.reads, chopper_out.reads, mask_out.reads)
    sums = CHECKSUM(all_reads.filter { _v, _l, rep, _r -> rep == (params.checksum_replicate as int) })

    sums.collectFile(name: "checksums.tsv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: true)
    prod_out.host.mix(pigz_out.host, chopper_out.host, mask_out.host)
        .collectFile(name: "hosts.tsv", storeDir: params.outdir, seed: "variant\tlabel\trep\tcore_wall_s\thost_vcpus\thost_mem_gib\tcpu_model\n", sort: true)
}
