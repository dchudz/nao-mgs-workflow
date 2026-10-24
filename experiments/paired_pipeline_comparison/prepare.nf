// Use the unchanged production cleaning module on the experimental screened input.
include { FASTP } from '../../modules/local/fastp'
workflow {
    reads = channel.of(tuple(params.sample, file(params.reads, checkIfExists:true)))
    cleaned = FASTP(reads, file(params.adapters, checkIfExists:true), true)
}
