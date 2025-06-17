process FASTQC {
    tag "FastQC on ${sample_id}"
    label 'standard'
    container "quay.io/biocontainers/fastqc:0.11.9--0"
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: "*_fastqc.{zip,html}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}