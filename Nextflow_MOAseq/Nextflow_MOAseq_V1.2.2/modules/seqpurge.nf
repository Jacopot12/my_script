process SEQPURGE {
    tag "SeqPurge on ${sample_id}"
    label 'high_resources'
    container "quay.io/biocontainers/ngs-bits:2025_03--py313h572c47f_0"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.trim.fq.gz"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*trim.fq.gz"), emit: trimmed_reads
    path "*.trim.stats", emit: trim_stats

    script:
    """
    SeqPurge -min_len 20 -threads ${task.cpus} -qcut 20 \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 ${sample_id}_1.trim.fq.gz -out2 ${sample_id}_2.trim.fq.gz \
    -summary ${sample_id}.trim.stats -compression_level 5
    """
}