process NGMERGE {
    tag "NGmerge on ${sample_id}"
    label 'standard'
    container "quay.io/biocontainers/ngmerge:0.3--0"
    publishDir "${params.outdir}/merged_reads", mode: 'copy', pattern: "*_merged.fq.gz"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}_merged.fq.gz", emit: merged_reads
    
    script:
    """
    NGmerge -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}_merged.fq.gz -p 0.2 -m 15 -d -e 30 -z -n ${task.cpus} -v
    """
}