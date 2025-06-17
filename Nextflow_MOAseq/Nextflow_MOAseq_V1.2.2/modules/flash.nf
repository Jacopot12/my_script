process FLASH {
    tag "FLASH on ${sample_id}"
    label 'standard'
    container "quay.io/biocontainers/flash2:2.2.00--h577a1d6_7"
    publishDir "${params.outdir}/merged_reads", mode: 'copy', pattern: "*.extendedFrags.fastq.gz"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}.extendedFrags.fastq.gz", emit: merged_reads
    
    script:
    """
    flash2 ${reads[0]} ${reads[1]} -t ${task.cpus} --cap-mismatch-quals -M 100 -o ${sample_id} -z
    """
}


