process BAM_PROCESSING_M {
    tag "BAM processing"
    label 'standard'
    container "docker.io/staphb/samtools:1.21"
    publishDir "${params.outdir}/bam_parsiong", mode: 'copy', pattern: "*_max80_255.bam"
    publishDir "${params.outdir}/bam_parsiong", mode: 'copy', pattern: "*_max80_255.bam.bai"
    publishDir "${params.outdir}/bam_parsiong", mode: 'copy', pattern: "*_max80_255.txt"
    
    input:
    path bams
    val sample_name
    
    output:
    path "${sample_name}_max80_255.bam", emit: filtered_bam
    path "${sample_name}_max80_255.bam.bai", emit: filtered_bam_bai
    path "stats_${sample_name}_max80_255.txt", emit: stats
    env AFL, emit: afl
    
    script:
    """
    samtools merge -f -@ ${task.cpus} ${sample_name}.bam ${bams}
    
    samtools view -@ ${task.cpus} -h ${sample_name}.bam | \
    awk 'length(\$10) < 81 || \$1 ~ /^@/' | \
    samtools sort -@ ${task.cpus} - | \
    samtools view -@ ${task.cpus} -bS - -o ${sample_name}_max80.bam
    
    samtools index -@ ${task.cpus} ${sample_name}_max80.bam
    
    samtools view -@ ${task.cpus} -q 255 ${sample_name}_max80.bam -o ${sample_name}_max80_255.bam
    
    samtools index -@ ${task.cpus} ${sample_name}_max80_255.bam
    
    samtools stats ${sample_name}_max80_255.bam > stats_${sample_name}_max80_255.txt
    
    AFL=\$(awk -v OFS='\t' '{if(\$2=="average" && \$3=="first" && \$4=="fragment"){print \$6}}' stats_${sample_name}_max80_255.txt)
    """
}