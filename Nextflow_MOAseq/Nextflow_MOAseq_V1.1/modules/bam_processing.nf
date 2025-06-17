process BAM_PROCESSING {
    tag "BAM processing"
    container "docker.io/staphb/samtools:1.21"
    publishDir "${params.outdir}/bam_parsing", mode: 'copy', pattern: "*_max80_255.bam"
    publishDir "${params.outdir}/bam_parsing", mode: 'copy', pattern: "*_max80_255.bam.bai"
    publishDir "${params.outdir}/bam_parsing", mode: 'copy', pattern: "*_max80_255.txt"
    
    input:
    path bams
    
    output:
    path "${name}_max80_255.bam", emit: filtered_bam
    path "${name}_max80_255.bam.bai", emit: filtered_bam_bai
    path "stats_${name}_max80_255.txt", emit: stats
    env AFL, emit: afl
    
    script:
    name = bams.baseName

    """
    
    samtools view -@ ${task.cpus} -h ${name}.bam | \
    awk 'length(\$10) < 81 || \$1 ~ /^@/' | \
    samtools sort -@ ${task.cpus} - | \
    samtools view -@ ${task.cpus} -bS - -o ${name}_max80.bam
    
    samtools index -@ ${task.cpus} ${name}_max80.bam
    
    samtools view -@ ${task.cpus} -q 255 ${name}_max80.bam -o ${name}_max80_255.bam
    
    samtools index -@ ${task.cpus} ${name}_max80_255.bam
    
    samtools stats ${name}_max80_255.bam > stats_${name}_max80_255.txt
    
    AFL=\$(awk -v OFS='\t' '{if(\$2=="average" && \$3=="first" && \$4=="fragment"){print \$6}}' stats_${name}_max80_255.txt)
    """
}