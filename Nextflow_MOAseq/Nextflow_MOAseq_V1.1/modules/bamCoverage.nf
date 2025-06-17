process BAM_COVERAGE {
    tag "bamCoverage"
    container "quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0"
    publishDir "${params.outdir}/bamCoverage", mode: 'copy'
    
    input:
    path bam
    path bai
    val egs
    
    output:
    path "${bam}.bigwig", emit: bigwig
    
    script:
    """
    bamCoverage \
    -b ${bam} \
    -o ${bam}.bigwig \
    -of bigwig  \
    -p ${task.cpus} \
    --effectiveGenomeSize ${egs}  \
    --minMappingQuality 0 \
    --normalizeUsing RPGC \
    --exactScaling \
    --smoothLength 0 \
    --binSize 1 \
    --minFragmentLength 20 \
    --maxFragmentLength 80
    """
}