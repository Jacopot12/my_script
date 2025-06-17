process MACS3_PEAK_CALLING {
    tag "MACS3 peak calling"
    container "quay.io/biocontainers/macs3:3.0.2--py310h397c9d8_2"
    publishDir "${params.outdir}/macs3", mode: 'copy'
    
    input:
    path bam
    val egs
    val afl
    
    output:
    path "*"
    
    script:
    max_gap = 2 * afl.toInteger()
    name = bam.baseName
    """
    macs3 callpeak \
    -t ${bam} \
    -n ${name} \
    -q 0.05 \
    -f BAM \
    -g ${egs} \
    -s ${afl} \
    --min-length ${afl} \
    --max-gap ${max_gap} \
    --nomodel \
    --extsize ${afl} \
    --keep-dup all \
    --buffer-size 10000000
    """
}