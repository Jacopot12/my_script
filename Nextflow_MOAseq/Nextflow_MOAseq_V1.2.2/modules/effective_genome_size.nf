process EFFECTIVE_GENOME_SIZE {
    tag "Effective genome size"
    label 'standard'
    container "quay.io/biocontainers/khmer:3.0.0a3--py311heabec7a_7"
    publishDir "${params.outdir}/bam_parsing", mode: 'copy', pattern: "egs_output.txt"
    
    input:
    path bam
    path genome
    val afl
    
    output:
    env EGS, emit: egs
    
    script:
    name = bam.baseName

    """
    unique-kmers.py -q -k ${afl} -R ${name}.egs_output.txt ${genome}
    EGS=\$(grep "number of unique k-mers" ${name}.egs_output.txt | awk '{print \$5}')
    """
}