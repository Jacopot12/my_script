process EFFECTIVE_GENOME_SIZE {
    tag "Effective genome size"
    container "quay.io/biocontainers/khmer:3.0.0a3--py311heabec7a_7"
    publishDir "${params.outdir}/bam_parsiong", mode: 'copy', pattern: "egs_output.txt"
    
    input:
    path genome
    val afl
    
    output:
    env EGS, emit: egs
    
    script:
    """
    unique-kmers.py -q -k ${afl} -R egs_output.txt ${genome}
    EGS=\$(grep "number of unique k-mers" egs_output.txt | awk '{print \$5}')
    """
}