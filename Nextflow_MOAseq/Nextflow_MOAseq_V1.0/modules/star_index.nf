process STAR_INDEX {
    tag "STAR index"
    container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    publishDir "${params.outdir}/star_index", mode: 'copy', pattern: "genomeDir"
    
    input:
    path genome
    path annotation
    
    output:
    path "genomeDir", emit: index
    
    script:
    def annotation_arg = annotation.name != 'NO_FILE' ? "--sjdbGTFfile $annotation" : ''
    def memory_bytes = task.memory.toBytes()

    """
    mkdir genomeDir
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeFastaFiles ${genome} \
         --genomeDir genomeDir \
         --genomeSAindexNbases 13 \
         --limitGenomeGenerateRAM ${memory_bytes} \
         ${annotation_arg}
    """
}

