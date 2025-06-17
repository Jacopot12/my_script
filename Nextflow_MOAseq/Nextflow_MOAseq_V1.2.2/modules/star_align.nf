process STAR_ALIGN {
    tag "STAR alignment on ${sample_id}"
    label 'high_resources'
    container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    publishDir "${params.outdir}/star_align", mode: 'copy', pattern: "*_{Aligned.sortedByCoord.out.bam,Log.final.out,Log.out,Log.progress.out,SJ.out.tab}"
    
    input:
    path reads
    path index
    
    output:
    path "${sample_id}_Aligned.sortedByCoord.out.bam", emit: bam
    path "${sample_id}_Log.final.out", emit: log_final
    path "${sample_id}_Log.out", emit: log_out
    path "${sample_id}_Log.progress.out", emit: log_progress
    path "${sample_id}_SJ.out.tab", emit: SJ_table
    
    script:
    sample_id = reads.simpleName
    
    """
    STAR --readFilesCommand zcat \
         --genomeDir ${index} \
         --runThreadN ${task.cpus} \
         --readFilesIn ${reads} \
         --outSAMmultNmax 1 \
         --outFilterMultimapNmax 2 \
         --winAnchorMultimapNmax 100 \
         --outMultimapperOrder Random \
         --runRNGseed 100 \
         --outFileNamePrefix ${sample_id}_ \
         --outBAMsortingBinsN 15 \
         --alignIntronMax 1 \
         --outSAMtype BAM SortedByCoordinate
    """
}