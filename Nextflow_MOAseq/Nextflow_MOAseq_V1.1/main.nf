#!/usr/bin/env nextflow

//nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { SEQPURGE } from './modules/seqpurge'
include { NGMERGE } from './modules/ngmerge'
include { STAR_INDEX } from './modules/star_index'
include { STAR_ALIGN } from './modules/star_align'

// rep not merged fork
include { BAM_PROCESSING } from './modules/bam_processing'
include { EFFECTIVE_GENOME_SIZE } from './modules/effective_genome_size'
include { MACS3_PEAK_CALLING } from './modules/macs3_peak_calling'
include { BAM_COVERAGE } from './modules/bamCoverage'

// rep merged fork
include { BAM_PROCESSING_M } from './modules/rep_merged/bam_processing'
include { EFFECTIVE_GENOME_SIZE_M } from './modules/rep_merged/effective_genome_size'
include { MACS3_PEAK_CALLING_M } from './modules/rep_merged/macs3_peak_calling'
include { BAM_COVERAGE_M } from './modules/rep_merged/bamCoverage'


workflow {
    log.info """
    ==============================================================================================
    
            .      __  _                           
         /|/| /  )/_| __   _ _ _    '  _ /'  _ 
        /   |(__/(  |    _) (-(/ /)//)(-(//)(- 
                      / /  /           "
                                    
    ==============================================================================================

    INPUT PARAMETERS:
        - samplesheet : ${params.samplesheet}
        - genome : ${params.genome}
        - annotation file (if any) : ${params.annotation}
        - genome index (if any) : ${params.genome_index}
        - sample name : ${params.sample_name}
        - sample merged : ${params.rep_merged}
        - output directory : ${params.outdir}

    ==============================================================================================
    """.stripIndent()

    // Set input data
    

    // Create a channel from the CSV file
    def read_pairs_ch = Channel
                    .fromPath(params.samplesheet)
                    .splitCsv(header:true)
                    .map{ row -> tuple(row.sampleId, [ file(row.forward_read), file(row.reverse_read) ]) }


    // FastQC on raw reads
    FASTQC_RAW(read_pairs_ch)

    // Trimming with SeqPurge
    SEQPURGE(read_pairs_ch)

    // FastQC on trimmed reads
    FASTQC_TRIMMED(SEQPURGE.out.trimmed_reads)

    // NGmerge
    NGMERGE(SEQPURGE.out.trimmed_reads)
    
    genome_ch = Channel.fromPath(params.genome).collect()

    annotation_ch = params.annotation 
        ? Channel.fromPath(params.annotation) 
        : Channel.empty()
    
    if (params.genome_index) {
        genome_index_ch = Channel.fromPath(params.genome_index).collect()
    } else {
        STAR_INDEX(genome_ch, annotation_ch)
            genome_index_ch = STAR_INDEX.out.index.collect()
    }

    STAR_ALIGN(NGMERGE.out.merged_reads, genome_index_ch)
    
    
    if (params.rep_merged) {

        BAM_PROCESSING_M(STAR_ALIGN.out.bam.collect(), params.sample_name)
    
        EFFECTIVE_GENOME_SIZE_M(genome_ch, BAM_PROCESSING_M.out.afl)

        BAM_COVERAGE_M(BAM_PROCESSING_M.out.filtered_bam, BAM_PROCESSING_M.out.filtered_bam_bai, EFFECTIVE_GENOME_SIZE_M.out.egs)
    
        MACS3_PEAK_CALLING_M(BAM_PROCESSING_M.out.filtered_bam, EFFECTIVE_GENOME_SIZE_M.out.egs, BAM_PROCESSING_M.out.afl)

    } else {

        BAM_PROCESSING(STAR_ALIGN.out.bam)

        EFFECTIVE_GENOME_SIZE(BAM_PROCESSING.out.filtered_bam, genome_ch, BAM_PROCESSING.out.afl)

        BAM_COVERAGE(BAM_PROCESSING.out.filtered_bam, BAM_PROCESSING.out.filtered_bam_bai, EFFECTIVE_GENOME_SIZE.out.egs)
    
        MACS3_PEAK_CALLING(BAM_PROCESSING.out.filtered_bam, EFFECTIVE_GENOME_SIZE.out.egs, BAM_PROCESSING.out.afl)


    }


    
}