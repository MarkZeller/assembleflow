#!/usr/bin/env nextflow

nextflow.enable.dsl=2
version = '1.1'

log.info"""\

A S S E M B L E F L O W - P I P E L I N E    
=========================================
reads                            : ${params.reads}
output folder                    : ${params.outdir}
Diamond database		 : ${params.diamond_db}
adapter sequences for trimming   : ${params.adapt}
"""



include { MERGE; TRIM } from './modules/preprocessing_reads'
include { ASSEMBLE; CIRCULAR; ALIGN; SAMTOOLS } from './modules/assembly'
include { DIAMOND; TAXONOMY } from './modules/assign_taxonomy'


workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    MERGE( reads )
    merge_adapt = MERGE.out.map { sample_id, unmerged_reads, merged_reads -> tuple(sample_id, unmerged_reads, merged_reads, file(params.adapt)) }
    TRIM( merge_adapt )
    ASSEMBLE( TRIM.out )
    CIRCULAR( ASSEMBLE.out )
    ASSEMBLE.out | join( TRIM.out, by: 0 ) | ALIGN | SAMTOOLS
    
    // Join assemblies and diamond_db for DIAMOND process
    ASSEMBLE.out.map { sample_id, assembly -> tuple(sample_id, assembly, file(params.diamond_db)) } \
        .set { diamond_inputs }

    DIAMOND ( diamond_inputs )
    taxonomy_input_ch = DIAMOND.out.join(SAMTOOLS.out.magnitudes, by: 0)
        .map { sample_id, diamond_out, magnitudes -> tuple(sample_id, diamond_out, magnitudes) }
    TAXONOMY ( taxonomy_input_ch )
}

workflow.onComplete {
    log.info "[AssembleFlow] Pipeline Complete"
}
