#!/usr/bin/env nextflow

nextflow.enable.dsl=2
version = '1.0'

log.info"""\

A S S E M B L E F L O W - P I P E L I N E    
=========================================
database                         : ${params.db}
reads                            : ${params.reads}
output folder                    : ${params.outdir}
minimal scaffold length          : ${params.scaf_len}
adapter sequences for trimming   : ${params.adapt}
"""



include { DEDUPLICATE; MERGE; TRIM } from './modules/preprocessing_reads'
include { ASSEMBLE; CIRCULAR; ALIGN } from './modules/assembly'
include { DIAMOND; LCA; TAXONOMY; SUMMARIZE } from './modules/assign_taxonomy'


workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    DEDUPLICATE( reads )
    MERGE( DEDUPLICATE.out )
    TRIM( MERGE.out )
    ASSEMBLE( TRIM.out )
    DIAMOND( ASSEMBLE.out )
    CIRCULAR( ASSEMBLE.out )
    ASSEMBLE.out | join( TRIM.out, by: 0 ) | ALIGN
    DIAMOND.out.join( ALIGN.out.magnitudes, by: 0 ) | LCA
    TAXONOMY( LCA.out.magnitudes )
    SUMMARIZE( TAXONOMY.out )
}

workflow.onComplete {
    log.info "[AssembleFlow] Pipeline Complete"
}
