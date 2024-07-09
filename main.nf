#!/usr/bin/env nextflow

nextflow.enable.dsl=2
version = '1.0'

log.info"""\

A S S E M B L E F L O W - P I P E L I N E    
=========================================
reads                            : ${params.reads}
output folder                    : ${params.outdir}
Diamond database		 : ${params.diamond_db}
KronaTools			 : ${params.krona_db}
minimal scaffold length          : ${params.scaf_len}
adapter sequences for trimming   : ${params.adapt}
"""



include { DEDUPLICATE; MERGE; TRIM } from './modules/preprocessing_reads'
include { ASSEMBLE; CIRCULAR; ALIGN; SAMTOOLS } from './modules/assembly'
include { DIAMOND; LCA; TAXONOMY; SUMMARIZE } from './modules/assign_taxonomy'


workflow {
    reads = Channel.fromFilePairs( params.reads, checkIfExists: true )

    DEDUPLICATE( reads )
    MERGE( DEDUPLICATE.out )
    merge_adapt = MERGE.out.map { sample_id, unmerged_reads, merged_reads -> tuple(sample_id, unmerged_reads, merged_reads, file(params.adapt)) }
    TRIM( merge_adapt )
    ASSEMBLE( TRIM.out )
    assemble_diamond_db  = ASSEMBLE.out.map { sample_id, assembly -> tuple(sample_id, assembly, file(params.diamond_db)) }
    DIAMOND( assemble_diamond_db )
    CIRCULAR( ASSEMBLE.out )
    ASSEMBLE.out | join( TRIM.out, by: 0 ) | ALIGN | SAMTOOLS
    diamond_samtools = DIAMOND.out.join( SAMTOOLS.out.magnitudes, by: 0 )
    diamond_krona_db  = diamond_samtools.map { sample_id, diamond_output, magnitudes -> tuple(sample_id, diamond_output, magnitudes, file(params.krona_db)) }
    LCA ( diamond_krona_db )
    TAXONOMY( LCA.out.magnitudes )
    SUMMARIZE( TAXONOMY.out )
}

workflow.onComplete {
    log.info "[AssembleFlow] Pipeline Complete"
}
