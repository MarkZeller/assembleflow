/*
 * Assemble reads into contigs
 */ 

process ASSEMBLE {
    publishDir "${params.outdir}/assembly/", mode: 'copy', pattern:"${sample_id}_assembly.fasta"
    container 'staphb/spades:3.15.4'
    tag { sample_id }

    input:
    tuple val(sample_id), path(unmerged_reads), path(merged_reads)
 
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:

    """
    metaspades.py -t ${params.cpus} -1 ${unmerged_reads[0]} -2 ${unmerged_reads[1]} -s ${merged_reads} -o .
    mv scaffolds.fasta ${sample_id}_assembly.fasta

    """
}

/*
 * Identify circular contigs longer than 500nt based on kmer overlap
 */ 

process CIRCULAR {
    publishDir "${params.outdir}/circular/", mode: 'copy', pattern: "${sample_id}_circular.fasta"
    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "*_circular.fasta"
    
    script:

    """
    getCircular.py -i ${assembly} -k 55 -l 500 -o ${sample_id}_circular.fasta

    """
}

/*
 * Process to perform BWA indexing and alignment to generate SAM files
 */

process ALIGN {
    publishDir "${params.outdir}/alignments/", mode: 'copy'

    container 'staphb/bwa'
    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly), path(unmerged_reads), path(merged_reads)

    output:
    tuple val(sample_id), path("${sample_id}_unmerged.sam"), path("${sample_id}_merged.sam"), path(assembly), emit: sam_files

    script:
    """
    bwa index ${assembly}
    bwa mem ${assembly} ${unmerged_reads[0]} ${unmerged_reads[1]} -t ${params.cpus} > ${sample_id}_unmerged.sam
    bwa mem ${assembly} ${merged_reads} -t ${params.cpus} > ${sample_id}_merged.sam
    """
}

/*
 * Process to perform SAM to BAM conversion, sorting, merging, and indexing with Samtools
 */

process SAMTOOLS {
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"

    container 'staphb/samtools'
    tag { sample_id }

    input:
    tuple val(sample_id), path("${sample_id}_unmerged.sam"), path("${sample_id}_merged.sam"), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.magnitudes"), emit: magnitudes
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    samtools view -u -@ ${params.cpus} -S ${sample_id}_unmerged.sam | samtools sort -@ ${params.cpus} -o ${sample_id}_unmerged.bam -
    samtools view -u -@ ${params.cpus} -S ${sample_id}_merged.sam | samtools sort -@ ${params.cpus} -o ${sample_id}_merged.bam -
    samtools merge ${sample_id}.bam ${sample_id}_unmerged.bam ${sample_id}_merged.bam
    samtools index ${sample_id}.bam ${sample_id}.bam.bai
    samtools idxstats -@ ${params.cpus} ${sample_id}.bam | cut -f1,3 > ${sample_id}.magnitudes
    """
}






