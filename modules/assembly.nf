/*
 * Assemble reads into contigs
 */ 

process ASSEMBLE {
    publishDir "${params.outdir}/assembly/", mode: 'copy', pattern:"${sample_id}_assembly.fasta"
    container 'staphb/spades:3.15.4'
    tag { sample_id }
    cpus 12

    input:
    tuple val(sample_id), path(unmerged_reads), path(merged_reads)
 
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:

    """
    metaspades.py -t ${task.cpus} -1 ${unmerged_reads[0]} -2 ${unmerged_reads[1]} -s ${merged_reads} -o .
    mv scaffolds.fasta ${sample_id}_assembly.fasta

    """
}

/*
 * Identify circular contigs longer than 500nt based on kmer overlap
 */ 

process CIRCULAR {
    publishDir "${params.outdir}/circular/", mode: 'copy', pattern: "${sample_id}_circular.fasta"
    tag { sample_id }
    cpus 1

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
 * Align reads to contigs
 */ 

process ALIGN {
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"  
    container 'staphb/bwa:0.7.17'
    container 'staphb/samtools:1.9'
    tag { sample_id }
    cpus 8 

    input:
    tuple val(sample_id), path(assembly), path(unmerged_reads), path(merged_reads)
 
    output:
    tuple val(sample_id), path("${sample_id}.magnitudes"), emit: magnitudes
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:

    """
    bwa index ${assembly}
    bwa mem ${assembly} ${unmerged_reads[0]} ${unmerged_reads[1]} -t ${task.cpus} | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o ${sample_id}_unmerged.bam -
    bwa mem ${assembly} ${merged_reads} -t $task.cpus | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o ${sample_id}_merged.bam -
    samtools merge ${sample_id}.bam ${sample_id}_unmerged.bam ${sample_id}_merged.bam 
    samtools index ${sample_id}.bam ${sample_id}.bam.bai
    samtools idxstats -@ $task.cpus ${sample_id}.bam | cut -f1,3 > ${sample_id}.magnitudes

    """
}







