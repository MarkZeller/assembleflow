/*
 * Assemble reads into contigs
 */ 

process ASSEMBLE {
    label 'high_cpu'
    publishDir "${params.outdir}/assembly/", mode: 'copy', pattern:"${sample_id}_assembly.fasta"
    conda "environment.yml"
    tag { sample_id }

    input:
    tuple val(sample_id), path(unmerged_reads), path(merged_reads)
 
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:

    """
    metaspades.py -t ${task.cpus} -1 ${unmerged_reads[0]} -2 ${unmerged_reads[1]} -s ${merged_reads} -o . --only-assembler
    mv scaffolds.fasta ${sample_id}_assembly.fasta

    """
}

/*
 * Identify circular contigs longer than 500nt based on kmer overlap
 */ 

process CIRCULAR {
    label 'low_cpu'   
    publishDir "${params.outdir}/circular/", mode: 'copy', pattern: "${sample_id}_circular.fasta"
    conda "environment.yml"
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
    label 'low_cpu'
    conda "environment.yml"
    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly), path(unmerged_reads), path(merged_reads)

    output:
    tuple val(sample_id), path("${sample_id}_unmerged.sam"), path("${sample_id}_merged.sam"), path(assembly)

    script:
    """
    bwa index ${assembly}
    bwa mem ${assembly} ${unmerged_reads[0]} ${unmerged_reads[1]} -t ${task.cpus} > ${sample_id}_unmerged.sam
    bwa mem ${assembly} ${merged_reads} -t ${task.cpus} > ${sample_id}_merged.sam
    """
}

/*
 * Process to perform SAM to BAM conversion, sorting, merging, and indexing with Samtools
 */

process SAMTOOLS {
    label 'low_cpu'
    publishDir "${params.outdir}/alignments/", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"
    conda "environment.yml"
    tag { sample_id }

    input:
    tuple val(sample_id), path(unmerged), path(merged), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.magnitudes"), emit: magnitudes
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    samtools view -u -@ ${task.cpus} -S ${sample_id}_unmerged.sam | samtools sort -@ ${task.cpus} -o ${sample_id}_unmerged.bam -
    samtools view -u -@ ${task.cpus} -S ${sample_id}_merged.sam | samtools sort -@ ${task.cpus} -o ${sample_id}_merged.bam -
    samtools merge ${sample_id}.bam ${sample_id}_unmerged.bam ${sample_id}_merged.bam
    samtools index ${sample_id}.bam ${sample_id}.bam.bai
    samtools idxstats -@ ${task.cpus} ${sample_id}.bam | awk 'BEGIN{OFS=","} {print \$1,\$3}' > ${sample_id}.magnitudes
    """
}

