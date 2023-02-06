/*
 * Remove optical duplicates from paired-end sequence data
 */ 

process DEDUPLICATE {
    container 'staphb/bbtools:39.01'
    tag { sample_id }
    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*_deduplicated.fastq.gz')

    script:

    """
    clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=${sample_id}_R1_deduplicated.fastq.gz out2=${sample_id}_R2_deduplicated.fastq.gz dedupe subs=0 reorder threads=$task.cpus
    """
}

/*
 * Merge overlapping paired-end reads
 */ 

process MERGE {
    container 'staphb/bbtools:39.01'
    tag { sample_id }
    cpus 4

    input:
    tuple val(sample_id), path(reads)
 
    output:
    tuple val(sample_id), path('*_unmerged.fastq.gz'), path('*_merged.fastq.gz')

    script:

    """
    bbmerge.sh in1=${reads[0]} in2=${reads[1]} out=${sample_id}_merged.fastq.gz outu1=${sample_id}_R1_unmerged.fastq.gz outu2=${sample_id}_R2_unmerged.fastq.gz vstrict threads=$task.cpus

    """
}

/*
 * Trim reads based on quality, length, and the presence of adapter sequences
 */ 

process TRIM {
    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy', pattern: "${sample_id}_*merged_trimmed.fastq.gz"
    container 'staphb/bbtools:39.01'
    tag { sample_id }
    cpus 4

    input:
    tuple val(sample_id), path(unmerged_reads), path(merged_reads)

    output:
    tuple val(sample_id), path('*_unmerged_trimmed.fastq.gz'), path('*_merged_trimmed.fastq.gz')

    script:

    """
    bbduk.sh in1=${unmerged_reads[0]} in2=${unmerged_reads[1]} out1=${sample_id}_R1_unmerged_trimmed.fastq.gz out2=${sample_id}_R2_unmerged_trimmed.fastq.gz ref=${params.adapt} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo threads=$task.cpus
    bbduk.sh in=${merged_reads} out=${sample_id}_merged_trimmed.fastq.gz ref=${params.adapt} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo threads=$task.cpus

    """
}
