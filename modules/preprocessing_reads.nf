/*
 * Merge overlapping paired-end reads
 */ 

process MERGE {
    publishDir "${params.outdir}/merged_reads/", mode: 'copy', pattern: "${sample_id}_*merged.fastq.gz"
    container 'staphb/flash:latest'
    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)
 
    output:
    tuple val(sample_id), path('*_unmerged.fastq.gz'), path('*_merged.fastq.gz')

    script:

    """
    read_length=\$(zcat ${reads[0]} | awk 'NR==2 {print length(\$0)}')

    flash ${reads[0]} ${reads[1]} -o ${sample_id} --output-prefix=${sample_id} -z -m 10 -M \$read_length -x 0.1

    mv ${sample_id}.notCombined_1.fastq.gz ${sample_id}_R1_unmerged.fastq.gz 
    mv ${sample_id}.notCombined_2.fastq.gz ${sample_id}_R2_unmerged.fastq.gz 
    mv ${sample_id}.extendedFrags.fastq.gz ${sample_id}_merged.fastq.gz

    """
}

/*
 * Trim reads based on quality, length, and the presence of adapter sequences
 */ 

process TRIM {
    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy', pattern: "${sample_id}_*merged_trimmed.fastq.gz"
    container 'staphb/fastp:latest'
    tag { sample_id }

    input:
    tuple val(sample_id), path(unmerged_reads), path(merged_reads), path(adapt)

    output:
    tuple val(sample_id), path('*_unmerged_trimmed.fastq.gz'), path('*_merged_trimmed.fastq.gz')

    script:

    """
    # Unmerged reads trimming
    fastp --in1=${unmerged_reads[0]} --in2=${unmerged_reads[1]} --out1=${sample_id}_R1_unmerged_trimmed.fastq.gz --out2=${sample_id}_R2_unmerged_trimmed.fastq.gz --adapter_fasta=${adapt} --length_required 50 --cut_right --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --thread=${params.cpus} 

    # Merged reads trimming
    fastp --in1=${merged_reads} --out1=${sample_id}_merged_trimmed.fastq.gz --adapter_fasta=${adapt} --length_required 50 --cut_right --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --thread=${params.cpus}
    """
}
