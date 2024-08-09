/*
 * Assign taxonomy to contigs that are longer than threshold length
 */ 

process DIAMOND {
    maxForks 7

    publishDir "${params.outdir}/diamond/", mode: 'copy', pattern: "${sample_id}.out"
    container 'staphb/diamond'
    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly), path(diamond_db)
 
    output:
    tuple val(sample_id), path("${sample_id}.out")

    script:

    """
    filter_sequences_by_length.sh -i ${assembly} -s ${params.scaf_len} -o ${sample_id}_scaffolds_${params.scaf_len}.fasta 
    diamond blastx -d ${diamond_db} -q ${sample_id}_scaffolds_${params.scaf_len}.fasta  -o ${sample_id}.out -p ${params.cpus} --outfmt 6 qseqid sseqid staxids pident length mismatch gapopen evalue bitscore --top 10 -e 10 -b 6 --more-sensitive --matrix BLOSUM45 --comp-based-stats 1 
    """
}


/*
 * Obtain lowest common ancestor from blast hits and taxonomy information for hits and LCAs
 */ 

process TAXONOMY {
    publishDir "${params.outdir}/taxonomy/", mode: 'copy', pattern: "*.csv"
    tag { sample_id }

    input:
    tuple val(sample_id), path(diamond_out), path(magnitudes)

    output:
    tuple val(sample_id), path("${sample_id}_table.csv"), emit: table
    tuple val(sample_id), path("${sample_id}_lca_summary.csv"), emit: lca_summary
    tuple val(sample_id), path("${sample_id}_hits_summary.csv"), emit: hits_summary
    tuple val(sample_id), path("${sample_id}_virus_hits_summary.csv"), optional: true, emit: virus_hits_summary

    
    script:

    """
    taxonomyOut.py -i ${diamond_out} -m ${magnitudes} -o ${sample_id}_table.csv -l ${sample_id}_lca_summary.csv -s ${sample_id}_hits_summary.csv -v ${sample_id}_virus_hits_summary.csv

    """
}
