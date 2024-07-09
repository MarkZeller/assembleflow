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
    diamond blastx -d ${diamond_db} -q ${sample_id}_scaffolds_${params.scaf_len}.fasta  -o ${sample_id}.out -p ${params.cpus} -f 6 --top 10 -e 10 -b 6 --more-sensitive --matrix BLOSUM45 --comp-based-stats 1 

    """
}

/*
 * Obtain the lowest common ancestor from top blast hits
 */ 

process LCA {
    publishDir "${params.outdir}/krona/", mode: 'copy', pattern: "${sample_id}.html"
    container 'nanozoo/krona:2.7.1--e7615f7'
    tag { sample_id }

    input:
    tuple val(sample_id), path(diamond_output), path(magnitudes),path(krona_db)

    output:
    tuple val(sample_id), path("${sample_id}.csv"), emit: magnitudes
    tuple val(sample_id), path("${sample_id}.html"), emit: krona

    script:

    """
    ktImportBLAST ${diamond_output} ${diamond_output}:${magnitudes} -o ${sample_id}.html -tax ${krona_db}
    ktClassifyBLAST -t 50 -p -o ${sample_id}.tsv ${diamond_output} -tax ${krona_db}
    mergeMagnitudes.py -l ${sample_id}.tsv -m ${magnitudes} -o ${sample_id}.csv

    """
}

/*
 * Fill in parental nodes in the taxonomy tree for each sequence
 */ 

process TAXONOMY {
    publishDir "${params.outdir}/taxonomy/", mode: 'copy', pattern: "${sample_id}_taxonomy.csv"
    tag { sample_id }

    input:
    tuple val(sample_id), path(lca_magnitudes)

    output:
    tuple val(sample_id), path("${sample_id}_taxonomy.csv")

    
    script:

    """
    getTaxonomy.py -i ${lca_magnitudes} -o ${sample_id}_taxonomy.csv

    """
}

/*
 * Summarize based on taxonomic level
 */ 

process SUMMARIZE {
    publishDir "${params.outdir}/taxonomy_summary/", mode: 'copy', pattern: "*.csv"
    tag { sample_id }

    input:
    tuple val(sample_id), path(taxonomy)

    output:
    path "*.csv"

    script:

    """
    summarizeTaxonomy.py -i ${taxonomy} -s ${sample_id}

    """
}
