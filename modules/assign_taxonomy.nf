/*
 * Assign taxonomy to contigs that are longer than threshold length
 */ 

process DIAMOND {
    publishDir "${params.outdir}/diamond/", mode: 'copy', pattern: "${sample_id}.out"
    container 'buchfink/diamond:version2.0.13'
    tag { sample_id }
    cpus 12

    input:
    tuple val(sample_id), path(assembly)
 
    output:
    tuple val(sample_id), path("${sample_id}.out")

    script:

    """
    filter_sequences_by_length.sh -i ${assembly} -s ${params.scaf_len} -o ${sample_id}_scaffolds_${params.scaf_len}.fasta 
    diamond blastx -d ${params.db} -q ${sample_id}_scaffolds_${params.scaf_len}.fasta  -o ${sample_id}.out -p ${task.cpus} -f 6 --top 10 -e 10 -b 6 --more-sensitive --matrix BLOSUM45 --comp-based-stats 1 

    """
}

/*
 * Obtain the lowest common ancestor from top blast hits
 */ 

process LCA {
    publishDir "${params.outdir}/krona/", mode: 'copy', pattern: "${sample_id}.html"
    container 'nanozoo/krona:2.7.1--e7615f7'
    tag { sample_id }
    cpus 1

    input:
    tuple val(sample_id), path(diamond_output), path(magnitudes)

    output:
    tuple val(sample_id), path("${sample_id}.csv"), emit: magnitudes
    tuple val(sample_id), path("${sample_id}.html"), emit: krona

    script:

    """
    ktImportBLAST ${diamond_output} ${diamond_output}:${magnitudes} -o ${sample_id}.html
    ktClassifyBLAST -t 50 -p -o ${sample_id}.tsv ${diamond_output}
    mergeMagnitudes.py -l ${sample_id}.tsv -m ${magnitudes} -o ${sample_id}.csv

    """
}

/*
 * Fill in parental nodes in the taxonomy tree for each sequence
 */ 

process TAXONOMY {
    publishDir "${params.outdir}/taxonomy/", mode: 'copy', pattern: "${sample_id}_taxonomy.csv"
    tag { sample_id }
    cpus 1

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
    cpus 1

    input:
    tuple val(sample_id), path(taxonomy)

    output:
    path "*.csv"

    script:

    """
    summarizeTaxonomy.py -i ${taxonomy} -s ${sample_id}

    """
}
