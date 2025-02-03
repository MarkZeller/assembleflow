/*
 * Assign taxonomy to contigs that are longer than threshold length
 */ 

process DIAMOND {
    label 'high_cpu'
    publishDir "${params.outdir}/diamond/", mode: 'copy', pattern: "${sample_id}.out"
    conda "environment.yml"
    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly), path(diamond_db)
 
    output:
    tuple val(sample_id), path("${sample_id}.out")

    script:

    """
    diamond blastx \
      -q ${assembly} \
      -d ${diamond_db} \
      -o ${sample_id}.out \
      -p ${task.cpus} \
      --masking 0 \
      --unal 1 \
      --sensitive \
      -l 1 \
      -k 1 \
      -b 6 \
      -f 6 qseqid qlen qstrand sseqid slen evalue bitscore staxids sscinames full_qseq
    """
}

/*
 * Obtain taxonomy information for hits
 */ 

process TAXONOMY {
    label 'low_cpu'
    publishDir "${params.outdir}/taxonomy/", mode: 'copy', pattern: "*.csv"
    conda "environment.yml"
    tag { sample_id }

    input:
    tuple val(sample_id), path(diamond_out), path(magnitudes)

    output:
    tuple val(sample_id), path("${sample_id}.csv")
    
    script:

    """
    getTaxonomy.py -i ${diamond_out} -r ${magnitudes} -o ${sample_id}.csv

    """
}
