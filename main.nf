#!/usr/bin/env nextflow
 
params.reads = './*_R{1,2}.fastq.gz'
params.adapt = ''
params.db = ''
params.scaf_len = ''
params.output = ''
params.cpus = 8
params.mem = '16GB'

log.info"""\

A S S E M B L E F L O W - P I P E L I N E    
===================================
database                : ${params.db}
reads                   : ${params.reads}
output folder           : ${params.output}
minimal scaffold length : ${params.scaf_len}
adapters for trimming   : ${params.adapt}
"""

Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Oops! Cannot find any file matching: ${params.reads}"  }
    .set { read_pairs }

adapters = file(params.adapt)
db = file(params.db)

process deduplicateReads {

    input:
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, '*deduplicated.fastq.gz' into deduplicated_reads

    script:

    """
    clumpify.sh in1=${reads[0]} in2=${reads[1]} out1=R1_deduplicated.fastq.gz out2=R2_deduplicated.fastq.gz dedupe subs=0 reorder
    """
}

process mergeReads {

    input:
    set pair_id, file(reads) from deduplicated_reads
 
    output:
    file 'merged.fastq.gz' into merged_reads
    set pair_id, '*unmerged.fastq.gz' into unmerged_reads


    script:

    """
    bbmerge.sh in1=${reads[0]} in2=${reads[1]} out=merged.fastq.gz outu1=R1_unmerged.fastq.gz outu2=R2_unmerged.fastq.gz vstrict

    """
}

process trimReads {

    input:
    file(merged) from merged_reads
    set pair_id, file(reads) from unmerged_reads
 
    output:
    set pair_id, '*_trimmed.fastq.gz' into (trimmed_reads_paired_for_alignment, trimmed_reads_paired_for_assembly)
    file 'trimmed.fastq.gz' into (trimmed_reads_single_for_alignment, trimmed_reads_single_for_assembly)

    script:

    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=R1_trimmed.fastq.gz out2=R2_trimmed.fastq.gz ref=${adapters} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo 
    bbduk.sh in=${merged} out=trimmed.fastq.gz ref=${adapters} ktrim=r k=23 mink=4 hdist=1 minlength=50 qtrim=w trimq=20 tpe tbo 

    """
}

process assembleReads {

    publishDir params.output, mode: 'copy', pattern: 'assembly.fasta'

    input:
    set pair_id, file(reads) from trimmed_reads_paired_for_assembly
    file(single_reads) from trimmed_reads_single_for_assembly
 
    output:
    file 'assembly.fasta' into (assembly_for_circular, assembly_for_alignment, assembly_for_diamond)

    script:

    """
    metaspades.py -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -s ${single_reads} -o spades_out
    mv spades_out/scaffolds.fasta ./assembly.fasta

    """
}

process circularReads {

    publishDir params.output, mode: 'copy', pattern: 'circular.fasta'

    input:
    file(assembly) from assembly_for_circular
 
    output:
    file 'circular.fasta' into circular_scaffolds
    
    script:

    """
    getCircular.py -i ${assembly} -k 55 -l 500 > circular.fasta

    """
}

process alignReads {
    publishDir params.output, mode: 'copy'

    input:
    set pair_id, file(reads) from trimmed_reads_paired_for_alignment
    file(single_reads) from trimmed_reads_single_for_alignment
    file(assembly) from assembly_for_alignment
 
    output:
    file 'alignment.bam' into aligned_reads_bam
    file 'alignment.bai' into aligned_reads_bai
    file 'alignment.magnitudes' into magnitudes_for_LCA

    script:

    """
    bwa index ${assembly}
    bwa mem ${assembly} ${reads[0]} ${reads[1]} -t ${task.cpus} | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o paired.bam -
    bwa mem ${assembly} ${single_reads} -t $task.cpus | samtools view -u -@ $task.cpus - | samtools sort -@ $task.cpus -o single.bam -
    samtools merge alignment.bam paired.bam single.bam 
    samtools index alignment.bam alignment.bai
    samtools idxstats -@ $task.cpus alignment.bam | cut -f1,3 > alignment.magnitudes

    """
}

process diamondBlast {
    publishDir params.output, mode: 'copy', pattern: 'diamond.out'

    input:
    file(assembly) from assembly_for_diamond
    val(db) from db
 
    output:
    file 'diamond.out' into diamond_for_LCA

    shell:

    '''
    awk '/^>/ {printf("\\n%s\\n",$0);next; } { printf("%s",$0);} END {printf("\\n");}' < !{assembly} | \
    awk -v var='!{params.scaf_len}' 'BEGIN {RS = ">" ; ORS = ""} length($2) >= var {print ">"$0}' > long_scaffolds.fasta
    diamond blastx -d !{db} -q long_scaffolds.fasta -o diamond.out -p !{task.cpus} -f 6 --top 10 -e 10 -b 6 --more-sensitive --matrix BLOSUM45 --comp-based-stats 1 

    '''
}

process getLCA {
    publishDir params.output, mode: 'copy', pattern: 'krona.html'

    input:
    file(diamond) from diamond_for_LCA
    file(magnitudes) from magnitudes_for_LCA
 
    output:

    file 'krona.html' into krona
    file 'diamond_magnitudes.tab' into lCA_for_taxonomy

    shell:

    '''
    ktImportBLAST !{diamond} !{diamond}:!{magnitudes} -o krona.html
    ktClassifyBLAST -t 50 -p -o diamond.tab !{diamond}
    awk 'NR==FNR { a[$1]=$2; next} $1 in a {print $0,"\\t"a[$1]}' !{magnitudes} diamond.tab > diamond_magnitudes.tab

    '''
}

process getTaxonomy {
    publishDir params.output, mode: 'copy', pattern: 'diamond_taxonomy.tab'

    input:
    file(lCA_magnitudes) from lCA_for_taxonomy

 
    output:
    file 'diamond_taxonomy.tab' into lCA_for_summarize

    
    script:

    """
    getTaxonomy.py -i ${lCA_magnitudes} > diamond_taxonomy.tab

    """
}

process summarizeTaxonomy {
    publishDir params.output, mode: 'copy'

    input:
    file(diamond_taxonomy) from lCA_for_summarize   
   
    output:
    file '*.tab'

    script:

    """
    summarizeTaxonomy.py -i ${diamond_taxonomy}

    """
}


workflow.onComplete {
    log.info "[AssembleFlow] Pipeline Complete"
}
