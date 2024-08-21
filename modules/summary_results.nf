/*
 * Merge LCA summary tables to produce single output table that includes all samples, then subset this to produce a separate table that is virus-specific
 */
 process SUMMARY {
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "*.csv"
    tag { sample_id }

    input:
    path lca_summary_files
    val taxonomy_complete

    output:
    path 'merged_lca_summary.csv'

    script:

    """
    mergeLCAs.py --output merged_lca_summary.csv ${lca_summary_files.join(' ')}

    """

 }

/*
 * Creates virus-specific summary table from the merged LCA summary table
 */
 process VIRUSES {
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: "virus_summary.csv"
    tag { sample_id }

    input:
    path merged_lca_summary

    output:
    path "virus_summary.csv"

    script:
    """
    filterViruses.py -i $merged_lca_summary -o virus_summary.csv

    """

 } 