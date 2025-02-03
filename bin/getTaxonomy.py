#!/usr/bin/env python3

import argparse
import csv
from ete3 import NCBITaxa

# Initialize NCBI taxonomy
ncbi = NCBITaxa()

# Taxonomic ranks to extract
ranks_to_extract = ["root", "cellular organisms", "superkingdom", "phylum", 
                    "class", "superfamily", "family", "genus", "species"]

def resolve_taxonomy(taxid):
    """Resolve taxonomy for a given TaxID."""
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        
        # Map ranks to their names
        rank_names = {ranks[taxid]: names[taxid] for taxid in lineage if ranks[taxid] in ranks_to_extract}
        
        # Fill the ranks_to_extract order
        resolved_ranks = [rank_names.get(rank, "NA") for rank in ranks_to_extract]
        return resolved_ranks
    except Exception as e:
        # Return NA for all ranks if an error occurs
        return ["NA"] * len(ranks_to_extract)

def load_read_counts(read_counts_file):
    """Load read counts from a CSV file into a dictionary."""
    read_counts = {}
    with open(read_counts_file, "r") as infile:
        reader = csv.reader(infile)
        for row in reader:
            if len(row) >= 2:
                sequence_name, count = row[0], row[1]
                read_counts[sequence_name] = int(count)
    return read_counts

def process_diamond_output(input_file, output_file, read_counts):
    """Process DIAMOND BLASTX output and generate a taxonomy CSV with read counts."""
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t", fieldnames=[
            "qseqid", "qlen", "qstrand", "sseqid", "slen", "evalue", 
            "bitscore", "staxids", "sscinames", "full_qseq"
        ])
        writer = csv.writer(outfile)
        
        # Write header to output CSV
        header = ["qseqid", "taxid", "evalue", "bitscore"] + ranks_to_extract + ["read_count"]
        writer.writerow(header)
        
        # Process each line of the DIAMOND output
        for row in reader:
            qseqid = row["qseqid"]
            staxids = row["staxids"]
            evalue = row["evalue"]
            bitscore = row["bitscore"]
            staxids = row["staxids"]
            
            # Handle cases with multiple TaxIDs (e.g., split by ';')
            taxid = int(staxids.split(";")[0]) if staxids else None
            if taxid:
                taxonomy = resolve_taxonomy(taxid)
            else:
                taxonomy = ["NA"] * len(ranks_to_extract)
            
            # Get read count
            read_count = read_counts.get(qseqid, "NA")
            
            # Write row to output CSV
            writer.writerow([qseqid, staxids, evalue, bitscore] + taxonomy + [read_count])

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process DIAMOND BLASTX output and generate a taxonomy CSV with read counts.")
    parser.add_argument("-i", "--input", required=True, help="Input DIAMOND BLASTX output file.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file.")
    parser.add_argument("-r", "--reads", required=True, help="Read counts CSV file.")
    args = parser.parse_args()

    # Load read counts
    read_counts = load_read_counts(args.reads)
    
    # Process the input file
    process_diamond_output(args.input, args.output, read_counts)
