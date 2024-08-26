#!/bin/bash

# Set variables
REFSEQ_DB_NAME="refseq_protein"
REFSEQ_DB_DIR="refseq_protein_db"
FASTA_FILE="refseq_protein.fasta"
DIAMOND_DB_NAME="refseq_protein_diamond"
TAXONMAP_FILE="prot.accession2taxid.FULL.gz"
TAXONNODES_FILE="nodes.dmp"
TAXONNAMES_FILE="names.dmp"
THREADS="20"

# Create a directory for the database
mkdir -p "$REFSEQ_DB_DIR"

# Download the RefSeq Protein BLAST database files
echo "Downloading RefSeq Protein BLAST database files..."
for file in $(wget -qO- "ftp://ftp.ncbi.nlm.nih.gov/blast/db/" | grep -oP "${REFSEQ_DB_NAME}.*?\.tar\.gz"); do
    wget -P "$REFSEQ_DB_DIR" "ftp://ftp.ncbi.nlm.nih.gov/blast/db/$file"
done

# Check if wget was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download RefSeq Protein database files."
    exit 1
fi

# Extract the downloaded files
echo "Extracting database files..."
tar -xzvf "$REFSEQ_DB_DIR/${REFSEQ_DB_NAME}.*.tar.gz" -C "$REFSEQ_DB_DIR"

# Check if tar was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract RefSeq Protein database files."
    exit 1
fi

# Use blastdbcmd to create a FASTA file
echo "Creating FASTA file..."
blastdbcmd -db "$REFSEQ_DB_DIR/$REFSEQ_DB_NAME" -out "$FASTA_FILE" -outfmt "%f" -entry all

# Check if blastdbcmd was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to create FASTA file from BLAST database."
    exit 1
fi

# Download the taxonomy files
echo "Downloading taxonomy files..."
wget -P "$REFSEQ_DB_DIR" "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/$TAXONMAP_FILE"
wget -P "$REFSEQ_DB_DIR" "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"

# Check if wget was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download taxonomy files."
    exit 1
fi

# Unzip the taxonomy files
echo "Extracting taxonomy files..."
unzip "$REFSEQ_DB_DIR/taxdmp.zip" -d "$REFSEQ_DB_DIR"

# Check if unzip was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract taxonomy files."
    exit 1
fi

# Ensure necessary taxonomy files are available
if [ ! -f "$REFSEQ_DB_DIR/$TAXONMAP_FILE" ]; then
    echo "Error: Taxonmap file not found."
    exit 1
fi

if [ ! -f "$REFSEQ_DB_DIR/$TAXONNODES_FILE" ]; then
    echo "Error: Taxonnodes file not found."
    exit 1
fi

if [ ! -f "$REFSEQ_DB_DIR/$TAXONNAMES_FILE" ]; then
    echo "Error: Taxonnames file not found."
    exit 1
fi

# Build the DIAMOND database with taxonomy information
echo "Building DIAMOND database with taxonomy..."
diamond makedb --in "$FASTA_FILE" -d "$DIAMOND_DB_NAME" --taxonmap "$REFSEQ_DB_DIR/$TAXONMAP_FILE" --taxonnodes "$REFSEQ_DB_DIR/$TAXONNODES_FILE" --taxonnames "$REFSEQ_DB_DIR/$TAXONNAMES_FILE" --threads "$THREADS"

# Check if DIAMOND database build was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to build DIAMOND database."
    exit 1
fi

echo "DIAMOND database build completed successfully with taxonomy features."
