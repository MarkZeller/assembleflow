# Assembleflow

Assembleflow is a metagenomic sequencing assembly pipeline for Illumina paired-end reads. It assembles reads into contigs and assigns a taxonomic rank to each contig. The pipeline is implemented in Nextflow and is designed for high-throughput analysis.

## Requirements
- Conda
- Nextflow

## Installation

Before running Assembleflow, you need to build a DIAMOND database for taxonomic assignment. You can do this by running the following script:

```sh
bash build_diamond_db.sh
```
**Note**: The DIAMOND database requires hundreds of gigabytes of disk space. Make sure you have sufficient storage before proceeding.

## Pipeline Overview

Assembleflow consists of the following steps:

1. **Preprocessing Reads**
   - Merging overlapping read pairs
   - Adapter trimming

2. **Assembly**
   - Assemble reads into contigs
   - Align reads back to contigs
   - Detect circular contigs

3. **Taxonomic Assignment**
   - Assign taxonomy using DIAMOND
   - Retrieve taxonomic lineage information

## Configuration & Usage

You need to update the `nextflow.config` file before running Assembleflow.

### Key Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `reads`   | Folder containing paired-end sequence reads | `./test/*{R1,R2}.fastq.gz` |
| `adapt`   | File containing adapter sequences (FASTA format) | `./bin/adapter_sequences.fasta` |
| `db`      | Path to the DIAMOND database | `./refseq_protein_db/refseq_protein_diamond.dmnd` |
| `outdir`  | Output directory for results | `./test/results/` |

## Running the Pipeline

Once your `nextflow.config` file is properly set up, run the pipeline with:

```sh
nextflow run main.nf
```
