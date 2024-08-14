# Assembleflow

Assembleflow is a metagenomics sequencing assembly pipeline for Illumina paired-end reads. It will assemble reads into contigs and assigns a taxonomic rank to each contig. 

## Table of Contents
* [Requirements](#requirements)
* [Installing](#installing)
* [Contents](#contents)
* [Usage](#usage)
* [Running Workflow](#running)

### Requirements
* Linux or MacOS
* Diamond
* KronaTools

## Installing
Assembleflow requires a Diamond database. Run the following script to build a protein refseq Diamond database:
```sh
bash ./scripts/download_and_build_diamond_db_with_taxonomy.sh
```

### Contents

The entire workflow has been implemented in NextFlow. An overview of each step in assembleflow:

* Preprocessing reads
	+ Deduplicating reads
	+ Merging overlapping read pairs
	+ Adapter trimming
* Assembly
	+ Assemble reads into contigs
	+ Align reads back to contigs
	+ Find circular contigs
* Assign taxonomy
	+ Assign taxonomy to contigs
	+ Summarize taxonomy for each taxonomic rank


### Usage
Assembleflow uses the `nextflow.config` file to provide parameters to the pipeline. To run the pipeline simply upodate the `nextflow.config` file:
Areas to configure: 
* reads: folder containing the paired-end reads
* adapt: location of the file that contains the adapter sequences in fasta format
* db: location of the database used by Diamond
* scaf_len: contig length filter (reads below this length will be ignored when assigning taxonomy)
* outdir: output directory


### Running
Once your workflow.sh has been configured you can start the workflow by simply running:
`nextflow run main.nf`
