# Multiquas

[![](https://img.shields.io/badge/nextflow-25.10.7-green)](https://www.nextflow.io) [![](https://img.shields.io/badge/docker-29.1.3-blue)](https://docs.docker.com/get-docker)

## Install

Requires [nextflow](https://www.nextflow.io) and [docker](https://docs.docker.com/get-docker)) installed. Please install those dependencies.

## No time to waste
If you just want to see the tool in action run the following:

```bash
# clone (or manually download and unpack if git is not installed) the repo like this:

git clone https://github.com/marcocacciabue/multiquas.git

# go to the directory

cd multiquas
# run the pipeline with a test dataset.
nextflow run main.nf --profile fast --input_csv test.csv 

```
This will download all the docker images needed and then run the pipeline on a test sample (could
take some time, only the first time around). 
If you wish to run a specific version of nextflow you can use the on the fly variable. For example
to run nextflow version 25.10.7:

```bash

NXF_VER=25.10.7 nextflow run main.nf --profile fast --input_csv test.csv 

```

## I have my samples, how do I run the pipeline on them?

You need to have the following in the working directory: 
- a reference file (fasta). 
- sample reads files (fastq).
- a cvs file with four columns (sample_id,fastq_1,fastq_2,input_ref).

Each line in the cvs file is a different sample of name sample_id. fastq_1 and fastq_2
are the paths for the corresponding reads. input_ref is the reference to use (can be different
between samples).






