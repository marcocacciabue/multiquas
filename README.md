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
between samples). See [test.csv](https://github.com/marcocacciabue/multiquas/blob/7dbce91bae963665a2db907a9f699dc1e8de66c0/test.csv) for an example.
If your sample files is called "samples.csv", you can then run:
```bash
nextflow run marcocacciabue/multiquas --profile fast --input_csv samples.csv 

```
## What is that profile variable?
Multiquas runs differents reconstructions softwares (Clique, Qure, Viquas, Haploflow, Savage) and with two modes, single  and multiple (using all references). The profiles option is an easy way to define how to run the pipeline. 
```bash
# profile complete will run ALL the reconstructions algorithms for each sample. Most accurate option but it may take
same time.
nextflow run marcocacciabue/multiquas --profile complete --input_csv samples.csv

# profile single will run ALL the reconstrucctions algorithms but only in single mode.
nextflow run marcocacciabue/multiquas --profile single --input_csv samples.csv

# profile multiple will run ALL the reconstrucctions algorithms but only in multiple mode.
nextflow run marcocacciabue/multiquas --profile multiple --input_csv samples.csv

# profile fast will run only the Clique algorithm.
nextflow run marcocacciabue/multiquas --profile multiple --input_csv samples.csv

# if no profile is selected Clique, Qure and Haploflow will be used.
nextflow run marcocacciabue/multiquas --input_csv samples.csv


```
Users can override this behaviour and set on or off specific algorithm. For example: the following adds the Viquas multiple step to the default behaivour
```bash

nextflow run marcocacciabue/multiquas --input_csv samples.csv --viquas_m ON


```

