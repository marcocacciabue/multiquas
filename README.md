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

