# Multiquas

[![](https://img.shields.io/badge/nextflow-25.10.7-green)](https://www.nextflow.io) [![](https://img.shields.io/badge/docker-29.1.3-blue)](https://docs.docker.com/get-docker)

## Install

Requires [nextflow](https://www.nextflow.io) and [docker](https://docs.docker.com/get-docker)) installed. Please install those dependencies.

## No time to waste
If you just want to see the tool in action run the following:
```bash

nextflow run marcocacciabue/multiquas -latest --input_csv data/test.csv 

```
This will clone this repo and all the docker images needed to run the pipeline on a test sample (could
take some time, only the first time around). 