# MultiQuas

[![](https://img.shields.io/badge/nextflow-25.10.7-green)](https://www.nextflow.io)
[![](https://img.shields.io/badge/docker-29.1.3-blue)](https://docs.docker.com/get-docker)

## Why Multiquas?

MultiQuas is all in one pipeline for the reconstruction of viral
quasispecies population from short read data. It comes with *multiple*
reconstructions algorithms and the user can use *multiple* references
(ideally sequences expected to be present in the population) to produce
better results.  
Also, the pipeline includes an evaluation step that gives an idea of how
plausible the quasispecies is for each of the reconstruction algorithms.
The higher the R-squared and the closer to 1 the slop are the more
plausible the reconstruction is considered. Nonetheless, users are
recommended to check and select the quasispecies that better represent
the underlying complexity.

## Requirements

Requires [nextflow](https://www.nextflow.io) and
[docker](https://docs.docker.com/get-docker). Please install those
dependencies.

## No time to waste

If you just want to see the tool in action run the following:

``` bash
# clone (or manually download and unpack if git is not installed) the repo like this:

git clone https://github.com/marcocacciabue/multiquas.git

# go to the directory

cd multiquas
# run the pipeline with a test dataset.
nextflow run main.nf -resume -profile fast --input_csv test.csv 
```

This will download all the docker images needed and then run the
pipeline on a test sample (could take some time, only the first time
around). If you wish to run a specific version of nextflow you can use
the on the fly variable. For example to run nextflow version 25.10.7:

``` bash

NXF_VER=25.10.7 nextflow run main.nf -resume -profile fast --input_csv test.csv 
```

## I have my samples, how do I run the pipeline on them?

You need to have the following in the working directory: - a reference
file (fasta). - sample reads files (fastq). - a cvs file with four
columns (sample_id,fastq_1,fastq_2,input_ref).

Each line in the cvs file is a different sample of name sample_id.
fastq_1 and fastq_2 are the relative paths for the corresponding reads.
input_ref is the reference to use (can be different between samples).
See
[test.csv](https://github.com/marcocacciabue/multiquas/blob/7dbce91bae963665a2db907a9f699dc1e8de66c0/test.csv)
for an example. If your sample files is called “samples.csv”, you can
then run:

``` bash
nextflow run marcocacciabue/multiquas -resume -profile fast --input_csv samples.csv 
```

Note: this command is different in the way that it will clone the repo
directly to a nextflow folder for you. It is the same \## What is that
profile variable? Multiquas runs differents reconstructions softwares
(Clique, Qure, Viquas, Haploflow, Savage) and with two modes, single and
multiple (using all references). The profiles option is an easy way to
define how to run the pipeline.

``` bash
# profile complete will run ALL the reconstructions algorithms for each sample. Most accurate option but it may take
#same time.
nextflow run marcocacciabue/multiquas -resume -profile complete --input_csv samples.csv

# profile single will run ALL the reconstrucctions algorithms but only in single mode.
nextflow run marcocacciabue/multiquas -resume -profile single --input_csv samples.csv

# profile multiple will run ALL the reconstrucctions algorithms but only in multiple mode.
nextflow run marcocacciabue/multiquas -resume -profile multiple --input_csv samples.csv

# profile fast will run only the Clique algorithm.
nextflow run marcocacciabue/multiquas -resume -profile multiple --input_csv samples.csv

# if no profile is selected Clique, Qure and Haploflow will be used.
nextflow run marcocacciabue/multiquas -resume --input_csv samples.csv
```

Users can override this behavior and set on or off specific algorithms.
For example: the following adds the Viquas multiple step to the default
behaivour

``` bash

nextflow run marcocacciabue/multiquas -resume --input_csv samples.csv --viquas_m ON
```

## Reconstruction algorithms

Multiquas includes a set of several reconstructions algorithms that may
increase over time (if there is some interest in it). Rigth now it comes
with:

- [CliqueSNV](https://github.com/vtsyvina/CliqueSNV) Knyazev *et al.*
  (2021).
  [![](https://zenodo.org/badge/DOI/10.1093/nar/gkab576.svg)](https://doi.org/10.1093/nar/gkab576)

- [QuRe](https://sourceforge.net/projects/qure/) Prosperi & Salemi
  (2012).
  [![](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btr627.svg)](https://doi.org/10.1093/bioinformatics/btr627)

- [ViQuaS](http://sourceforge.net/projects/viquas/) Jayasundara *et al.*
  (2015).
  [![](https://zenodo.org/badge/DOI/10.1093/nar/gkab576.svg)](https://doi.org/10.1093/nar/gkab576)

- [Haploflow](https://github.com/hzi-bifo/Haploflow) Fritz *et al.*
  (2021).
  [![](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btu754.svg)](https://doi.org/10.1093/bioinformatics/btu754)

- [Savage](https://github.com/HaploConduct/HaploConduct) Baaijens &
  Schönhuth (2019).
  [![](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btz255.svg)](https://doi.org/10.1093/bioinformatics/btz2554)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-baaijens_overlap_2019" class="csl-entry">

Baaijens J.A. & Schönhuth A. (2019). Overlap graph-based generation of
haplotigs for diploids and polyploids. Bioinformatics 35 (21):
4281–4289. <https://doi.org/10.1093/bioinformatics/btz255>.

</div>

<div id="ref-fritz_haploflow_2021" class="csl-entry">

Fritz A., Bremges A., Deng Z.-L., Lesker T.R., Götting J., Ganzenmueller
T., Sczyrba A., Dilthey A., Klawonn F. & McHardy A.C. (2021). Haploflow:
Strain-resolved de novo assembly of viral genomes. Genome Biology 22
(1): 212. <https://doi.org/10.1186/s13059-021-02426-8>.

</div>

<div id="ref-jayasundara_viquas_2015" class="csl-entry">

Jayasundara D., Saeed I., Maheswararajah S., Chang B.C., Tang S.-L. &
Halgamuge S.K. (2015). ViQuaS: An improved reconstruction pipeline for
viral quasispecies spectra generated by next-generation sequencing.
Bioinformatics 31 (6): 886–896.
<https://doi.org/10.1093/bioinformatics/btu754>.

</div>

<div id="ref-knyazev_accurate_2021" class="csl-entry">

Knyazev S., Tsyvina V., Shankar A., Melnyk A., Artyomenko A., Malygina
T., Porozov Y.B., Campbell E.M., Switzer W.M., Skums P., Mangul S. &
Zelikovsky A. (2021). Accurate assembly of minority viral haplotypes
from next-generation sequencing through efficient noise reduction.
Nucleic Acids Research 49 (17): e102–e102.
<https://doi.org/10.1093/nar/gkab576>.

</div>

<div id="ref-prosperi_qure_2012" class="csl-entry">

Prosperi M.C.F. & Salemi M. (2012). QuRe: Software for viral
quasispecies reconstruction from next-generation sequencing data.
Bioinformatics 28 (1): 132–133.
<https://doi.org/10.1093/bioinformatics/btr627>.

</div>

</div>
