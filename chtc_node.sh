#!/bin/bash

# clone the pipeline code
git clone https://github.com/nrminor/sra2samrefiner.git .

# get pixi
curl -fsSL https://pixi.sh/install.sh | sh

# solve the locked environment
pixi shell --frozen

# modify path
PATH=$PATH:$(pwd)/.pixi/bin

# run the current accessions
nextflow run . \
--accession_list sra_queue.txt \
--ref_gbk ./assets/sars2.gbk \
--ref_fasta ./assets/sars2.fasta \
--end_trim_bases 30 \
-resume
