#!/bin/bash
set -euo pipefail

sras_list=$1
output_folder=$2

# clone the pipeline code
git clone https://github.com/nrminor/sra2samrefiner.git
cd sra2samrefiner

# set $HOME
export HOME=$PWD
export PIXI_HOME=$PWD/.pixi

# install Pixi
curl -fsSL https://pixi.sh/install.sh | bash

# add Pixi to PATH using full install location
export PATH="$HOME/.pixi/bin:$PATH"

# solve the locked environment
pixi install --frozen

# put the environment on the path
export PATH="$HOME/.pixi/envs/default/bin:$PATH"

# run the current accessions
nextflow run . \
--accession_list ${sras_list} \
--ref_gbk ./assets/sars2.gbk \
--ref_fasta ./assets/sars2.fasta \
--end_trim_bases 30 \
--results ${output_folder}

