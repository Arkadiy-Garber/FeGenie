#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

gzip -d test_dataset.tar.gz
tar xf test_dataset.tar
rm test_dataset.tar

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

## creating GToTree environment and installing dependencies
conda create -n fegenie hmmer diamond prodigal blast --yes

## activating environment
source activate fegenie

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding FeGenie bin path and HMM_dir variable:
echo '#!/bin/sh'" \

export PATH=\"$(pwd):"'$PATH'\"" \

export rscripts=\"$(pwd)/rscripts\"

export iron_hmms=\"$(pwd)/hmms/iron\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
source activate fegenie


printf "\n        ${GREEN}DONE!${NC}\n\n"
