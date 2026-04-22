## FeGenie

Please see the Wiki page for an introduction and tutorial on how to use this tool.

### Citing FeGenie

Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA and Merino N (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. *Front. Microbiol.* 11:37. doi: [10.3389/fmicb.2020.00037](https://www.frontiersin.org/articles/10.3389/fmicb.2020.00037/full)

Special thanks to [Michael Lee](https://github.com/AstrobioMike) for helping to put together the Conda environment for FeGenie. Thanks to [Natasha Pavlovikj](https://github.com/npavlovikj) for creating the Conda recipe for FeGenie. Thanks to [Michał Sitko](https://github.com/note) for creating a Dockerfile for FeGenie. Thanks to [Michelle Hallenbeck](https://microspacebiologist.com/) for helping to modernize the installation process.

### Installation with Mamba / Conda

FeGenie requires several external dependencies, including Python, R packages, HMMER, DIAMOND, BLAST, Prodigal, and MetaBAT2.

Create the environment with:

    mamba create -n fegenie \
      -c conda-forge -c bioconda \
      --strict-channel-priority \
      python=3.10 \
      r-base \
      r-ggplot2 \
      r-stringi \
      r-ggpubr \
      r-reshape \
      r-reshape2 \
      r-tidyverse \
      r-argparse \
      r-ggdendro \
      r-pvclust \
      hmmer \
      diamond \
      prodigal \
      blast \
      metabat2 \
      -y

Activate the environment:

    conda activate fegenie

After activating the environment, set up the required environment variables:

    echo ${CONDA_PREFIX}
    mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
    echo '#!/bin/sh
    export PATH="'$(pwd)':$PATH"
    export rscripts="'$(pwd)'/rscripts"
    export iron_hmms="'$(pwd)'/hmms/iron"' > ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

You can then confirm that the HMM path is available with:

    echo ${iron_hmms}

And check that FeGenie is available with:

    FeGenie.py -h

When you are done using FeGenie, deactivate the environment with:

    conda deactivate

### Installation without Conda

    git clone https://github.com/Arkadiy-Garber/FeGenie.git
    cd FeGenie
    bash setup.sh
    ./FeGenie.py -h

### Quick start

Run FeGenie on a directory of genome bins:

    FeGenie.py -bin_dir /directory/of/bins/ -bin_ext fasta -t 16

The argument for `-bin_ext` should match the filename extension of the FASTA files you want analyzed (for example: `fa`, `fasta`, `fna`).

### Quick start (if installed using the no-Conda setup)

    ./FeGenie.py -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie

The `hmms/iron` directory can be found within the main FeGenie repository.

The `-t` argument sets the number of threads used for HMMER and BLAST. For example, `-t 8` uses 8 threads. If your system has fewer available threads, set this number lower. The default is 1.

## Tutorial (Binder)

FeGenie introductory slideshow:

[Content](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-4/FeGenie%20intro%20and%20tutorial.pdf) | [Video presentation](https://www.youtube.com/watch?v=sp5ZDcHaYOc&t=24s)

FeGenie video tutorial:

[Content](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-4/README.md) | [Video presentation](https://www.youtube.com/watch?v=WV0GAGSD4kc)

To start the tutorial, hit the "launch binder" button below, and follow the commands in "Walkthrough".

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arkadiy-Garber/bvcn-binder-FeGenie/master?urlpath=lab)

(Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the [Binder](https://mybinder.org/) team.)

### Walkthrough

Enter the main FeGenie directory:

    cd FeGenie

Print the FeGenie help menu:

    FeGenie.py -h

Run FeGenie on the test dataset:

    FeGenie.py -bin_dir genomes/ -bin_ext fna -out fegenie_out

Go into the output directory and inspect the output files:

    cd fegenie_out
    less FeGenie-geneSummary-clusters.csv

Run FeGenie on gene calls:

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs

Run FeGenie on gene calls and use a reference database (RefSeq sub-sample) for cross-validation:

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs -ref refseq_db/refseq_nr.sample.faa

### Running with Docker

If running FeGenie with Docker, the only dependency you need installed is Docker itself ([installation guide](https://docs.docker.com/install/)).

With Docker installed, you can run FeGenie like this:

    docker run -it -v $(pwd):/data --env iron_hmms=/data/hmms/iron --env rscripts=/data/rscripts note/fegenie-deps ./FeGenie.py -bin_dir /data/test_dataset -bin_ext txt -out fegenie_out -t $(nproc)

`./FeGenie.py ...` follows the normal non-Dockerized flow of arguments.

Be aware that you need to mount the directories containing the files FeGenie is supposed to read. If you are not familiar with Docker, run the `docker run` command from the directory into which you cloned the FeGenie repository. If all the files you pass to FeGenie are inside this directory and you use relative file paths (for example `hmms/iron`), everything should work as expected.

### Upcoming updates (suggestions welcome via Issues)

1. Ability to accept previously annotated genomes and gene calls  
2. Include Cytochrome 579 (and possible rusticyanin)  
3. Improve delineation between MtrA and MtoA for better resolution of iron reduction vs. iron oxidation  
4. Option to report absolute values for gene counts rather than normalized gene counts  
5. Include option to release all results regardless of whether reporting rules were met  
6. Identification of iron-sulfur proteins
