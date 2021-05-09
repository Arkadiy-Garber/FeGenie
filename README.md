## FeGenie


Please see the Wiki page for introduction and tutorial on how to use this tool.

### Citing FeGenie:

Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA and Merino N (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Front. Microbiol. 11:37. doi: [10.3389/fmicb.2020.00037](https://www.frontiersin.org/articles/10.3389/fmicb.2020.00037/full)

Special thanks to Michael Lee (https://github.com/AstrobioMike) for helping to put together the setup.sh script, which signficantly eases installation.


### Easy Installation (if you have Conda installed)
    conda create -n fegenie -c conda-forge -c bioconda -c defaults fegenie=1.0 --yes
    conda activate fegenie
    FeGenie.py -h
    
    conda deactivate # when you are done using FeGenie and would like to deactivate the Conda environment for FeGenie

### Installation (if you don't have Conda)
    git clone https://github.com/Arkadiy-Garber/FeGenie.git
    cd FeGenie
    bash setup_noconda.sh
    ./FeGenie.py -h

### Quick-start
    FeGenie.py -bin_dir /directory/of/bins/ -bin_ext fasta -t 16
The argument for -bin_ext needs to represent the filename extension of the FASTA files in the selected directory that you would like analyzed (e.g. fa, fasta, fna, etc).


### Quick-start (if you installed using the 'setup_noconda.sh' script)
    ./FeGenie.py -hmm_lib hmms/iron -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie
`hmms/iron` directory can be found within FeGenie's main repository
-t 8 means that 8 threads will be used for HMMER and BLAST. If you have less than 16 available on your system, set this number lower (default = 1)

## Tutorial (Binder)

FeGenie introductory slideshow:

[Content](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-4/FeGenie%20intro%20and%20tutorial.pdf) | [Video presentation](https://www.youtube.com/watch?v=sp5ZDcHaYOc&t=24s)


FeGenie video tutorial:

[Content](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-4/README.md) | [Video presentation](https://www.youtube.com/watch?v=WV0GAGSD4kc)


To start the tutorial, hit the 'launch binder' button below, and follow the commands in 'Walkthrough'

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arkadiy-Garber/bvcn-binder-FeGenie/master?urlpath=lab)
(Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!)


### Walkthrough

Enter the main FeGenie directory

    cd FeGenie

print the FeGenie help menu

    FeGenie -h

run FeGenie on test dataset

    FeGenie.py -bin_dir genomes/ -bin_ext fna -out fegenie_out

Go into the output directory and check out the output files

    cd fegenie_out
    less FeGenie-geneSummary-clusters.csv

run FeGenie on gene calls

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs

run FeGenie on gene calls, and use reference database (RefSeq sub-sample) for cross-validation

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs -ref refseq_db/refseq_nr.sample.faa


### Running with docker

In case of running `FeGenie` with docker the only dependency you need to have installed is docker itself ([installation guide](https://docs.docker.com/install/)).

With docker installed you can run `FeGenie` in the following way:

    docker run -it -v $(pwd):/data note/fegenie-deps ./FeGenie.py -bin_dir test_dataset -bin_ext txt -out fegenie_out -hmm_lib hmms/iron -t $(nproc)

`./FeGenie.py ...` follows normal, non-dockerized flow of arguments.

Beware that you need to mount directories which contain files `FeGenie` is supposed to read. If you are not familiar with docker then run `docker run` command from the directory into which you cloned `FeGenie` repository. If all the files you pass to `FeGenie` are in inside this directory and you use relative filepaths (like e.g. `hmms/iron`) everything will work just fine.

### Upcoming Updates (we welcome more suggestions, which can be submitted as an Issue)
1) Ability to accept previously-annotated genomes and gene-calls.
2) Include Cytochrome 579 (and possible rusticyanin)
3) Improve dilineation between MtrA and MtoA for better resolution with respect to identification of iron reduction and iron oxidation, respectively.
5) Option to report absolute values for gene counts (rather than normalized gene counts)
6) Include option to release all results (regardless of whether rules for reporting were met)
7) Identification of iron-sulfur proteins.

