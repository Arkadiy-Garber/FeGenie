## FeGenie

Please see the Wiki page for introduction and tutorial on how to use this tool.

Special thanks to Michael Lee (https://github.com/AstrobioMike) for helping to put together the setup.sh script, which signficantly eases installation.

### Easy Installation (if you have Conda installed)
    git clone https://github.com/Arkadiy-Garber/FeGenie.git
    cd FeGenie
    ./setup.sh
    conda activate fegenie
    FeGenie.py -h

### Easy Installation (if you don't have Conda)
    git clone https://github.com/Arkadiy-Garber/FeGenie.git
    cd FeGenie
    ./setup_noconda.sh
    ./FeGenie.py -h

### Quick-start (if you installed using the 'setup.sh' script)
    FeGenie.py -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie
The argument for -bin_ext needs to represent the filename extension of the FASTA files in the selected directory that you would like analyzed (e.g. fa, fasta, fna, etc).

### Quick-start (if you installed using the 'setup_noconda.sh' script)
    ./FeGenie.py -hmm_lib HMM-lib/ -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie
HMM-lib directory can be found within FeGenie's main repository
-t 8 means that 8 threads will be used for HMMER and BLAST. If you have less than 16 available on your system, set this number lower (default = 1)

### Upcoming Updates (we welcome more suggestions, which can be submitted as an Issue)
1) Ability to accept previously-annotated genomes and gene-calls.
2) Include Cytochrome 579 (and possible rusticyanin)
3) Improve dilineation between MtrA and MtoA for better resolution with respect to identification of iron reduction and iron oxidation, respectively.
5) Option to report absolute values for gene counts (rather than normalized gene counts)
6) Include option to release all results (regardless of whether rules for reporting were met)
7) Identification of iron-sulfur proteins.

### Citing FeGenie:

Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA and Merino N (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Front. Microbiol. 11:37. doi: 10.3389/fmicb.2020.00037
