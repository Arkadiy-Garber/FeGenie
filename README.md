## FeGenie

Please see the Wiki page for introduction and tutorial on how to use this tool.

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
2) Identification of iron-sulfur proteins.
3) Inlcude Cytochrome 579.
4) Improve dilineation between MtrA and MtoA for better resolution with respect to identification of iron reduction and iron oxidation, respectively.

### Citing FeGenie:
This project is involved in a publication currently in review: "Garber, A.I., Nealson, K.H., Okamoto, A., McAllister, S.M., Chan, C.S., Barco, B.A., and Merino, N. FeGenie: a comprehensive tool for the identification of iron genes and iron gene neighborhoods in genomes and metagenome assemblies".

If it was useful for your work, you can cite the bioRxiv pre-print: https://www.biorxiv.org/content/10.1101/777656v1


Please also cite various dependencies used by FeGenie.
