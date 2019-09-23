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

### Quick-start (if you installed using the '-setup.sh' script)
    FeGenie.py -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie

### Quick-start (if you installed using the '-setup_noconda.sh' script)
    ./FeGenie.py -hmm_lib HMM-lib/ -bin_dir /directory/of/bins/ -bin_ext fasta -t 16 -out output_fegenie
HMM-lib directory can be found within FeGenie's main repository
-t 8 means that 8 threads will be used for HMMER and BLAST. If you have less than 16 available on your system, set this number lower (default = 1)



## Citing FeGenie:
FeGenie is developed by Arkadiy I. Garber Kenneth H. Nealson, Akihiro Okamoto, Sean M. McAllister, Clara S. Chan, Roman. A. Barco, and Nancy Merino


This project is still a work in progress, and is involved in a publication currently in preparation: "Garber, A.I., Nealson, K.H., Okamoto, A., McAllister, S.M., Chan, C.S., Barco, B.A., and Merino, N. FeGenie: a new database and tool for identification of iron genes and iron gene neighborhoods in genomes and metagenome assemblies". 

If it was useful for your work, you can cite it as: Garber, A.I., Nealson, K.H., Okamoto, A., McAllister, S.M., Chan, C.S., Barco, B.A., and Merino, N. 2018: FeGenie, GitHub repository: https://github.com/Arkadiy-Garber/FeGenie/.


Please also cite various dependencies used by FeGenie.
