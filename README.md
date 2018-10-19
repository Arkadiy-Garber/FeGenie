# FeGenie

Nearly all forms of life are dependent on iron for basic cellular maintenance and growth. Microorganisms are often tasked with survival and growth in environments where iron is scantily available. To overcome this challenge, a variety of mechanisms have evolved to allow organisms to compete for whatever small concentration of iron is available. These mechanisms have been widely studied over the past several decades, and while the genetic underpinnings of these mechanisms have been confirmed in many model organisms, the advent of next-generation sequencing and the accumulation of genomic data from novel and uncultivated microorganisms necessitates the development of a systemetized database of genes and gene clusters related to iron utilization. This is the primary focus of this program. We developed a library of profile hidden Markov models (pHMMs) that are representative of most known iron-related cellular functions, including iron acquisition, iron storage, iron gene regulation, magnetosome formation, and iron reduction/oxidation.

The input to the program is simply a folder of FASTA files. The FASTA files should contain contigs from a single genome or a metagenomic assembly. There are two other inputs to this program, which are provided here. These inputs are 1) The pHMM library and 2) a text file containing calibrated bitscore cutoffs for each HMM. The workflow of this program is described in the PDF titled "workflow". it is pfarily straightforward, but if you, the user, have any questions, feel free to shoot me an email (arkadiyg@usc.edu).

Also, a part of this program includes an optional cross-validation of the identified putative iron genes against NCBI's nr database. This option is exercized by simply providing the script with the lcation of the nr database, which should be one large file (~100GB) called "nr.faa". The latest release of this databse can be downloaded by running:

        wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

Make sure you have enough room where you will be downloading it, and that your WiFi is good, otherwise, it may take a very long time!

Two things reagrding this optional cross-validation: first, this step greatly increases the computational load and takes about 100 times longer to complete, compared to running the program without cross-validation. So what would have been a 5 minute analysis of a dozen genomes may take 10 hours, and if you are analyzing large metagenome assemblies, it may take several days to complete. However, the identification of the closest homolog in NCBI to your identified iron genes may be incredbily informative, escpially because our HMM library isn't perfect and false positives are a possibility (as they are with most annotation tools). Second, the part of the algorithm that is dedicated to the cross-validation step is largely untested. So by exercizing this optional parameter you are, in effect, acting as a beta tester for our program. So feel free to start issues on GitHub, or yell at me via email, if there are any snafus with the program or its output when the nr database is provided.  

## Dependencies:

### Diamond
### BLAST
### HMMER
### Prodigal

## Sample command:

    python3 FeGenie.v.5.py -DB HMM-lib/ -bin_dir your/genomes/ -bin_ext fa -contigs_source single -bit HMM-bitcutoffs.txt -d 10 -out fegenie-out -inflation 100 -t 4
