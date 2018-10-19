# FeGenie

Nearly all forms of life are dependent on iron for basic cellular maintenance and growth. Microorganisms are often tasked with survival and growth in environments where iron is scantily available. To overcome this challenge, a variety of mechanisms have evolved to allow organisms to compete for whatever small concentration of iron is available. These mechanisms have been widely studied over the past several decades, and while the genetic underpinnings of these mechanisms have been confirmed in many model organisms, the advent of next-generation sequencing and the accumulation of genomic data from novel and uncultivated microorganisms necessitates the development of a systemetized database of genes and gene clusters related to iron utilization. This is the primary focus of this program. We developed a library of profile hidden Markov models (pHMMs) that are representative of most known iron-related cellular functions, including iron acquisition, iron storage, iron gene regulation, magnetosome formation, and iron reduction/oxidation.

The input to the program is simply a folder of FASTA files. The FASTA files should contain contigs from a single genome or a metagenomic assembly. There are two other inputs to this program, which are provided here. These inputs are 1) The pHMM library and 2) a text file containing calibrated bitscore cutoffs for each HMM. The workflow of this program is described in the PDF titled "workflow". it is pfarily straightforward, but if you, the user, have any questions, feel free to shoot me an email (arkadiyg@usc.edu).

## Dependencies:

### Diamond
### BLAST
### HMMER
### Prodigal

## Sample command:

    python3 FeGenie.py
