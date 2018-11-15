# FeGenie

Please see the Wiki page for introduction and tutorial on how to use this tool.

## Dependencies:

### Python (version 3.6 or higher)
### Diamond (only necessary if you are doing the cross-validation against nr)
### BLAST
### HMMER
### Prodigal
###

## Obtaining NCBI's nr database for cross-validation

Part of this program includes an optional cross-validation of the identified putative iron genes against NCBI's nr database. This option is exercized by simply providing the script with the lcation of the nr database, which should be one large file (~100GB) called "nr.faa". The latest release of this databse can be downloaded by running:

        wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

Make sure you have enough room where you will be downloading it, and that your WiFi is good, otherwise, it may take a very long time! And really, you can provide the program with any reference file of proteins in FASTA format, and the program will blast the identified iron genes against that.

Two things regarding this optional cross-validation: first, this step greatly increases the computational load and takes about 100 times longer to complete, compared to running the program without cross-validation. So what would have been a 5 minute analysis of a dozen genomes may take 10 hours, and if you are analyzing large metagenome assemblies, it may take several days to complete. However, the identification of the closest homolog in NCBI to your identified iron genes may be incredbily informative, escpially because our HMM library isn't perfect and false positives are a possibility (as they are with most annotation tools). Second, the part of the algorithm that is dedicated to the cross-validation step is largely untested. So by exercizing this optional parameter you are, in effect, acting as a beta tester for our program. So feel free to start issues on GitHub, or yell at me via email, if there are any snafus with the program or its output when the nr database is provided.

# Citing taxonsluice:
FeGenie is developed by Nancy Merino, Arkadiy Garber, and Kenneth Nealson, University of Southern California, Los Angeles, CA, USA.

This project is still a work in progress, and is involved in a publication currently in preparation: "Garber, A.I., Nealson, K.H., Merino, N. FeGenie: a new database and tool for identification of iron genes and iron gene clusters in genomes and metagenome assemblies". If it was useful for your work, you can cite it as: Garber, A.I., and Merino, N. and Nealson, K.H. 2018: FeGenie, GitHub repository: https://github.com/Arkadiy-Garber/FeGenie/.


Please also cite various dependencies used by FeGenie.
