#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD CYTOCHROME 579 HMM (next release)
# TODO: ADD BITSCORE CUTOFF FOR EACH HMM HIT


def unique(ls, ls2):
    unqlist = []
    for i in ls:
        if i not in unqlist and i in ls2:
            unqlist.append(i)
    return len(unqlist)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def remove2(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    # outString = "".join(emptyList)
    return emptyList


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="FeGenie.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber and Nancy Merino;
    University of Southern California, Earth Sciences
    Please send comments and inquiries to arkadiyg@usc.edu

    *******************************************************
    '''))

parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

parser.add_argument('-contigs_source', type=str, help="are the provided contigs from a single organism (single)"
                                                     "or are you providing this program with metagenomic/metatranscriptomic assemblies (meta)? "
                                                     "(default=single)", default="single")

parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                         "This number should be an integer and should reflect the maximum number of "
                                         "genes in between putative iron-related genes identified by the HMM database "
                                         "(default=10)", default=10)

parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format", default="NA")

parser.add_argument('-out', type=str, help="name of output file; please provide full path (default=fegenie_out)", default="fegenie_out")

parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                    default=1000)

parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                         "(default=1, max=16)", default=1)

parser.add_argument('-makeplots', type=str, help="Would you like FeGenie to make some figures from your data? y = yes, n = no (default = n). "
                                                 "If so, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                                                 "Warning: this part of the program is currently under beta-testing, and if there are any problems running Rscript, "
                                                 "or installing any of the required packages, you may get a bunch of error messages at the end. "
                                                 "This will not crash the program, however, and you can still expect to "
                                                 "see the main output (CSV files) from FeGenie.", default="n")

# CHECKING FOR CONDA INSTALL
os.system("echo ${HMM_dir}/HMM-bitcutoffs.txt > HMMlib.txt")
file = open("HMMlib.txt")
for i in file:
    location = i.rstrip()

try:
    bits = open(location)
    conda = 1
except FileNotFoundError:
    conda = 0

if conda == 0:
    parser.add_argument('-hmm_lib', type=str, help='HMM database; directory titled \'HMM-lib\', can be found in the FeGenie folder', default="NA")

    parser.add_argument('-R', type=str,
                        help="location of R scripts directory (note: this optional argument requires Rscript to be "
                             "installed on your system). The R scripts directory is in the same directory as the "
                             "FeGenie python code", default="NA")

args = parser.parse_args()


# ************** Checking for the required arguments ******************* #
cwd = os.getcwd()
print("checking arguments")
if conda == 0:
    if args.hmm_lib != "NA":
        print(".")
    else:
        print("You have not provided the location of the HMM library via the -hmm_lib argument. Please do so, and try "
              "again. The HMM library is found within the same directory as the FeGenie executable.")
        print("Exiting")
        raise SystemExit

if args.bin_dir != "NA":
    binDir = args.bin_dir + "/"
    binDirLS = os.listdir(args.bin_dir)
    print(".")
else:
    print("Looks like you did not provide a directory of genomes/bins or assemblies.")
    print("Exiting")
    raise SystemExit

if args.makeplots == 'y':
    if args.R != "NA":
        print(".")
    else:
        if conda == 0:
            print('Looks like you told FeGenie to automatically generate R plots. '
                  'However, you have not provided the location of the directory that contains the R scripts '
                  '(as required of you because you did not go through the conda-based installation.')
            print("Exiting")
            raise SystemExit


if args.bin_ext != "NA":
    print(".")
else:
    print('Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
          ' which files in the provided directory are fasta files that you would like analyzed.')
    print("Exiting")
    raise SystemExit

try:
    os.listdir(args.out)
    print("Looks like you already have a directory with the name: " + args.out)
    print("To avoid overwriting potentially valuable files, FeGenie will now exit. "
          "Please delete or rename said directory and try running again.")
    print("Exiting")
    raise SystemExit
except FileNotFoundError:
    print(".")
    os.system("mkdir %s" % args.out)
    outDirectory = "%s" % args.out
    outDirectoryLS = os.listdir("%s" % args.out)

print("All required arguments provided!")
print("")


# *************** MAKE NR A DIAMOND DB AND READ THE FILE INTO HASH MEMORY ************************ #
if args.ref != "NA":
    try:
        testFile = open(args.ref + ".dmnd")

    except FileNotFoundError:
        print("Making diamond database out of provided reference file")
        os.system("diamond makedb --in %s -d %s" % (args.ref, args.ref))

    ref = open(args.ref, "r")
    refDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    print("Reading reference file into memory")
    for i in ref:
        if re.match('>', i):
            accession = i.rstrip().split(" ")[0][1:]
            header = replace(i.rstrip(), [","], ";")
            header = header[1:]
            refDict[accession] = header
else:
    pass


# *************** CALL ORFS FROM BINS AND READ THE ORFS INTO HASH MEMORY ************************ #
BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in binDirLS:
    if lastItem(i.split(".")) == args.bin_ext:
        cell = i
        try:
            testFile = open("%s%s-proteins.faa" % (binDir, i), "r")
            print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))

        except FileNotFoundError:
            print("Finding ORFs for " + cell)
            if args.contigs_source == "single":
                os.system("prodigal -i %s%s -a %s%s-proteins.faa -o %s%s-prodigal.out -q" % (binDir, i, binDir, i, binDir, i))
            elif args.contigs_source == "meta":
                os.system("prodigal -i %s%s -a %s%s-proteins.faa -o %s%s-prodigal.out -p meta -q" % (binDir, i, binDir, i, binDir, i))
            else:
                print("WARNING: you did not specify whether the provided FASTA files are single genomes or "
                      "metagenome/metatranscriptome assemblies. By default, FeGenie is assuming that these are "
                      "single genomes, and running Prodigal accordingly. Just an FYI.")
                os.system("prodigal -i %s%s -a %s%s-proteins.faa -o %s%s-prodigal.out -q" % (binDir, i, binDir, i, binDir, i))

        file = open(binDir + i + "-proteins.faa", "r")
        file = fasta(file)
        for j in file.keys():
            orf = j.split(" # ")[0]
            BinDict[cell][orf] = file[j]


# ******************** READ BITSCORE CUT-OFFS INTO HASH MEMORY ****************************** #
if conda == 1:
    os.system("echo ${HMM_dir} > HMMlib.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = (i.rstrip())
    os.system("rm HMMlib.txt")
else:
    HMMdir = args.hmm_lib


HMMbits = HMMdir + "/HMM-bitcutoffs.txt"
meta = open(HMMbits, "r")
metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in meta:
    ls = i.rstrip().split("\t")
    metaDict[ls[0]] = ls[1]

# ******************* BEGINNING MAIN ALGORITHM **********************************))))
print("starting main pipeline...")
HMMdirLS = os.listdir(HMMdir)
for FeCategory in HMMdirLS:
    if not re.match(r'\.', FeCategory) and FeCategory != "HMM-bitcutoffs.txt":
        print("")
        print(".")
        print("Looking for following iron-related functional category: " + FeCategory)
        hmmDir = "%s/%s/" % (HMMdir, FeCategory)
        hmmDirLS2 = os.listdir("%s/%s" % (HMMdir, FeCategory))

        HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
        for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
            if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                # print("analyzing: " + i)
                FASTA = open(binDir + i + "-proteins.faa", "r")
                FASTA = fasta(FASTA)
                os.system(
                    "mkdir " + binDir + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

                print("")
                count = 0
                for hmm in hmmDirLS2:  # ITERATING THROUGH ALL THE HMM FILES IN THE HMM DIRECTORY
                    count += 1
                    perc = (count / len(hmmDirLS2)) * 100
                    sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                    sys.stdout.flush()
                    if len(metaDict[hmm.split(".")[0]]) == 0:
                        bit = 0
                    else:
                        bit = metaDict[hmm.split(".")[0]]


                        # print("performing HMMSEARCH")
                        os.system(
                            "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/%s-proteins.faa"
                            % (int(args.t), float(bit), binDir, i, hmm, binDir, i, hmm, hmmDir, hmm, binDir, i)
                        )

                        # REMOVING THE STANDARD OUTPUT FILE
                        os.system(
                            "rm " + binDir + "/" + i + "-HMM/" + hmm + ".txt"
                        )

                        # READING IN THE HMMSEARCH RESULTS (TBLOUT) OUT FILE
                        hmmout = open(binDir + i + "-HMM/" + hmm + ".tblout", "r")

                        # COLLECTING SIGNIFICANT HMM HITS IN THE FILE
                        for line in hmmout:
                            if not re.match(r'#', line):
                                ls = delim(line)
                                evalue = float(ls[4])
                                bit = float(ls[5])
                                orf = ls[0]
                                if evalue < float(1E-1):  # FILTERING OUT BACKGROUND NOISE
                                    # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS

                                    if orf not in HMMdict[i]:
                                        HMMdict[i][orf]["hmm"] = hmm
                                        HMMdict[i][orf]["evalue"] = evalue
                                        HMMdict[i][orf]["bit"] = bit
                                    else:
                                        # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                                        if bit > HMMdict[i][orf]["bit"]:
                                            HMMdict[i][orf]["hmm"] = hmm
                                            HMMdict[i][orf]["evalue"] = evalue
                                            HMMdict[i][orf]["bit"] = bit

                os.system("rm -rf " + binDir + "/" + i + "-HMM")

        out = open(outDirectory + "/%s-summary.csv" % (FeCategory), "w")
        out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "\n")
        for key in HMMdict.keys():
            for j in HMMdict[key]:
                out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                          str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "\n")

        out.close()

print("\n")
print("Consolidating summary files into one master summary file")
out = open(outDirectory + "/FinalSummary.csv", "w")
if args.ref == "NA":
    out.write("category" + "," + "cell" + "," + "orf" + "," + "related_hmm" + "," + "HMM-bitscore" + "\n")

resultsDir = os.listdir(outDirectory)
for i in resultsDir:
    if lastItem(i.split("-")) == "summary.csv":
        result = open(outDirectory + "/" + i, "r")
        for j in result:
            ls = j.rstrip().split(",")
            cell = ls[0]
            orf = ls[1]
            hmm = ls[2]
            bit = ls[4]

            if cell != "cell":
                out.write(i.split("-summary")[0] + "," + cell + "," + orf + "," + hmm + "," + str(bit) + "\n")

out.close()


# ****************************************** DEREPLICATION *********************************************************
summary = open(outDirectory + "/FinalSummary.csv", "r")
SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
for i in summary:
    ls = i.rstrip().split(",")
    if ls[0] != "category" and ls[0] != "FeGenie":
        if len(ls) > 0:
            category = ls[0]
            cell = ls[1]
            orf = ls[2]
            hmm = ls[3]
            hmmBit = ls[4]

            if cell not in SummaryDict.keys():
                SummaryDict[cell][orf]["hmm"] = hmm
                SummaryDict[cell][orf]["hmmBit"] = hmmBit
                SummaryDict[cell][orf]["category"] = category

            else:
                if orf not in SummaryDict[cell]:
                    SummaryDict[cell][orf]["hmm"] = hmm
                    SummaryDict[cell][orf]["hmmBit"] = hmmBit
                    SummaryDict[cell][orf]["category"] = category

                else:
                    if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                        SummaryDict[cell][orf]["hmm"] = hmm
                        SummaryDict[cell][orf]["hmmBit"] = hmmBit
                        SummaryDict[cell][orf]["category"] = category


# ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
print("Identifying genomic proximities and putative operons")
CoordDict = defaultdict(lambda: defaultdict(list))
for i in SummaryDict.keys():
    if i != "category":
        for j in SummaryDict[i]:
            contig = allButTheLast(j, "_")
            numOrf = lastItem(j.split("_"))
            CoordDict[i][contig].append(int(numOrf))

counter = 0
print("Clustering ORFs...")
print("")
out = open(outDirectory + "/FinalSummary-dereplicated-clustered.csv", "w")
for i in CoordDict.keys():
    print(".")
    for j in CoordDict[i]:
        LS = (CoordDict[i][j])
        clusters = (cluster(LS, args.d))
        for k in clusters:
            if len(RemoveDuplicates(k)) == 1:
                orf = j + "_" + str(k[0])

                if args.ref != "NA":
                    ncbiHomolog = SummaryDict[i][orf]["NCBImatch"]
                    if ncbiHomolog != "NA":
                        ncbiHomolog = ncbiHomolog.split("]")[0] + "]"

                    out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," +
                              str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "," + ncbiHomolog + "," + str(SummaryDict[i][orf]["NCBIeval"]) + "\n")

                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1

                else:
                    out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1

            else:
                if args.ref != "NA":
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)
                        ncbiHomolog = SummaryDict[i][orf]["NCBImatch"]
                        if ncbiHomolog != "NA":
                            ncbiHomolog = ncbiHomolog.split("]")[0] + "]"

                        out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"]
                            + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "," + ncbiHomolog + "," + str(SummaryDict[i][orf]["NCBIeval"]) + "\n")

                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1
                else:
                    for l in RemoveDuplicates(k):
                        orf = j + "_" + str(l)

                        out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1
out.close()


# ************************** BLAST-BASED METHODS/LOOKING FOR UNMODELED MARKERS ********************************
thermincola = "%s/iron_reduction/non-aligned/TherJR_SLCs.faa" % HMMdir
geobacter = "%s/iron_reduction/non-aligned/geobacter_PCCs.faa" % HMMdir

print("Looking for Thermincola S-layer cytochromes and Geobacter porin-cytochromes")

for i in binDirLS:
    if lastItem(i.split(".")) == args.bin_ext:
        os.system(
            "makeblastdb -dbtype prot -in %s/%s-proteins.faa -out %s/%s-proteins.faa -logfile %s/makedbfile.txt" % (binDir, i, binDir, i, binDir))
        os.system("rm %s/makedbfile.txt" % binDir)

        os.system(
            "blastp -query %s -db %s/%s-proteins.faa -num_threads %s -outfmt 6 -out %s/%s-thermincola.blast -evalue 1E-10"
            % (thermincola, binDir, i, args.t, outDirectory, i))

        os.system(
            "blastp -query %s -db %s/%s-proteins.faa -num_threads %s -outfmt 6 -out %s/%s-geobacter.blast -evalue 1E-10"
            % (geobacter, binDir, i, args.t, outDirectory, i))


geoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
geo = open("%s/iron_reduction/non-aligned/geobacter_PCCs.faa" % HMMdir)
geo = fasta(geo)
for i in geo.keys():
    id = i.split(" ")[0]
    type = (i.split(" ")[2])
    type = type[1:len(type) - 1]
    geoDict[id]["type"] = type
    geoDict[id]["header"] = i

thermDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
therm = open("%s/iron_reduction/non-aligned/TherJR_SLCs.faa" % HMMdir)
therm = fasta(therm)
for i in therm.keys():
    id = i.split(" ")[0]
    thermDict[id]["header"] = i


out = open(outDirectory + "/GeoThermin.csv", "w")
for blastresult in os.listdir(args.out):
    blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    if re.findall(r'geobacter\.blast', blastresult):
        blast = open(outDirectory + "/" + blastresult, "r")

        for i in blast:
            ls = i.rstrip().split("\t")
            if ls[1] not in blastDict.keys():
                blastDict[ls[1]]["query"] = ls[0]
                blastDict[ls[1]]["e"] = ls[10]
            else:
                if float(ls[10]) < float(blastDict[ls[1]]["e"]):
                    blastDict[ls[1]]["query"] = ls[0]
                    blastDict[ls[1]]["e"] = ls[10]
                else:
                    pass

        operonDict = defaultdict(list)
        for i in blastDict.keys():
            orf = (i)
            contig = allButTheLast(i, "_")
            orfCall = lastItem(orf.split("_"))
            id = (blastDict[i]["query"])
            type = geoDict[id]["type"]
            operonDict[contig].append(int(orfCall))

        count = 0
        operonDict2 = defaultdict(lambda: defaultdict(list))
        for i in operonDict.keys():
            contig = (i)
            clu = cluster(operonDict[i], 2)
            for j in clu:
                if len(j) > 1:
                    operon = (j)
                    count += 1
                    for k in operon:
                        orf = str(contig) + "_" + str(k)
                        id = blastDict[orf]["query"]
                        header = geoDict[id]["header"]
                        type = geoDict[id]["type"]
                        operonDict2["operon" + str(count)]["types"].append(type)
                        operonDict2["operon" + str(count)]["orfs"].append(orf)
                        operonDict2["operon" + str(count)]["headers"].append(header)

        for i in operonDict2.keys():
            if "porin" in operonDict2[i]["types"] and (
                    "pc" in operonDict2[i]["types"] or "omc" in operonDict2[i]["types"]):
                genome = blastresult.split("-geobacter.blas")[0]
                category = "iron_reduction"
                for j in range(0, len(operonDict2[i]["types"])):
                    orf = (operonDict2[i]["orfs"][j])
                    evalue = blastDict[orf]["e"]
                    header = (operonDict2[i]["headers"][j])
                    out.write(category + "," + genome + "," + orf + "," + replace(header, [","], ";") + "," + "evalue: " + str(evalue) + "," + str(counter) + "\n")
                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                counter += 1
                # print("\n\n")

    if re.findall(r'thermincola\.blast', blastresult):
        blast = open(outDirectory + "/" + blastresult, "r")

        for i in blast:
            ls = i.rstrip().split("\t")
            if ls[1] not in blastDict.keys():
                blastDict[ls[1]]["query"] = ls[0]
                blastDict[ls[1]]["e"] = ls[10]
            else:
                if float(ls[10]) < float(blastDict[ls[1]]["e"]):
                    blastDict[ls[1]]["query"] = ls[0]
                    blastDict[ls[1]]["e"] = ls[10]
                else:
                    pass

        operonDict = defaultdict(lambda: defaultdict(list))
        for i in blastDict.keys():
            orf = (i)
            contig = allButTheLast(i, "_")
            id = (blastDict[i]["query"])
            header = thermDict[id]["header"]
            operonDict[blastresult]["headers"].append(header)
            operonDict[blastresult]["orfs"].append(orf)

        count = 0
        operonDict2 = defaultdict(lambda: defaultdict(list))
        for i in operonDict.keys():
            gene1 = "646797728 YP_003639887 TherJR_1122 cytochrome C family protein [Thermincola sp. JR: NC_014152]"
            gene2 = "646799199 YP_003641333 TherJR_2595 hypothetical protein [Thermincola sp. JR: NC_014152]"
            gene3 = "646796949 YP_003639120 TherJR_0333 hypothetical protein [Thermincola sp. JR: NC_014152]"
            if gene1 in operonDict[i]["headers"] and gene2 in operonDict[i]["headers"] and gene3 in operonDict[i][
                "headers"]:
                genome = blastresult.split("-thermincola.blas")[0]
                category = "iron_reduction"
                for j in range(0, len(operonDict[i]["orfs"])):
                    orf = (operonDict[i]["orfs"][j])
                    header = (operonDict[i]["headers"][j])
                    evalue = blastDict[orf]["e"]
                    out.write(category + "," + genome + "," + orf + "," + replace(header, [","], ";") + "," + "evalue: " + str(evalue) + "," + str(counter) + "\n")
                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                counter += 1
out.close()

os.system("cat %s/FinalSummary-dereplicated-clustered.csv %s/GeoThermin.csv > %s/FinalSummary-dereplicated-clustered-blast.csv" % (outDirectory, outDirectory, outDirectory))


# ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************
# print("Filtering out likely false positives")
clusterDict = defaultdict(lambda: defaultdict(list))
summary = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory, "r")
for i in summary:
    if not re.match(r'#', i):
        ls = i.rstrip().split(",")
        clusterDict[ls[5]]["line"].append(ls)
        clusterDict[ls[5]]["gene"].append(ls[3])
        clusterDict[ls[5]]["category"].append(ls[0])


out = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % outDirectory, "w")
for i in (clusterDict.keys()):
    ls = (clusterDict[i]["gene"])
    if "EetA.hmm" in ls or "EetB.hmm" in ls or "Ndh2.hmm" in ls or "FmnB.hmm" in ls or "FmnA.hmm" in ls or "DmkA.hmm" in ls or "DmkB.hmm" in ls or "PplA.hmm" in ls:
        fleet = ["EetA.hmm", "EetB.hmm", "Ndh2.hmm", "FmnB.hmm", "FmnA.hmm", "DmkA.hmm", "DmkB.hmm", "PplA.hmm"]

        if unique(ls, fleet) < 5:  # If there are less than 5 FLEET genes in the cluster
            if len(remove2(ls, fleet)) < 1:  # If FLEET genes are the only ones in the cluster
                pass
            else:  # If there are other genes in the cluster that are not FLEET
                for j in clusterDict[i]["line"]:
                    if j[3] not in fleet:  # avoiding the fleet genes
                        out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        else:  # if there are 5 or more of the FLEET genes present within cluster
            for j in clusterDict[i]["line"]:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    elif "MamB.hmm" in ls or "MamC.hmm" in ls or "MamD.hmm" in ls or "MamF.hmm" in ls or "MamG.hmm" in ls or "MamJ.hmm" \
            in ls or "MamK.hmm" in ls or "MamM.hmm" in ls or "MamP.hmm" in ls or "MamW.hmm" in ls or "MamX.hmm" in ls or "MamY.hmm" in ls:
        mam = ["MamB.hmm", "MamC.hmm", "MamD.hmm", "MamF.hmm", "MamG.hmm", "MamJ.hmm", "MamK.hmm", "MamM.hmm", "MamP.hmm", "MamW.hmm", "MamX.hmm", "MamY.hmm"]

        if unique(ls, mam) < 7:
            if len(remove2(ls, mam)) < 1:
                pass
            else:
                for j in clusterDict[i]["line"]:
                    if j[3] not in mam:
                        out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        else:
            for j in clusterDict[i]["line"]:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    elif "MtoA.hmm" in ls or "MtrA.hmm" in ls or "MtrC_TIGR03507.hmm" in ls or "MtrB_TIGR03509.hmm" in ls:
        if "MtoA.hmm" in ls and "MtrB_TIGR03509.hmm" in ls:
            for j in clusterDict[i]["line"]:
                if j[3] in ["MtrB_TIGR03509.hmm", "MtoA.hmm"]:
                    out.write("iron_oxidation" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                else:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        if "MtrA.hmm" in ls and "MtrB_TIGR03509.hmm" in ls:
            for j in clusterDict[i]["line"]:
                if j[3] in ["MtrA.hmm", "MtrB_TIGR03509.hmm"]:
                    out.write("iron_reduction" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                else:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "MtrB_TIGR03509.hmm" not in ls:
            pass

    elif "FoxA.hmm" in ls or "FoxB.hmm" in ls or "FoxC.hmm" in ls:
        foxabc = ["FoxA.hmm", "FoxB.hmm", "FoxC.hmm"]
        if unique(ls, foxabc) < 2:
            pass
        else:
            for j in clusterDict[i]["line"]:
                if j[3] not in foxabc:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    elif "FoxE.hmm" in ls or "FoxY.hmm" in ls or "FoxZ.hmm" in ls:
        foxeyz = ["FoxE.hmm", "FoxY.hmm", "FoxZ.hmm"]
        if unique(ls, foxeyz) < 2:
            pass
        else:
            for j in clusterDict[i]["line"]:
                if j[3] not in foxeyz:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    elif "Cyc1.hmm" in ls:
        if "Cyc2_repCluster3.hmm" not in ls and "Cyc2_repCluster2.hmm" not in ls and "Cyc2_repCluster1.hmm" not in ls:
            pass

    elif "CymA.hmm" in ls:
        if "MtrB_TIGR03509.hmm" not in ls and "MtrA.hmm" not in ls and "MtoA.hmm" not in ls and "MtrC_TIGR03507.hmm" not in ls:
            pass

    elif "iron_aquisition-siderophore_synthesis" in clusterDict[i]["category"] or \
                    "iron_aquisition-siderophore_transport" in clusterDict[i]["category"] or \
                    "iron_aquisition-iron_transport" in clusterDict[i]["category"] or "iron_aquisition-heme_transport" in clusterDict[i]["category"]:

        if len(ls) > 1:
            for j in clusterDict[i]["line"]:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        else:
            if "FutA1-iron_ABC_transporter_iron-binding-rep.hmm" in ls or "FutA2-iron_ABC_transporter_iron-binding-rep.hmm" in ls \
                    or "FutC-iron_ABC_transporter_ATPase-rep.hmm" in ls or "LbtU-LvtA-PiuA-PirA-RhtA.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls \
                    or "LbtU_LbtB-legiobactin_receptor_2.hmm" in ls or "IroC-salmochelin_transport-rep.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                pass

    else:
        linels = (clusterDict[i]["line"])
        for j in linels:
            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

        out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

out.close()


# ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
# print("Writing heatmap-compatible CSV")
cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport", "iron_aquisition-heme_oxygenase", "iron_aquisition-siderophore_synthesis",
        "iron_aquisition-siderophore_transport", "iron_gene_regulation", "iron_oxidation", "iron_reduction",
        "iron_storage", "magnetosome_formation"]

Dict = defaultdict(lambda: defaultdict(list))
final = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % args.out, "r")
for i in final:
    ls = (i.rstrip().split(","))
    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
        if not re.match(r'#', i):
            process = ls[0]
            cell = ls[1]
            orf = ls[2]
            gene = ls[3]
            Dict[cell][process].append(gene)

normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in os.listdir(args.bin_dir):
    if lastItem(i.split(".")) == args.bin_ext:
        file = open("%s/%s-proteins.faa" % (args.bin_dir, i), "r")
        file = fasta(file)
        normDict[i] = len(file.keys())


outHeat = open("%s/FeGenie-heatmap-data.csv" % args.out, "w")
outHeat.write("X" + ',')
for i in sorted(Dict.keys()):
    outHeat.write(i + ",")
outHeat.write("\n")

for i in cats:
    outHeat.write(i + ",")
    for j in sorted(Dict.keys()):
        if not re.match(r'#', j):
            outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(args.inflation)) + ",")
    outHeat.write("\n")

outHeat.close()


# REMOVING FILES
os.system("rm %s/GeoThermin.csv" % args.out)
os.system("rm %s/*summary*" % args.out)
os.system("rm %s/FinalSummary-dereplicated-clustered-blast.csv" % args.out)
os.system("rm %s/*blast" % args.out)
os.system("rm %s/FinalSummary.csv" % args.out)
os.system("rm %s/FinalSummary-dereplicated-clustered.csv" % args.out)
os.system("mv %s/FinalSummary-dereplicated-clustered-blast-filtered.csv %s/FeGenie-summary.csv" % (args.out, args.out))


# CROSS-VALIDATION AGAINST REFERENCE DATABASE
if args.ref != "NA":
    print("")
    print("Performing Diamond BLASTx search of putative iron genes against reference database")
    summary = open("%s/FeGenie-summary.csv" % args.out, "r")
    out = open("%s/FeGenie-summary.fasta" % args.out, "w")
    for i in summary:
        if not re.match(r'#', i):
            ls = (i.rstrip().split(","))
            seq = (BinDict[ls[1]][ls[2]])
            header = (">" + ls[1] + "|" + ls[2])
            out.write(header + "\n")
            out.write(seq + "\n")
    out.close()

    os.system("diamond blastp --db %s.dmnd --out "
              "%s/FeGenie-summary.dmndout --max-target-seqs 1 --outfmt 6 "
              "--threads %s --query %s/FeGenie-summary.fasta --quiet" % (args.ref, args.out, str(args.t), args.out))

    dmndblast = open("%s/FeGenie-summary.dmndout" % args.out)
    dmndblastDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'No_Hits')))
    for hit in dmndblast:
        lst = hit.rstrip().split("\t")
        evalue = lst[10]
        cell = lst[0].split("|")[0]
        orf = lst[0].split("|")[1]
        target = refDict[lst[1]]
        dmndblastDict[cell][orf]["e"] = evalue
        dmndblastDict[cell][orf]["target"] = target


summary = open("%s/FeGenie-summary.csv" % args.out, "r")
out = open("%s/FeGenie-summary-blasthits.csv" % args.out, "w")
if args.ref != "NA":
    out.write("category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore_ratio" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
else:
    out.write("category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore_ratio" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "protein_sequence" + "\n")


counter = 1
for i in summary:
    if not re.match(r'#', i):
        ls = i.rstrip().split(",")
        seq = BinDict[ls[1]][ls[2]]
        hemes = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                + len(re.findall(r'C(...............)CH', seq))
        if args.ref != "NA":
            blasthit = dmndblastDict[ls[1]][ls[2]]["target"]
            e = dmndblastDict[ls[1]][ls[2]]["e"]
            try:
                out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(float(metaDict[ls[3].split(".")[0]]) / float(ls[4])) + "," + str(counter) + "," + str(hemes) + "," + blasthit + "," + str(e) + "," + seq + "\n")
            except TypeError:
                out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(hemes) + "," + blasthit + "," + str(e) + "," + seq + "\n")

        else:
            try:
                out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(float(metaDict[ls[3].split(".")[0]]) / float(ls[4])) + "," + str(counter) + "," + str(hemes) + "," + seq + "\n")

            except TypeError:
                out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(hemes) + "," + seq + "\n")

    else:
        counter += 1
        out.write(i)

out.close()


# REMOVING FILES
if args.ref != "NA":
    os.system("rm %s/FeGenie-summary.dmndout" % args.out)
    os.system("rm %s/FeGenie-summary.fasta" % args.out)

os.system("rm %s/FeGenie-summary.csv" % args.out)
os.system("mv %s/FeGenie-summary-blasthits.csv %s/FeGenie-summary.csv" % (args.out, args.out))


# ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
if args.makeplots == "y":
    if conda == 0:
        Rdir = args.R
    else:
        os.system("echo ${rscripts} > r.txt")
        file = open("r.txt")
        for i in file:
            Rdir = (i.rstrip())
        os.system("rm r.txt")

    os.system("Rscript -e 'install.packages(\"ggplot2\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"reshape\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"reshape2\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"tidyverse\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"argparse\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"ggdendro\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"ggpubr\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"grid\", repos = \"http://cran.us.r-project.org\")\'")

    os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-heatmap-data.csv %s" % (Rdir, args.out, args.out))
    os.system("Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s" % (Rdir, args.out, args.out))
    print("\n\n\n")
    print("...")

    # ******** CHECKING ON SUCCESS OF RSCRIPTS ***********
    DIR = os.listdir(args.out)
    count = 0
    for i in DIR:
        if not re.match(r'\.', i):
            count += 1

    if count == 2:
        print("Looks like Rscript has not performed succesfully. This, unfortunately, is a very finicky part of the pipeline. "
          "The CSV files have, nonetheless, been successfully created, so you can take that data and plot if manually as you wish. "
          "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

    if count > 2 and count < 5:
        print("Looks like at least one plot was generated by Rscript, but there was likely an error with one of the scripts. "
              "The main CSV output should be present, however, so that you can plot the data as you wish on your own. "
              "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

    if count == 5:
        print("Looks like Rscript ran succesfully! Congrats on this. Hopefully, the resulting plots are of use to you.")


print("")
print("Pipeline finished without crashing!!! Thanks for using :)")



