#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def main():
    def Strip(ls):
        outList = []
        for i in ls:
            gene = i.split("|")[0]
            outList.append(gene)
        return outList


    def unique(ls, ls2):
        unqlist = []
        for i in ls:
            if i not in unqlist and i in ls2:
                unqlist.append(i)
        return len(unqlist)


    def Unique(ls):
        unqList = []
        for i in ls:
            if i not in unqList:
                unqList.append(i)
        return unqList


    def Unique2(ls):
        unqList = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in unqList:
                unqList.append(hmm)
        return unqList


    def checkFe(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_reduction", "iron_oxidation"]:
                    count += 1
        return count


    def checkDFE1(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if hmm in ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"]:
                    count += 1
        return count


    def checkDFE2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if hmm in ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"]:
                    count += 1
        return count


    def check1(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_transport", "iron_aquisition-heme_transport"]:
                    count += 1
        return count


    def check1_2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count


    def check2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"]:
                    count += 1
        return count


    def check3(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count


    def checkReg(ls):
        count = 0
        for i in ls:
            hmm = i.split("|")[0]
            if re.findall(r'aquisition', geneToCatDict[hmm]):
                count += 1
        return count


    def checkMam(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] == "magnetosome_formation":
                    count += 1
        return count


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
        for i in range(0, length - 1):
            x += iterable.split(delim)[i]
            x += delim
        return x[0:len(x) - 1]


    def secondToLastItem(ls):
        x = ''
        for i in ls[0:len(ls) - 1]:
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


    def fastaRename(fasta_file):
        counter = 0
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[header] = seq
                    header = i[1:]
                    header = header.split(" ")[0]
                    counter += 1
                    header = header + "_" + str(counter)
                    seq = ''
                else:
                    header = i[1:]
                    header = header.split(" ")[0]
                    counter += 1
                    header = header + "_" + str(counter)
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
    
            )`-.--.  )\.---.     )\.-.    )\.---.   )\  )\  .'(   )\.---.  
            ) ,-._( (   ,-._(  ,' ,-,_)  (   ,-._( (  \, /  \  ) (   ,-._( 
            \ `-._   \  '-,   (  .   __   \  '-,    ) \ (   ) (   \  '-,   
             ) ,_(    ) ,-`    ) '._\ _)   ) ,-`   ( ( \ \  \  )   ) ,-`   
            (  \     (  ``-.  (  ,   (    (  ``-.   `.)/  )  ) \  (  ``-.  
             ).'      )..-.(   )/'._.'     )..-.(      '.(    )/   )..-.(                                                                                    
                                  %(?/////////&//%                                                
              .,,.                   (%((&@@@#/*.                      .,,.        
              .,,.                     @(((/&@@@#///**                  ...        
                                         #&((///////////////*/@                                
                                                             #*@.                             
                                      ()                   * )//*
                                      <^^>             *     (/*   .
                                     .-""-.                  *)
                          .---.    ."-....-"-._     _...---''`/. '
                         ( (`\ \ .'            ``-''    _.-"'`
                          \ \ \ : :.                 .-'
                           `\`.\: `:.             _.'
                           (  .'`.`            _.'
                            ``    `-..______.-'
                                      ):.  (
                                    ."-....-".
                                  .':.        `.
                                  "-..______..-"
    
        Image design: Nancy Merino (2018);
        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=5)", default=5)

    parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format",
                        default="NA")

    parser.add_argument('-out', type=str, help="name output directory (default=fegenie_out)",
                        default="fegenie_out")

    parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                        default=1000)

    parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--orfs', type=str, help="include this flag if you are providing bins as open-reading frames or genes in FASTA amino-acid format", const=True,
                        nargs="?")

    parser.add_argument('--meta', type=str, help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies", const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each iron gene category to be normalized to "
                             "the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, FeGenie will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, FeGenie will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--makeplots', type=str,
                        help="include this flag if you would like FeGenie to make some figures from your data?. "
                             "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Please be sure to install all the required R packages as instrcuted in the FeGenie Wiki: "
                             "https://github.com/Arkadiy-Garber/FeGenie/wiki/Installation. "
                             "If you see error or warning messages associated with Rscript, you can still expect to "
                             "see the main output (CSV files) from FeGenie.", const=True, nargs="?")

    # CHECKING FOR CONDA INSTALL
    os.system("echo ${iron_hmms}/HMM-bitcutoffs.txt > HMMlib.txt")
    file = open("HMMlib.txt")
    for i in file:
        location = i.rstrip()

    try:
        bits = open(location)
        conda = 1
    except FileNotFoundError:
        conda = 0

    if conda == 0:
        parser.add_argument('-hmm_lib', type=str,
                            help='HMM database; directory titled \'HMM-lib\', can be found in the FeGenie folder',
                            default="NA")

        parser.add_argument('-R', type=str,
                            help="location of R scripts directory (note: this optional argument requires Rscript to be "
                                 "installed on your system). The R scripts directory is in the same directory as the "
                                 "FeGenie python code", default="NA")

    args = parser.parse_known_args()[0]

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

    if args.makeplots:
        if conda == 0:
            if args.R != "NA":
                print(".")
            else:
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
        # raise SystemExit
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
            if not args.gbk:

                if args.orfs:
                    testFile = open("%s/%s" % (binDir, i), "r")
                    print("ORFS for %s found. Skipping Prodigal, and going with %s" % (i, i))
                    for line in testFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \|")
                                print(
                                    "Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                else:
                    try:
                        testFile = open("%s/%s-proteins.faa" % (binDir, i), "r")
                        print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                        for line in testFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print("Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                    except FileNotFoundError:
                        binFile = open("%s/%s" % (binDir, i), "r")
                        for line in binFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print("Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                        print("Finding ORFs for " + cell)
                        if args.meta:
                            os.system("prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -p meta -q" % (
                            binDir, i, binDir, i, binDir, i))
                        else:
                            os.system(
                                "prodigal -i %s/%s -a %s/%s-proteins.faa -o %s/%s-prodigal.out -q" % (
                                binDir, i, binDir, i, binDir, i))
            else:
                os.system('gtt-genbank-to-AA-seqs -i %s/%s -o %s/%s.faa' % (binDir, i, binDir, i))

                faa = open("%s/%s.faa" % (binDir, i))
                faa = fasta(faa)

                gbkDict = defaultdict(list)
                counter = 0

                count = 0
                gbk = open("%s/%s" % (binDir, i))
                for gbkline in gbk:
                    ls = gbkline.rstrip()
                    if re.findall(r'/locus_tag', ls):
                        count += 1

                if count > 0:
                    print(i)
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])
                            counter += 1

                        if counter > 0:
                            gbkDict[locus].append(locusTag)
                            counter = 0
                else:
                    # print(i)
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)
                            counter += 1

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])

                        if counter > 0:
                            gbkDict[locus].append(altContigName)
                            counter = 0

                idxOut = open("%s/%s-proteins.idx" % (binDir, i), "w")
                faaOut = open("%s/%s-proteins.faa" % (binDir, i), "w")

                for gbkkey1 in gbkDict.keys():
                    # print(i)
                    counter = 0
                    for gbkey2 in gbkDict[gbkkey1]:
                        print(gbkey2)
                        counter += 1
                        if len(faa[gbkey2]) > 0:
                            newOrf = gbkkey1 + "_" + str(counter)
                            # print(j + "\t\t" + newOrf + "\t\t" + str(faa[j]))
                            idxOut.write(gbkey2 + "," + newOrf + "\n")
                            faaOut.write(">" + newOrf + "\n")
                            faaOut.write(str(faa[gbkey2]) + "\n")
                    # print("")

                idxOut.close()
                faaOut.close()

                # orfs = open("%s/%s" % (binDir, i))
                # orfs = fastaRename(orfs)
                # out = open("%s/%s-proteins.faa" % (binDir, i), "w")
                # for key in orfs.keys():
                #     out.write(">" + key + "\n")
                #     out.write(orfs[key] + "\n")

            if args.orfs:
                file = open("%s/%s" % (binDir, i))
            else:
                file = open("%s/%s-proteins.faa" % (binDir, i))
            file = fasta(file)
            for j in file.keys():
                orf = j.split(" # ")[0]
                BinDict[cell][orf] = file[j]

    # ******************** READ BITSCORE CUT-OFFS INTO HASH MEMORY ****************************** #
    if conda == 1:
        os.system("echo ${iron_hmms} > HMMlib.txt")
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
        if not re.match(r'\.', FeCategory) and FeCategory not in ["HMM-bitcutoffs.txt", "FeGenie-map.txt"]:
            print("")
            print(".")
            print("Looking for following iron-related functional category: " + FeCategory)
            hmmDir = "%s/%s/" % (HMMdir, FeCategory)
            hmmDirLS2 = os.listdir("%s/%s" % (HMMdir, FeCategory))

            HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
            for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
                # if lastItem(i.split(".")) == args.bin_ext and not re.findall(r'proteins.faa', i):  # FILTERING OUT ANY NON-BIN-RELATED FILES
                if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                    # print("analyzing: " + i)
                    # FASTA = open(binDir + i + "-proteins.faa", "r")
                    # FASTA = fasta(FASTA)
                    os.system(
                        "mkdir " + binDir + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

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
                            if args.orfs:
                                os.system(
                                    "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/%s"
                                    % (int(args.t), float(bit), binDir, i, hmm, binDir, i, hmm, hmmDir, hmm, binDir, i)
                                )
                            else:
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

                    print("")
                    os.system("rm -rf " + binDir + "/" + i + "-HMM")

            out = open(args.out + "/%s-summary.csv" % (FeCategory), "w")
            out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "\n")
            for key in HMMdict.keys():
                for j in HMMdict[key]:
                    out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                              str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "\n")

            out.close()

    print("\n")
    print("Consolidating summary files into one master summary file")
    out = open(args.out + "/FinalSummary.csv", "w")
    if args.ref == "NA":
        out.write("category" + "," + "cell" + "," + "orf" + "," + "related_hmm" + "," + "HMM-bitscore" + "\n")

    resultsDir = os.listdir(args.out)
    for i in resultsDir:
        if lastItem(i.split("-")) == "summary.csv":
            result = open(args.out + "/" + i, "r")
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
    summary = open(args.out + "/FinalSummary.csv", "r")
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
                print(j)
                CoordDict[i][contig].append(int(numOrf))

    counter = 0
    print("Clustering ORFs...")
    print("")
    out = open(args.out + "/FinalSummary-dereplicated-clustered.csv", "w")
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

                        out.write(
                            SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"] + "," +
                            str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "," + ncbiHomolog + "," + str(
                                SummaryDict[i][orf]["NCBIeval"]) + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1

                    else:
                        out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf][
                            "hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1

                else:
                    if args.ref != "NA":
                        for l in RemoveDuplicates(k):
                            orf = j + "_" + str(l)
                            ncbiHomolog = SummaryDict[i][orf]["NCBImatch"]
                            if ncbiHomolog != "NA":
                                ncbiHomolog = ncbiHomolog.split("]")[0] + "]"

                            out.write(
                                SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"]
                                + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(
                                    counter) + "," + ncbiHomolog + "," + str(SummaryDict[i][orf]["NCBIeval"]) + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1
                    else:
                        for l in RemoveDuplicates(k):
                            orf = j + "_" + str(l)

                            out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf][
                                "hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1
    out.close()

    # ************************** BLAST-BASED METHODS/LOOKING FOR UNMODELED MARKERS ********************************
    thermincola = "%s/iron_reduction/non-aligned/TherJR_SLCs.faa" % HMMdir
    geobacter = "%s/iron_reduction/non-aligned/geobacter_PCCs.faa" % HMMdir

    print("Looking for Thermincola S-layer cytochromes and Geobacter porin-cytochromes")

    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext:
            if args.orfs:
                os.system(
                    "makeblastdb -dbtype prot -in %s/%s -out %s/%s -logfile %s/makedbfile.txt" % (
                        binDir, i, binDir, i, binDir))
                os.system("rm %s/makedbfile.txt" % binDir)

                os.system(
                    "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-thermincola.blast -evalue 1E-10"
                    % (thermincola, binDir, i, args.t, args.out, i))

                os.system(
                    "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-geobacter.blast -evalue 1E-10"
                    % (geobacter, binDir, i, args.t, args.out, i))

            else:
                os.system(
                    "makeblastdb -dbtype prot -in %s/%s -out %s/%s -logfile %s/makedbfile.txt" % (
                    binDir, i, binDir, i, binDir))
                os.system("rm %s/makedbfile.txt" % binDir)

                os.system(
                    "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-thermincola.blast -evalue 1E-10"
                    % (thermincola, binDir, i, args.t, args.out, i))

                os.system(
                    "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-geobacter.blast -evalue 1E-10"
                    % (geobacter, binDir, i, args.t, args.out, i))

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

    out = open(args.out + "/GeoThermin.csv", "w")
    for blastresult in os.listdir(args.out):
        blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        if re.findall(r'geobacter\.blast', blastresult):
            blast = open(args.out + "/" + blastresult, "r")

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
                        out.write(category + "," + genome + "," + orf + "," + replace(header, [","],
                                                                                      ";") + "," + "evalue: " + str(
                            evalue) + "," + str(counter) + "\n")
                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1
                    # print("\n\n")

        if re.findall(r'thermincola\.blast', blastresult):
            blast = open(args.out + "/" + blastresult, "r")

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
                        out.write(category + "," + genome + "," + orf + "," + replace(header, [","],
                                                                                      ";") + "," + "evalue: " + str(
                            evalue) + "," + str(counter) + "\n")
                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    counter += 1
    out.close()

    summary = open("%s/FinalSummary-dereplicated-clustered.csv" % args.out)
    out = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % args.out, "w")
    DeRepDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in summary:
        ls = i.rstrip().split(",")
        unq = ls[1] + "|" + ls[2]
        DeRepDict[unq] = ls[0]
        out.write(i.rstrip() + "\n")

    blastHits = open("%s/GeoThermin.csv" % args.out)
    for i in blastHits:
        ls = i.rstrip().split(",")
        unq = ls[1] + "|" + ls[2]
        if unq not in DeRepDict.keys():
            out.write(i.rstrip() + "\n")

    out.close()

    # ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************
    # print("Filtering out likely false positives")
    fleet = ["EetA.hmm", "EetB.hmm", "Ndh2.hmm", "FmnB.hmm", "FmnA.hmm", "DmkA.hmm", "DmkB.hmm", "PplA.hmm"]
    mam = ["MamA.hmm", "MamB.hmm", "MamE.hmm", "MamK.hmm", "MamP.hmm", "MamM.hmm", "MamP.hmm", "MamQ.hmm", "MamI.hmm",
           "MamL.hmm", "MamO.hmm"]
    foxabc = ["FoxA.hmm", "FoxB.hmm", "FoxC.hmm"]
    foxeyz = ["FoxE.hmm", "FoxY.hmm", "FoxZ.hmm"]

    clusterDict = defaultdict(lambda: defaultdict(list))
    summary = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % args.out, "r")
    for i in summary:
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            clusterDict[ls[5]]["line"].append(ls)
            clusterDict[ls[5]]["gene"].append(ls[3])
            clusterDict[ls[5]]["category"].append(ls[0])

    out = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % args.out, "w")
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

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:  # if there are 5 or more of the FLEET genes present within cluster
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "MamA.hmm" in ls or "MamB.hmm" in ls or "MamE.hmm" in ls or "MamK.hmm" in ls or "MamM.hmm" in ls or "MamO.hmm" \
                in ls or "MamP.hmm" in ls or "MamQ.hmm" in ls or "MamI.hmm" in ls or "MamL.hmm" in ls:
            mam = ["MamA.hmm", "MamB.hmm", "MamE.hmm", "MamK.hmm", "MamP.hmm", "MamM.hmm", "MamP.hmm", "MamQ.hmm",
                   "MamI.hmm", "MamL.hmm", "MamO.hmm"]
            if unique(ls, mam) < 5:
                if len(remove2(ls, mam)) < 1:
                    pass
                else:
                    for j in clusterDict[i]["line"]:
                        if j[3] not in mam:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "MtoA.hmm" in ls or "MtrA.hmm" in ls or "MtrC_TIGR03507.hmm" in ls or "MtrB_TIGR03509.hmm" in ls:
            if "MtoA.hmm" in ls and "MtrB_TIGR03509.hmm" in ls and "MtrC_TIGR03507.hmm" not in ls:
                for j in clusterDict[i]["line"]:
                    if j[3] in ["MtrB_TIGR03509.hmm", "MtoA.hmm", "CymA.hmm"]:
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
                if len(remove2(ls, foxabc)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[3] not in foxabc:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
            else:

                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "FoxE.hmm" in ls or "FoxY.hmm" in ls or "FoxZ.hmm" in ls:
            foxeyz = ["FoxE.hmm", "FoxY.hmm", "FoxZ.hmm"]
            if "FoxE.hmm" not in ls:
                if len(remove2(ls, foxeyz)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[3] not in foxeyz:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")
                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
            else:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "DFE_0448.hmm" in ls or "DFE_0449.hmm" in ls or "DFE_0450.hmm" in ls or "DFE_0451.hmm" in ls:
            DFE1 = ["DFE_0448.hmm", "DFE_0449.hmm", "DFE_0450.hmm", "DFE_0451.hmm"]

            if unique(ls, DFE1) < 3:
                if len(remove2(ls, DFE1)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[3] not in DFE1:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
            else:

                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "DFE_0461.hmm" in ls or "DFE_0462.hmm" in ls or "DFE_0463.hmm" in ls or "DFE_0464.hmm" in ls or "DFE_0465.hmm" in ls:
            DFE2 = ["DFE_0461.hmm", "DFE_0462.hmm", "DFE_0463.hmm", "DFE_0464.hmm", "DFE_0465"]

            if unique(ls, DFE2) < 3:
                if len(remove2(ls, DFE2)) < 1:
                    pass

                else:
                    for j in clusterDict[i]["line"]:
                        if j[3] not in DFE2:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
            else:

                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "Cyc1.hmm" in ls:
            if "Cyc2_repCluster3.hmm" not in ls and "Cyc2_repCluster2.hmm" not in ls and "Cyc2_repCluster1.hmm" not in ls:
                pass

            else:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

        elif "CymA.hmm" in ls:
            if "MtrB_TIGR03509.hmm" not in ls and "MtrA.hmm" not in ls and "MtoA.hmm" not in ls and "MtrC_TIGR03507.hmm" not in ls:
                pass

        elif "iron_aquisition-siderophore_synthesis" in clusterDict[i]["category"] or \
                        "iron_aquisition-siderophore_transport" in clusterDict[i]["category"] or \
                        "iron_aquisition-iron_transport" in clusterDict[i][
                    "category"] or "iron_aquisition-heme_transport" in clusterDict[i]["category"]:

            if len(Unique(ls)) > 1:
                for j in clusterDict[i]["line"]:
                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            else:
                if "FutA1-iron_ABC_transporter_iron-binding-rep.hmm" in ls or "FutA2-iron_ABC_transporter_iron-binding-rep.hmm" in ls \
                        or "FutC-iron_ABC_transporter_ATPase-rep.hmm" in ls or "LbtU-LvtA-PiuA-PirA-RhtA.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls \
                        or "LbtU_LbtB-legiobactin_receptor_2.hmm" in ls or "IroC-salmochelin_transport-rep.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls:
                    for j in clusterDict[i]["line"]:
                        out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write(
                        "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                else:
                    pass

        else:
            linels = (clusterDict[i]["line"])
            for j in linels:
                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

    out.close()

    # ****************************** REMOVING SINGLETONS *************************************
    '''
    cluDict = defaultdict(list)
    cluDictIndex = defaultdict(list)
    summary = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % outDirectory, "r")
    for i in summary:
        if not re.match(r'#', i):
            ls = (i.rstrip().split(","))
            cluDict[ls[5]].append(i.rstrip())
            cluDictIndex[ls[5]].append(ls)
    
    acceptableSingltons = ["sulfocyanin.hmm", "Cyc2_repCluster1.hmm", "Cyc2_repCluster2.hmm", "Cyc2_repCluster3.hmm",
                           "Transferrin_TbpB_binding_protein_Haemophilus_influenzae_P44971.hmm",
                           "PF13668_Ferritin_like_domain.hmm", "PF00210-Ferritin_like_domain.hmm",
                           "Sid_YqjI_regulator_for_YqjH_P64588_Escherichia_coli_180606.hmm",
                           "Sid_PvdS_regulator_Paeruginosa_PA2426_180620.hmm",
                           "Sid_PchR_pyochelin_regulator_Pseudomonas_aeruginosa_PA4227_180623.hmm",
                           "Sid_FpvI_regulator_PA2387_Paeruginosa_PAO1_180620.hmm",
                           "PF09012_sub-FeoC_like_transcriptional_regulator.hmm", "PF04773_FecR.hmm",
                           "PF02742-Iron_dependent_repressor-dtxR_family_metal_binding_domain.hmm",
                           "PF01475-Iron_dependent_repressor-fur_family.hmm",
                           "PF01325-Iron_dependent_repressor-dtxR_family_N.hmm"]
    
    out = open("%s/FinalSummary-dereplicated-clustered-blast-filtered2.csv" % outDirectory, "w")
    for i in cluDict.keys():
        if (len(cluDict[i])) == 1:
            ls = cluDictIndex[i][0]
            if ls[3] in acceptableSingltons:
                for j in cluDict[i]:
                    out.write(j + "\n")
                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        else:
            for j in cluDict[i]:
                out.write(j + "\n")
            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
    
    out.close()
    '''

    # REMOVING FILES
    os.system("rm %s/GeoThermin.csv" % args.out)
    os.system("rm %s/*summary*" % args.out)
    os.system("rm %s/FinalSummary-dereplicated-clustered-blast.csv" % args.out)
    os.system("rm %s/*blast" % args.out)
    os.system("rm %s/FinalSummary.csv" % args.out)
    os.system("rm %s/FinalSummary-dereplicated-clustered.csv" % args.out)
    # os.system("rm %s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % args.out)
    os.system("mv %s/FinalSummary-dereplicated-clustered-blast-filtered.csv %s/FeGenie-summary.csv" % (args.out, args.out))

    # OPTIONAL CROSS-VALIDATION AGAINST REFERENCE DATABASE
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
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
    else:
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "protein_sequence" + "\n")

    # SUMMARIZING CROSS-REFERNCE RESULTS AND COUNTING HEME-BINDING MOTIFS
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
                    out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                        metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(hemes) + "," + blasthit + "," + str(
                        e) + "," + seq + "\n")
                except TypeError:
                    out.write(
                        ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(
                            hemes) + "," + blasthit + "," + str(e) + "," + seq + "\n")

            else:
                try:
                    out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                        metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(hemes) + "," + seq + "\n")

                except TypeError:
                    out.write(
                        ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(
                            hemes) + "," + seq + "\n")

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

    # FILTERING OUT FALSE POSITIVES FOR SIDEROPHORE GENES
    MAP = open(HMMdir + "/FeGenie-map.txt", "r")
    mapDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in MAP:
        ls = i.rstrip().split("\t")
        mapDict[ls[0]] = ls[1]

    out = open(args.out + "/FeGenie-summary-fixed.csv", "w")
    geneToCatDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    memoryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
    clusterDict = defaultdict(list)
    infile = open(args.out + "/FeGenie-summary.csv")
    for i in infile:
        if not re.match(r'#', i):
            ls = i.rstrip().split(",")
            if ls[6] != "cluster":
                if not re.findall(r'defaultdict', ls[5]):
                    clu = ls[6]
                    cat = ls[0]
                    dataset = ls[1]
                    orf = ls[2]
                    hmm = allButTheLast(ls[3], ".")
                    clusterDict[clu].append(hmm + "|" + dataset + "|" + orf)
                    geneToCatDict[hmm] = cat
                    hmm = allButTheLast(ls[3], ".")
                    memoryDict[dataset][orf]["cat"] = ls[0]
                    memoryDict[dataset][orf]["gene"] = ls[3]
                    memoryDict[dataset][orf]["bit"] = ls[4]
                    memoryDict[dataset][orf]["cutoff"] = ls[5]
                    memoryDict[dataset][orf]["clu"] = clu
                    memoryDict[dataset][orf]["heme"] = ls[7]
                    memoryDict[dataset][orf]["seq"] = ls[8]
                else:
                    cat = ls[0]
                    dataset = ls[1]
                    orf = ls[2]
                    clu = ls[7]
                    hmm = ls[3]
                    memoryDict[dataset][orf]["cat"] = ls[0]
                    memoryDict[dataset][orf]["gene"] = ls[3]
                    memoryDict[dataset][orf]["bit"] = ls[4]
                    memoryDict[dataset][orf]["cutoff"] = "evalue-cutoff: 1E-10"
                    memoryDict[dataset][orf]["clu"] = ls[7]
                    memoryDict[dataset][orf]["heme"] = ls[8]
                    memoryDict[dataset][orf]["seq"] = ls[9]
                    geneToCatDict[hmm] = cat
                    clusterDict[ls[7]].append(hmm + "|" + dataset + "|" + orf)
            else:
                out.write(i.rstrip())

    for i in clusterDict.keys():
        out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        for j in clusterDict[i]:
            hmm = j.split("|")[0]
            dataset = j.split("|")[1]
            orf = j.split("|")[2]
            cat = memoryDict[dataset][orf]["cat"]

            if cat in ["iron_aquisition-siderophore_transport", "iron_aquisition-heme_transport"]:
                if len(Unique2(clusterDict[i])) < 2:
                    break
                elif check1(clusterDict[i]) < 2:
                    pass
                else:

                    out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                              memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                              memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                              memoryDict[dataset][orf]["seq"] + "\n")

            elif cat in ["iron_aquisition-siderophore_synthesis"]:
                if len(Unique2(clusterDict[i])) < 3:
                    break
                elif check1_2(clusterDict[i]) < 3:
                    pass
                else:

                    out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                              memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                              memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                              memoryDict[dataset][orf]["seq"] + "\n")

            elif cat in ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"]:
                if len(Unique2(clusterDict[i])) < 2:
                    break
                elif check2(clusterDict[i]) < 2:
                    pass
                else:
                    out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                              memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                              memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                              memoryDict[dataset][orf]["seq"] + "\n")

            elif cat == "iron_gene_regulation":
                if checkReg(clusterDict[i]) < 1:
                    pass
                else:
                    out.write(
                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                        memoryDict[dataset][orf]["bit"] + "," +
                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                        memoryDict[dataset][orf]["heme"] + "," +
                        memoryDict[dataset][orf]["seq"] + "\n")

            elif cat == "magnetosome_formation":
                if checkMam(clusterDict[i]) < 5:
                    pass
                else:
                    out.write(
                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                        memoryDict[dataset][orf]["bit"] + "," +
                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                        memoryDict[dataset][orf]["heme"] + "," +
                        memoryDict[dataset][orf]["seq"] + "\n")

            elif cat == "iron_oxidation":
                if hmm == "Cyc1":
                    if checkFe(clusterDict[i]) < 2:
                        pass
                    else:
                        out.write(
                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")

                elif hmm in ["MtoA", "MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507"]:
                    operon = clusterDict[i]
                    operon = Strip(operon)
                    # print("iron oxidation")
                    # print(operon)
                    # print(hmm)
                    # print("")
                    if "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon and "MtrA" in operon:
                        out.write(
                            "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction")

                    elif "MtoA" in operon and "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction")

                    elif "MtrB_TIGR03509" in operon and "MtrA" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction_or_oxidation")

                    elif "MtoA" in operon and "MtrB_TIGR03509" in operon:
                        out.write(
                            "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction_or_oxidation")

                    elif "MtrC_TIGR03507" in operon and "MtrB_TIGR03509" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction_or_oxidation")

                    else:
                        pass

                else:
                    out.write(
                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                        memoryDict[dataset][orf]["bit"] + "," +
                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                        memoryDict[dataset][orf]["heme"] + "," +
                        memoryDict[dataset][orf]["seq"] + "\n")

            elif cat == "iron_reduction":
                if hmm in ["DFE_0465", "DFE_0464", "DFE_0463", "DFE_0462", "DFE_0461"]:
                    if checkDFE1(clusterDict[i]) < 3:
                        pass
                    else:
                        out.write(
                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")

                elif hmm in ["DFE_0451", "DFE_0450", "DFE_0449", "DFE_0448"]:
                    if checkDFE2(clusterDict[i]) < 3:
                        pass
                    else:
                        out.write(
                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")

                elif hmm in ["MtrC_TIGR03507", "MtrA", "MtrB_TIGR03509", "MtoA"]:
                    operon = clusterDict[i]
                    # print("iron reduction")
                    # print(operon)
                    operon = Strip(operon)
                    # print(operon)
                    # print(hmm)
                    # print("")
                    if "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon and "MtrA" in operon:
                        out.write(
                            "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction")

                    elif "MtoA" in operon and "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction")

                    elif "MtrB_TIGR03509" in operon and "MtrA" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction_or_oxidation")

                    elif "MtoA" in operon and "MtrB_TIGR03509" in operon:
                        out.write(
                            "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction")

                    elif "MtrC_TIGR03507" in operon and "MtrB_TIGR03509" in operon:
                        out.write(
                            "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                        # print("writing iron_reduction_or_oxidation")

                    else:
                        pass

                else:
                    out.write(
                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                        memoryDict[dataset][orf]["bit"] + "," +
                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                        memoryDict[dataset][orf]["heme"] + "," +
                        memoryDict[dataset][orf]["seq"] + "\n")

            else:
                out.write(
                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                    memoryDict[dataset][orf]["bit"] + "," +
                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                    memoryDict[dataset][orf]["heme"] + "," +
                    memoryDict[dataset][orf]["seq"] + "\n")

    out.close()
    os.system("mv %s/FeGenie-summary-fixed.csv %s/FeGenie-summary.csv" % (args.out, args.out))

    # ****************************** PRE-FINAL ALTERATION OF THE OUTPUT FILE ***************************************
    clu = 0
    summaryDict = defaultdict(list)
    summary = open(args.out + "/FeGenie-summary.csv")

    out = open(args.out + "/FeGenie-summary-altered.csv", "w")
    if args.ref != "NA":
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
    else:
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_binding_motifs" + "," + "protein_sequence" + "\n")

    for i in summary:
        if re.search(r'#', i):
            clu += 1
        else:
            summaryDict[clu].append(i.rstrip())

    if args.gbk:
        idxDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for idxfile in binDirLS:
            if lastItem(idxfile.split(".")) == "idx":
                print(idxfile)
                idxfileopen = open("%s/%s" % (binDir, idxfile))
                for idxline in idxfileopen:
                    ls = idxline.rstrip().split(",")
                    newOrf = ls[1]
                    oldOrf = ls[0]
                    idxDict[newOrf] = oldOrf

        for i in summaryDict.keys():
            if len(summaryDict[i]) > 0:
                for j in summaryDict[i]:
                    ls = j.split(",")
                    if args.ref != "NA":
                        out.write(
                            ls[0] + "," + ls[1] + "," + str(idxDict[ls[2]]) + "," + ls[3] + "," + ls[4] + "," + ls[
                                5] + "," + ls[6] + "," + ls[7] + "," + ls[8] + "," + ls[9] + "," + ls[10] + "\n")
                    else:
                        out.write(
                            ls[0] + "," + ls[1] + "," + str(idxDict[ls[2]]) + "," + ls[3] + "," + ls[4] + "," + ls[
                                5] + "," + ls[6] + "," + ls[7] + "," + ls[8] + "\n")
                out.write(
                    "#####################################################################################################"
                    "#####################################################################################################\n")
    else:
        for i in summaryDict.keys():
            if len(summaryDict[i]) > 0:
                for j in summaryDict[i]:
                    out.write(j + "\n")
                out.write(
                    "#####################################################################################################"
                    "#####################################################################################################\n")
    out.close()

    os.system("mv %s/FeGenie-summary-altered.csv %s/FeGenie-geneSummary-clusters.csv" % (args.out, args.out))
    os.system("rm %s/FeGenie-summary.csv" % args.out)

    # ****************************** REMOVING #'S ***************************************
    summary = open("%s/FeGenie-geneSummary-clusters.csv" % args.out, "r")
    out = open("%s/FeGenie-geneSummary.csv" % args.out, "w")
    if args.ref != "NA":
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
    else:
        out.write(
            "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_binding_motifs" + "," + "protein_sequence" + "\n")

    for i in summary:
        if not re.search(r'#', i):
            out.write(i.rstrip() + "\n")

    out.close()

    # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
    cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport", "iron_aquisition-heme_oxygenase",
            "iron_aquisition-siderophore_synthesis",
            "iron_aquisition-siderophore_transport", "iron_gene_regulation", "iron_oxidation",
            "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
            "iron_reduction", "iron_storage", "magnetosome_formation"]

    Dict = defaultdict(lambda: defaultdict(list))
    final = open("%s/FeGenie-geneSummary-clusters.csv" % args.out, "r")
    for i in final:
        ls = (i.rstrip().split(","))
        if not re.search(r'#', i):
            if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                if not re.match(r'#', i):
                    process = ls[0]
                    cell = ls[1]
                    orf = ls[2]
                    gene = ls[3]
                    Dict[cell][process].append(gene)

    normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in os.listdir(args.bin_dir):
        if lastItem(i.split(".")) == args.bin_ext and not re.findall(r'-proteins.faa', i):
            if args.orfs:
                file = open("%s/%s" % (args.bin_dir, i), "r")
            else:
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
                if args.norm:
                    outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(args.inflation)) + ",")
                else:
                    outHeat.write(str(len(Dict[j][i])) + ",")
        outHeat.write("\n")

    outHeat.close()

    # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
    if args.makeplots:
        if conda == 0:
            Rdir = args.R
        else:
            os.system("echo ${rscripts} > r.txt")
            file = open("r.txt")
            for i in file:
                Rdir = (i.rstrip())
            os.system("rm r.txt")

        if args.norm:
            os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-heatmap-data.csv %s" % (Rdir, args.out, args.out))
            os.system("Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s" % (Rdir, args.out, args.out))
        else:
            os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/FeGenie-heatmap-data.csv %s" % (Rdir, args.out, args.out))
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
            print(
                "Looks like Rscript has not performed succesfully. This, unfortunately, is a very finicky part of the pipeline. "
                "The CSV files have, nonetheless, been successfully created, so you can take that data and plot if manually as you wish. "
                "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

        if count > 2 and count < 5:
            print(
                "Looks like at least one plot was generated by Rscript, but there was likely an error with one of the scripts. "
                "The main CSV output should be present, however, so that you can plot the data as you wish on your own. "
                "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

        if count == 5:
            print("Looks like Rscript ran succesfully! Congrats on this. Hopefully, the resulting plots are of use to you.")

    print("")
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()





