#!/usr/bin/env python
import os, sys
import re
import argparse 
import logging
import collections
from Bio import SeqIO

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description='Tool to identify homopolymers in FASTA and FASTQ sequence.')
    parser.add_argument("-i","--input", dest="input", action='store', required=True, help="Input file in either FASTA or FASTQ format [Required]")
    parser.add_argument("--fastq", dest="ft", action='store_const', const='fq', required=False, help="Specify FASTQ format") 
    parser.add_argument("--fasta", dest="ft", action='store_const', const='fa', required=False, help="Specify FASTA format") 
    parser.add_argument('-a', dest='append', action='store_true', help='This flag will cause the output dataset to be appended too.')
    parser.add_argument("--plot", dest="plt", action='store', required=False, help="Name of the plot (png) if you would like plot output.") 
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Out File", metavar="OUT_FILE")
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')

def selectFormat(fname):
    filetype = ''
    fileName, fileExtension = os.path.splitext(fname)
    if fileExtension == '.fq' or fileExtension == '.fastq':
        filetype = 'fq'
    elif fileExtension == '.fa' or fileExtension == '.fasta':
        filetype = 'fa'
    else:
        raise ValueError('Please specify if the input file is "fastq" or "fasta"')
    return filetype

def writeOutput(handle, myList):
    """Function to write Output"""
    handle.write(','.join(str(x) for x in myList) + "\n")

def plotDist(dataDict):
    X = list()
    Y = list()
    for key in sorted(dataDict, key=lambda key: len(key)):
        X.append(len(key))
        Y.append(dataDict[key])
    return X,Y


def main():
    args = getOptions()                         # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)
    if not args.ft:                             # Figure out what kind of file we have if the user did not tell us.
        args.ft = selectFormat(args.input)

    logging.info("Begining to analyze homopolymers.")
    with open(args.input, "r") as FH:
        if args.ft == 'fa':
            seqIter = SeqIO.parse(FH, 'fasta')
        elif args.ft == 'fq':
            seqIter = SeqIO.parse(FH, 'fastq')

        # initialize lists()
        Ares, Tres, Cres, Gres = ([] for i in range(4))
        # initialize counts
        total_read = 0
        homopoly = 0
        for record in seqIter:
            if args.ft == 'fq': total_read += 1
            A = [x.group() for x in re.finditer('[Aa]{6,}',str(record.seq))]
            T = [x.group() for x in re.finditer('[Tt]{6,}',str(record.seq))]
            C = [x.group() for x in re.finditer('[Cc]{6,}',str(record.seq))]
            G = [x.group() for x in re.finditer('[Gg]{6,}',str(record.seq))]

            if A or T or C or G:
                if args.ft == 'fq': homopoly += 1
                Ares.extend(A)
                Tres.extend(T)
                Cres.extend(C)
                Gres.extend(G)

        Acnt = collections.Counter(Ares)
        Tcnt = collections.Counter(Tres)
        Ccnt = collections.Counter(Cres)
        Gcnt = collections.Counter(Gres)

    if args.ft == 'fq':
        myOutHeader = ["File_Name", "Num_reads", "Num_reads_w_homopolymers", "per_reads_w_homopolymers", "Num_homopolymer_A", "Num_homopolymer_T", "Num_homopolymer_C", "Num_homopolymer_G"]
        per_homopoly = float(homopoly) / float(total_read) * 100
        myOut = [args.input, total_read, homopoly, per_homopoly, sum(Acnt.values()), sum(Tcnt.values()), sum(Ccnt.values()), sum(Gcnt.values())]
    else:
        myOutHeader = ["File_Name", "Num_homopolymer_A", "Num_homopolymer_T", "Num_homopolymer_C", "Num_homopolymer_G"]
        myOut = [args.input, sum(Acnt.values()), sum(Tcnt.values()), sum(Ccnt.values()), sum(Gcnt.values())]

    logging.info("Wrigting Output")

    print args.out
    if args.append:
        if os.path.exists(args.out):
            try:
                with open(args.out, 'a') as handle:
                    writeOutput(handle,myOut)
            except:
                logging.error("Could not open output file, it must be busy")
        else:
            try:
                with open(args.out, 'a') as handle:
                    writeOutput(handle,myOutHeader)
                    writeOutput(handle,myOut)
            except:
                logging.error("Could not open output file")
    else:
        try:
            with open(args.out, 'w') as handle:
                writeOutput(handle,myOutHeader)
                writeOutput(handle,myOut)
        except:
            logging.error("Could not open output file")

    if args.plt:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        logging.info("Creating Plots")
        X,Y = plotDist(Acnt)
        plt.plot(X,Y, label="AAAAAA+")
        X,Y = plotDist(Tcnt)
        plt.plot(X,Y, label="TTTTTT+")
        X,Y = plotDist(Ccnt)
        plt.plot(X,Y, label="CCCCCC+")
        X,Y = plotDist(Gcnt)
        plt.plot(X,Y, label="GGGGGG+")
        plt.legend(loc='upper right')
        plt.xlabel("Size of Homopolymer Repeat", fontsize=18)
        plt.ylabel("Number of Reads", fontsize=18)
        plt.savefig(args.plt)

    logging.info("Script Finished.")

if __name__=='__main__':
    main()
