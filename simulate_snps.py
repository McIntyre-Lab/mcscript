#!/usr/bin/env python
import argparse
import logging
from Bio import SeqIO
import numpy as np
import random

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool to simulate SNPs into a FASTA Reference')
    parser.add_argument('-i','--input',dest='fname', action='store', required=True, help='FASTA file [Required]')
    parser.add_argument('-o','--out', dest='oname', action='store', required=True, help='FASTA file with SNPs incorporated [Required]')
    parser.add_argument('-n','--num', dest='nsnp', action='store', required=True, help='Number of SNPs to incorporate incorporated [Required]')
    parser.add_argument('-t','--table', dest='tname', action='store', required=False, help='Table to store SNP locations [Optional]')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
    args = parser.parse_args()
    #args = parser.parse_args(['-i', '/home/jfear/tmp/fa/fb551_si_fusions.fa', '-o', '/home/jfear/tmp/mut_test.fa', '-n', '160000', '-t', '/home/jfear/tmp/mut_test.csv'])
    return(args)

def generateChromTable(fname):
    """ Indexes a FASTA file for quick access to each chromosome """
    seqDict = SeqIO.index(fname, 'fasta')
    return(seqDict)
        
def selectBase(base):
    """ Randomly select a different base than the oldBase """
    while True:
        newBase = random.choice('ATCG')
        if base != newBase:
            break
    return(newBase)

def mutChrom(seq, count, table, chrom):
    """ Take a sequence and mutate it 'count' number of times. Record
    information in a table """
    # Convert the Bio.Seq.seq object into a mutable form
    seqMut = seq.tomutable()

    # Generate 'count' number of random mutations in the sequence
    randMut = np.random.randint(0, len(seq), count)
    for mut in randMut:
        oldBase = seqMut[mut]
        newBase = selectBase(oldBase)
        seqMut[mut] = newBase
        pos = mut + 1 # coordinate convert for our output table
        table.append([chrom, pos , oldBase, newBase])
    return(seqMut, table)

def writeTable(tname, table):
    """ Write table output """
    header = ['chrom', 'pos', 'ref', 'alt']
    table.sort()
    with open(tname, 'w') as TOUT:
        TOUT.write(','.join(str(x) for x in header) + "\n")
        for row in table:
            TOUT.write(','.join(str(x) for x in row) + "\n")

def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    table = list()

    logging.info("Indexing FASTA file and creating dictionary")
    seqDict = generateChromTable(args.fname)

    # Randomly select chromosomes to mutate
    logging.info("Generating a random list of %s mutation" % args.nsnp)
    randChrom = np.random.randint(0, len(seqDict), int(args.nsnp))

    # count the number of times each chromosome was in our list or mutations
    logging.info("Summarizing mutation count for each chromosome/fusion")
    binCount = np.bincount(randChrom)

    logging.info("Mutating Sequences")
    with open(args.oname, 'w') as OUT:
        # Loop through chromosomes and simulate mutation
        for num,chrom in enumerate(seqDict):
            try:
                count = binCount[num]
            except:
                count = 0

            curr = seqDict[chrom]
            if count > 0:
                curr.seq, table = mutChrom(curr.seq, count, table, chrom)
            SeqIO.write(curr, OUT, 'fasta')

    if args.tname:
        logging.info("Writing table output")
        writeTable(args.tname, table)

if __name__ == '__main__':
    main()
    logging.info("Script complete")
