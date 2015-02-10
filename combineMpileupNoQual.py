#!/usr/bin/env python

import logging
import argparse 
from collections import defaultdict

def getOptions():
    """ Function to pull in arguments """

    parser = argparse.ArgumentParser(description="This script takes mpileup files from replicates and combine them into a single mpileup where the count is the sum of the counts across replicates.")
    parser.add_argument("-i", "--input", dest="fname", nargs='+',action='store', required=True, help="List of MPILEUP files to combine [Required]")
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Output MPILEUP file [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    #args = parser.parse_args(['-i','/home/jfear/tmp/mpileup/r101_M1.mpileup', '/home/jfear/tmp/mpileup/r101_M2.mpileup', '/home/jfear/tmp/mpileup/r101_M3.mpileup', '-o', '/home/jfear/tmp/mpileup/r101_combine.mpileup'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    # I am using a dictionary to store all of the mpileup data, this will allow
    # me to easily add positions together. This has the potential to be memory
    # intensive but it should be pretty fast and easy to code.
    mdict = defaultdict(dict)
    for fname in args.fname:
        # Parse each pileup file and add its values to the storage dictionary mdict
        with open(fname, 'r') as PILE:
            for row in PILE:
                cols = row.rstrip().split('\t')
                count = int(cols[3])
                if count == 0:
                    continue
                try:
                    curr = mdict[cols[0]][cols[1]] 
                    curr[1] = curr[1] + count
                    curr[2][0] = curr[2][0]
                    curr[2][1] = curr[2][1]
                except:
                    mdict[cols[0]][cols[1]] = [cols[2], int(cols[3])]

    with open(args.oname, 'w') as OUT:
        # Write out mdict to a new mpileup file
        for chrom in mdict:
            for pos in mdict[chrom]:
                base = mdict[chrom][pos][0]
                count = mdict[chrom][pos][1]
                myout = [chrom, pos, base, count] 
                OUT.write('\t'.join(str(x) for x in myout) + "\n")

if __name__=='__main__':
    main()
    logging.info("Script complete.")

