#!/usr/bin/env python

import logging
import argparse 
from Bio import SeqIO

def getOptions():
    """ Function to pull in arguments """

    parser = argparse.ArgumentParser(description="Converts a FASTQ file into a FASTA file.")
    parser.add_argument("-i", "--input", dest="fname", action='store', required=True, help="Name of input FASTQ file [Required]")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Name of output FASTA file [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
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

    logging.info("Converting '%s' to '%s'" % (args.fname,args.out))
    SeqIO.convert(args.fname, "fastq", args.out, "fasta")
    logging.info("Finished converting '%s' to '%s'" % (args.fname,args.out))


if __name__=='__main__':
    main()
    logging.info("Script complete.")
