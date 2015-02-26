#!/usr/bin/env python

# Built-in packages
import argparse 
from argparse import RawDescriptionHelpFormatter
import os.path
import sys
import logging

# Imort 3rd party packages
import numpy as np
from numpy.random import choice
from Bio import SeqIO

# McLab Packages
import mclib

def getOptions():
    """ Function to pull in arguments """
    description="""This script takes a FASTA file and samples records randomly with replacement then outputs a new FASTA file.

    NOTE: This script requires numpy and Biopython.
    """
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-f", dest="fname", action='store', required=True, help="Name of FASTA file to be sampled from [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of output FASTA file [Required]")
    parser.add_argument("-n", dest="nrec", action='store', type=int, required=True, help="Number of FASTA records to sample without replacement [Required]")
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    args = parser.parse_args()
    #args = parser.parse_args(['-f', '/home/jfear/tmp/dmel-all-chromosome-r5.51.fasta', '-o', '/home/jfear/tmp/output.fa','-n', '2', '--debug'])
    return(args)

def countRecords(fname):
    """ Count the number of reqcords in the original FASTA file. 
    Arguments
    =========
    fname = file name of a FASTA file.
    
    Returns: the number of records in a FASTA file
    """

    logger.info('Counting the number of records in FASTA file: {}'.format(fname))
    count = 0
    with open(fname, 'r') as FA:
        for row in FA:
            if row.startswith('>'):
                count += 1
    logger.info('FASTA file contained {} records'.format(count))

    return count

def sampleRecrods(args, count):
    """ Randomly sample without replacement numpy records.
    Arguments
    ========
    args = object with all command line arguments
    count = the number of records in the FASTA file.
    
    Writes output to file
    """

    # Use numpy to randomly choose record numbers
    if count >= args.nrec:
        sampleNum = choice(count, args.nrec, replace=False)
        logger.debug('Sampled record indices: {}'.format(sampleNum))
    else:
        logger.error('There are only {} records in your original file and you are tyring to sample {}.'.format(count, args.nrec))
        sys.exit(1)

    logger.info('Sampling {} records from {} and writing to {}'.format(args.nrec, args.fname, args.oname))
    with open(args.oname, 'w') as OUT:
        with open(args.fname, 'r') as FA:
            cnt = 0
            for record in SeqIO.parse(FA, 'fasta'):
                if cnt in sampleNum:
                    SeqIO.write(record, OUT, 'fasta')
                cnt += 1

def main(args):
    ## Count the number of records in the file
    numRecords = countRecords(args.fname)

    ## Randomly sample records
    sampleRecrods(args, numRecords)

if __name__=='__main__':
    # Turn on Logging if option -g was given
    args = getOptions()

    # Turn on logging
    logger = logging.getLogger()
    if args.debug:
        mclib.logger.setLogger(logger, args.log, 'debug')
    else:
        mclib.logger.setLogger(logger, args.log)

    # Output git commit version to log, if user has access
    mclib.git.git_to_log(__file__)

    # Run Main part of the script
    main(args)
    logger.info("Script complete.")
