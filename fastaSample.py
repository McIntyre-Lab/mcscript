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

def setLogger(args):
    """ Function to set up log handling 
    Arguments
    ========
    args = object with all command line arguments
    
    Returns: a logger that will output general info/debug/warnings to STDOUT
    and all errors to STDERR
    """
    # If debug flag is given output additional debug info
    if args.debug:
        _lvl = logging.DEBUG
    else:
        _lvl = logging.INFO
    logger = logging.getLogger()
    logger.setLevel(_lvl)

    ## Set logging level for the different output handlers. 
    ## ERRORs to STDERR, and everything else to STDOUT
    loghandler = logging.StreamHandler(stream=sys.stdout)
    errhandler = logging.StreamHandler(stream=sys.stderr)
    loghandler.setLevel(_lvl)
    errhandler.setLevel(logging.ERROR)

    # Format the log handlers
    logfmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    loghandler.setFormatter(logfmt)
    errhandler.setFormatter(logfmt)

    # Add handler to the main logger
    logger.addHandler(loghandler)
    logger.addHandler(errhandler)

    return logger

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

    logger.info('Sampling records from {} and writing to {}'.format(args.fname, args.oname))
    with open(args.oname, 'w') as OUT:
        with open(args.fname, 'r') as FA:
            cnt = 0
            for record in SeqIO.parse(FA, 'fasta'):
                if cnt in sampleNum:
                    SeqIO.write(record, OUT, 'fasta')
                cnt += 1

if __name__=='__main__':
    # Pull in command line options
    args = getOptions()

    # Turn on logging
    logger = setLogger(args)

    # Run Main part of the script
    ## Count the number of records in the file
    numRecords = countRecords(args.fname)

    ## Randomly sample records
    sampleRecrods(args, numRecords)

    logger.info("Script complete.")
