#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
from collections import defaultdict

# Add-on packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# McLab packages
import mclib

def getOptions():
    """ Function to pull in arguments """
    description = """ This script takes an input fasta file of fusions and identifies all of the distinct fusions."""
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    group1 = parser.add_argument_group(description="Input Files")
    group1.add_argument('-f','--fasta', dest='fname', action='store', required=True, help='A fasta file [Required]')

    group2 = parser.add_argument_group(description="Output Files")
    group2.add_argument("-o", dest="oname", action='store', required=True, help="Name of output FASTA. [Required]")
    group2.add_argument("--log", dest="log", action='store', required=False, help="Name of the LOG file [Optional]") 

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 

    args = parser.parse_args()
    return(args)

def readFA(db):
    """ Iterates over a FASTA file and builds a dictionary where the keys are
    sequences and the values is a list of IDs. Since a dictionary key can only
    bd used once, it will create a distinct list of sequences
    
    Arguments:
    ----------
    db (obj) = a SeqIO.index object, which is a dictionary of SeqRecords

    
    Returns:
    --------
    seqDict (dict) = Dictionary where the keys are sequences and the values is a list of fusion_id
    
    """

    # Use defaultdict to make a dictionary whose values are a list
    seqDict = defaultdict(list)

    # Iterate through FASTA file and build a dictionary of distinct sequences
    for record in db:
        #TODO: replace 'pass' with code to build dictionary where 'sequence' is
        #      key and append 'fusion_id' to value
        seqDict[str(record.seq).append(record.fusion_id)]
        


    #TODO: Add logging information to count how many fusions we start with and
    #      how many sequences are in the dictionary
    logging.info("Reading the FASTA file, creating dictionary.")
    return(seqDict)

def writeFA(oname,mydict,db):
    """ Iterate over a list of distinct sequences. For each sequence grab the
    first fusion ID. Pull this fusion out of the database and write it a FASTA
    file. 

    Arguments:
    ----------
    oname (str) = output file
    mydict (dict) = Dictionary where the keys are sequences and the values is a list of fusion_id
    db (obj) = a SeqIO.index object, which is a dictionary of SeqRecords

    Returns:
    --------
    Writes FASTA file.
    """

    # open the output file
    with open(oname, 'w') as OUT:
        # Iterate through distinct sequences and grab a representative fusion_id
        for record in mydict:
            #TODO: record is going to be a sequence, so you want to pull out
            #      the first fusion_id
            print record.fusion_id[0]

            # Pull the SeqRecord out of database
            #TODO: To pull record, hint use fusion_id as the key for the db dictionary
            seqrecord = db['']

            # Write Sequence to FASTA file
            #TODO: figure out how to write
            #      http://biopython.org/wiki/SeqIO
            SeqIO.write(record,OUT,"fasta")

def main(args):
    """ Main Script """
    # Create FASTA Random Access Database
    logger.info("Reading FASTA file")
    db = SeqIO.index(args.fname, 'fasta')

    # Create a dictionary of distinct sequences
    #TODO: add logging info
    distinct = readFA(db)

    # Write a list of distinct sequences out to a fasta
    #TODO: add logging info
    writeFA(args.oname, distinct, db)


if __name__ == '__main__':
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
