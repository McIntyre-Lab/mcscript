#!/usr/bin/env python

import sys, csv, argparse, logging, pprint, re
from Bio import SeqIO
logger = logging.getLogger()

def getOptions():
    description="""This script converts a FASTA file to a 3-column bed file."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-f", "--fasta", dest="fasta", action='store', required=True, help="fasta file [Required]",metavar="FASTA_FILE")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="output file name [Required]", metavar="OUT_FILE")
    parser.add_argument("-g", "--log", dest="log", action='store', help="Log File", metavar="LOG_FILE")
    args = parser.parse_args()
    return(args)

def setLogger(fname):
    global logger
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    flog = logging.FileHandler(fname)
    logger.setLevel(logging.INFO)
    flog.setLevel(logging.INFO)
    flog.setFormatter(formatter)
    logger.addHandler(flog)

def main():
    args = getOptions()
    if args.log:        # If logging is turned on define log file
        setLogger(args.log)

    with open(args.out, 'w') as OUT:
        writer = csv.writer(OUT,delimiter='\t')
        logger.info("Reading the '%s' file" % args.fasta)

        with open(args.fasta, 'r') as FA:
            for record in SeqIO.parse(FA, "fasta"):    # Read each FASTA block and create bed output
                writer.writerow([record.id,'0',len(record.seq)])

if __name__=='__main__':
    main()
