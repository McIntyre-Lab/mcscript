#!/usr/bin/env python
import os
import logging
import argparse 
from Bio import SeqIO

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Parses a BLAT table and creates Flags as well as a summary table.")
    parser.add_argument("-b", "--blat", dest="bname", action='store', required=True, help="Name of the BLAT file [Required]")
    parser.add_argument("-f", "--fastq", dest="fname", action='store', required=True, help="Name of the FASTQ file that BLAT used [Required]")
    parser.add_argument("--outdir", dest="out", action='store', required=True, help="Name of the directory to create output tables [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def pslParse(bname):
    with open(bname,'r') as PSL:
        for row in PSL:
            # skip all of the header lines
            if row.startswith('-'):
                break
        idList = list()
        for row in PSL:
            idList.append(row.split('\t')[9])

    idSet = set(idList)
    return(idSet)

def fastqParse(fname):
    indexName = os.path.splitext(fname)[0] + '.idx'
    dbObj = SeqIO.index_db(indexName, fname, 'fastq')
    return(dbObj)


def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    bname = '/home/jfear/tmp/psl/r1.psl'
    fname = '/home/jfear/tmp/fq/r1.fq'
    oname = '/home/jfear/tmp/bob.csv'

    pslSet = pslParse(bname)
    fqDB = fastqParse(fname)

    with open(oname,'w') as OUT:
        for read in fqDB.itervalues():
            if read.id in pslSet:
                OUT.write(read.id +',1\n')
            else:
                OUT.write(read.id +',0\n')
            




if __name__=='__main__':
    main()
    logging.info("Script complete.")
