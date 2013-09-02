#!/usr/bin/env python
import os
import argparse
import logging
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a fastq file and creates a table containing time and x y coordinate information.")
    parser.add_argument("-i", "--input", dest="input", action='store', required=True, help="Name of FASTQ file [Required]")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Name of the output CSV table [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def headerSplit(header):
    if header.count(' '):
        # Header was created using CASAVA 1.8+ 
        (mhead, suphead) = header.split(' ')
        (machine, runid, flowid, lane, tile, xcoord, ycoord) = mhead.split(':')
        read = suphead[0]

    else:
        # Header was created using older versions of CASAVA
        (mhead, suphead) = header.split('#')
        (machine, lane, tile, xcoord, ycoord) = mhead.split(':')
        if suphead.endswith('1') or suphead.endswith('2'):
            read = suphead[-1]
        else:
            read = 1
    plane = tile[0]
    swath = tile[1]
    tileNum = tile[2:]
    return([machine, lane, plane, swath, tileNum, xcoord, ycoord, read])

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    with open(args.input, 'r') as FQ:
        with open(args.out, 'w') as OUT:
            OUT.write(','.join(str(x) for x in ['machine', 'lane', 'plane', 'swath', 'tileNum', 'x-coord', 'y-coord', 'readNum', 'sequence']) + "\n")
            for (header, sequence, quality) in FastqGeneralIterator(FQ):
                headerInfo = headerSplit(header)
                OUT.write(','.join(str(x) for x in headerInfo) + ',' + sequence + "\n")

if __name__=='__main__':
    main()
    logging.info("Script complete.")
