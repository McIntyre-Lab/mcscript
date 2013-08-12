#!/usr/bin/env python
import cProfile
import argparse
import logging
from Bio import SeqIO
import collections
import re
import numpy
import os

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool for identifying duplicates and creating various useful output')
    parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A list of fq file [Required]')
    parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
    parser.add_argument('--pixel', dest='pix', action='store', default=100, required=False, help='Number of pixels to consider a read as an optical duplicate [Default:100]')
    parser.add_argument('-a', dest='append', action='store_true', help='This flag will cause the output dataset to be appended too.')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """Function for handling error logging"""
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')


def readFASTQ(fname,count):
    """Read a fastq file and store information in a dictionary where the key is
       the sequence and the value is a list of headers. Also returns a total
       number of reads."""

    logging.info("Reading the FASTQ file.")
    myDict = collections.defaultdict(list)
    recTuple = collections.namedtuple('recTuple', 'name lane tile x y read')
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            match = parseHeader(record.name)
            if len(match) == 5:
                mytup = recTuple(record.name,*match)
            elif len(match) == 4:
                mytup = recTuple(record.name,*match,read=1)
            else:
                raise ValueError('Check regex, parsing header is broken')
            myDict[str(record.seq)].append(mytup)
            count += 1
    logging.info("Finished reading the FASTQ file.")
    return(myDict,count)

def parseHeader(fqName):
    """Function to parse the FASTQ header line
        Where [0] = Lane
              [1] = Tile
              [2] = X-coord
              [3] = Y-coord
              [4] = Read number """
    match = re.search('.*:?([0-9]):([0-9]+):([0-9]+):([0-9]+).*\/?([12])*$',fqName)
    matchList = filter(None,match.groups())     # Removes any empty strings, ie if there is no PE information
    intList = [int(x) for x in matchList]
    return(intList)

def checkOptical(myList, pix):
    """Given a list of FASTQ headers, this function will check if reads are
       within a certain 'args.pix' distance of each other. If there are within
       this distance, then they will be called optical duplicates."""
    # Create a dictionary to store optical duplicate information.
    myDict = dict()
    # Create a list of sets, where each set contains headers that are within
    # args.pix of each other.
    listOfSets = identifySets(myList, pix)
    # reduce the set list so that each set is a group of optical duplicates.
    redSetList = reduceSet(listOfSets)
    return(redSetList)

def identifySets(myList, pix):
    """This function steps through a list of headers and creates sets of those
       that fall within the args.pix range. These need to be reduced."""
    # Create list to store results
    setList = list()
    # Compare each header and create sets of headers that overlap.
    for item1 in myList:
        # Parse the first header to grab coordinate information.
        item1Set = {item1}
        coord1 = numpy.array([item1.x,item1.y])
        for item2 in myList:
            coord2 = numpy.array([item2.x,item2.y])
            # Don't compare a header to itself.
            if item1.name != item2.name: 
                # Parse the second header to grab coordinate information.
                # Test if headers are on the same lane and tile. If they are
                # then determine if they are within args.pix of each other
                # using numpy's norm function.
                if item1.lane == item2.lane and item1.tile == item2.tile:
                    eucDist = numpy.linalg.norm(coord1-coord2)
                    if eucDist <= pix:
                        item1Set.add(item2)
        setList.append(item1Set)
    return(setList)

def reduceSet(setList):
    """Function to step through a list of sets and combine sets that overlap.
       This will create a unique list of sets, where each set contains the headers
       of reads that are within args.pix of each other"""
    setList2 = [setList[0]]
    for item1 in setList[1:]:
        inSetList2 = False
        for item2 in setList2:
            if item1 & item2:
                item2 |= item1
                inSetList2 = True
                break
        if not inSetList2:
            setList2.append(item1)
    return(setList2)
    
def dupCounts(setList,cnt1,cnt2,cnt3):
    """This function calculates various counts and returns a list of counts."""
    pdup_cnt = len(setList)
    opdup_cnt = 0
    opdupset_cnt = 0
    for item in setList:
        if len(item) > 1:
            opdup_cnt += len(item)
            opdupset_cnt += 1
    cnt1 = pdup_cnt - opdupset_cnt
    cnt2 = opdup_cnt
    cnt3 = opdupset_cnt
    return(cnt1,cnt2,cnt3)

def writeOutput(handle, myList):
    """Function to write Output"""
    handle.write(','.join(str(x) for x in myList) + "\n")

def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    # Initialize main counters
    total_read = 0     # number of reads
    uniq_read = 0      # number of unique reads
    uniq_seq = 0       # number of unique sequences
    pdupCnt = 0        # number of PCR duplicates
    opdupCnt = 0       # number of optical duplicates
    opdupSetCnt = 0    # number of sets of optical duplicates

    # Read in the FASTQ file and return a dictionary with SEQ as the key and
    # list(HEADERS) as the values. Also return the total number of reads.
    logging.info("Starting to read in FASTQ file and creating dictionary.")
    myFASTQ, total_read = readFASTQ(args.fq,total_read)
    logging.info("Finished reading FASTQ file.")

    # Calculate the number of unique SEQUENCES. 
    uniq_seq = len(myFASTQ)

    # Loop through each sequence in the dictionary and examine it for
    # duplicates and optical duplicates.
    logging.info("Starting to identify duplicates.")
    for key in myFASTQ:
        # Copy the dictionary value using list() so that I don't modify the
        # dictionary value later in the script.
        myValue = list(myFASTQ[str(key)])   
        myLength = len(myValue)

        if myLength == 1:
            # If there is only one header for a given sequence, then this is a unique read.
            uniq_read += 1
        elif myLength > 1:
            # If there are more than one header for a given, then these are reads.
            total_dup = myLength
            # Identify if the duplicates are optical dups or something else.
            readSet = checkOptical(myValue, args.pix)
            (pdupCnt, opdupCnt, opdupSetCnt) = dupCounts(readSet, pdupCnt, opdupCnt, opdupSetCnt)
        else:
            logging.error("There is something wrong with the fastq dictionary at key %s" % key)
    logging.info("Finished identifying duplicates.")

    per_uniq_read = round(float(uniq_read)/total_read * 100, 2)
    myOutHeader= ["File_Name","Total_Reads", "Num_Unique_Reads", "Per_Uniq_Reads", "Num_Unique_Seq", "Num_PCR_Dups", "Num_OP_Dups", "Num_OP_Dup_Sets"]
    myOut = [args.fq,total_read, uniq_read, per_uniq_read, uniq_seq, pdupCnt, opdupCnt, opdupSetCnt]

    # Test if we want to append the file, if so handle the headers correctly.
    logging.info("Writing output.")
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
                logging.error("Could not open output file, it must be busy")
    else:
        try:
            with open(args.out, 'w') as handle:
                writeOutput(handle,myOutHeader)
                writeOutput(handle,myOut)
        except:
            logging.error("Could not open output file, it must be busy")
    logging.info("Script Finished.")

if __name__ == '__main__':
    main()
