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

def readFASTQ(fname):
    """Read a fastq file and store information in a dictionary where the key is
       the sequence and the value is a list of headers. Also returns a total
       number of reads."""

    logging.info("Reading the FASTQ file.")
    mySeqDict = collections.defaultdict(list)
    myFqDict = dict()
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            mySeqDict[str(record.seq)].append(record.name)
            match = parseHeader(record.name)
            if len(match) == 4: # if there is no read information, append a 1
                match.append(1)
            myReadDict[record.name] = {'lane':match[0], 'tile':match[1],'x':match[2],'y':match[3],'read':match[4]}
    logging.info("Finished reading the FASTQ file.")
    return(mySeqDict,myReadDict)

def parseHeader(fqName):
    """Function to parse the FASTQ header line
                [0] = Lane
                [1] = tile
                [2] = x-coord
                [3] = y-coord
                [4] = read number"""
    match = re.search('.*:?([0-9]):([0-9]+):([0-9]+):([0-9]+).*\/?([1-2])*',fqName)
    matchList = filter(None,match.groups())     # Removes any empty strings, ie if there is no PE information
    matchIntList = [int(x) for x in matchList]
    return(matchIntList)

def checkOptical(myList, readDict, pix):
    """Given a list of FASTQ headers, this function will check if reads are
       within a certain 'args.pix' distance of each other. If there are within
       this distance, then they will be called optical duplicates."""
    # Create a dictionary to store optical duplicate information.
    myDict = dict()
    # Create a list of sets, where each set contains headers that are within
    # args.pix of each other.
    listOfSets = identifySets(myList, readDict, pix)
    # reduce the set list so that each set is a group of optical duplicates.
    redSetList = reduceSet(listOfSets)
    return(redSetList)

def identifySets(myList, readDict, pix):
    """This function steps through a list of headers and creates sets of those
       that fall within the args.pix range. These need to be reduced."""
    # Create list to store results
    setList = list()
    # Compare each header and create sets of headers that overlap.
    for item1 in myList:
        # Parse the first header to grab coordinate information.
        lane1 = readDict[item1]['lane']
        tile1 = readDict[item1]['tile']
        coord1 = numpy.array(readDict[item1]['x'],readDict[item1]['y'])
        item1Set = {item1}
        for item2 in myList:
            if item1 != item2: 
                # Don't compare a header to itself.
                lane2 = readDict[item2]['lane']
                tile2 = readDict[item2]['tile']
                coord2 = numpy.array(readDict[item2]['x'],readDict[item2]['y'])

                if lane1 == lane2 and tile1 == tile2:
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
    
def dupCounts(setList):
    """This function calculates various counts and returns a list of counts."""
    pdup_cnt = len(setList)
    opdup_cnt = 0
    opdupset_cnt = 0
    for item in setList:
        if len(item) > 1:
            opdup_cnt += len(item)
            opdupset_cnt += 1
    pdup_cnt -= opdupset_cnt
    return(pdup_cnt,opdup_cnt,opdupset_cnt)

def writeOutput(handle, myList):
    """Function to write Output"""
    handle.write(','.join(str(x) for x in myList) + "\n")


def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    # Read in the FASTQ file and return a dictionary with SEQ as the key and
    # list(HEADERS) as the values. Also return the total number of reads.
    logging.info("Starting to read in FASTQ file and creating dictionary.")
    mySeq, myRead = readFASTQ(args.fq,total_read)
    logging.info("Finished reading FASTQ file.")

    # Simple Counts 
    total_read =  len(myRead)    # total number of reads
    uniq_seq = len(mySeq)        # number of uniq sequences

    # Initialize unique read count
    uniq_read = 0

    # Loop through each sequence in the dictionary and examine it for
    # duplicates and optical duplicates.
    logging.info("Starting to identify duplicates.")
    for key in mySeq:
        # Copy the dictionary value using list() so that I don't modify the
        # dictionary value later in the script.
        myValue = list(mySeq[str(key)])   
        myLength = len(myValue)

        if myLength == 1:
            # If there is only one header for a given sequence, then this is a unique read.
            uniq_read += 1
        elif myLength > 1:
            # If there are more than one header for a given, then these are reads.
            total_dup = myLength
            # Identify if the duplicates are optical dups or something else.
            readSet = checkOptical(myValue, myRead, args.pix)
            (pdupCnt, opdupCnt, opdupSetCnt) = dupCounts(readSet)
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
