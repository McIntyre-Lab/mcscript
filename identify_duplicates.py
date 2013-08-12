#!/usr/bin/env python
import cProfile
import argparse
import logging
from Bio import SeqIO
import collections
import re
import math
import os
import itertools

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool for identifying duplicates and creating various useful output')
    parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A list of fq file [Required]')
    parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
    parser.add_argument('--pixel', dest='pix', action='store', default=100, required=False, help='Number of pixels to consider a read as an optical duplicate [Default:100]')
    parser.add_argument('-a', dest='append', action='store_true', help='This flag will cause the output dataset to be appended too.')
    parser.add_argument('-t','--table', dest='table', action='store', required=False, help='Output table with the duplicate counts for each uniq sequence')
    parser.add_argument('-f','--fqOut', dest='fqOut', action='store', required=False, help='Output a FASTQ file optical duplicates are reduced')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
    args = parser.parse_args()
    return(args)


def setLogger(fname,loglevel):
    """Function for handling error logging"""
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')


def readFASTQ(fname):
    """Read a fastq file and store information into two different dictionaries.
       The first is mySeqDict which has the sequence as the key and a list of
       headers with that have that sequence as the values. The second is
       myReadDict, which has read name as the key and coordinate infromation as
       values."""
    logging.info("Reading the FASTQ file.")
    mySeqDict = collections.defaultdict(list)
    myReadDict = dict()
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            # Create mySeqDict
            mySeqDict[str(record.seq)].append(record.name)

            # Parse Header and create myReadDict
            match = parseHeader(record.name)
            if len(match) == 4: # if there is no PE read information, append a 1
                match.append(1)
            myReadDict[record.name] = {'lane':match[0], 'tile':match[1],'coord':(match[2],match[3]),'read':match[4]}

    logging.info("Finished reading the FASTQ file.")
    return(mySeqDict,myReadDict)


def parseHeader(fqName):
    """Function to parse the FASTQ header line. This does a simple regex to
       pull out different parts of the header line. Of importance are:
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
       within a certain 'args.pix' distance of each other. If reads are within
       this distance, then they will be called optical duplicates."""
    # Create a storage dictionary
    myDict = dict()

    # Create a list of sets, where each set contains headers that are optical duplicates.
    listOfSets = identifySets(myList, readDict, pix)

    # reduce the set list so that each set is a group of optical duplicates.
    redSetList = reduceSet(listOfSets)
    return(redSetList)


def identifySets(myList, readDict, pix):
    """This function steps through a list of headers and creates sets of those
       that fall within the args.pix range. The resulting sets may overlap so
       they need reduced."""
    # Create list to store results
    setList = list()

    # Compare each header and create sets of headers that overlap.
    for index, item1 in enumerate(myList):
        # Grab coordinate information from readDict
        item1Set = {item1}
        lane1 = readDict[item1]['lane']
        tile1 = readDict[item1]['tile']
        coord1 = readDict[item1]['coord']

        for item2 in myList[index+1:]:
            # Grab coordinate information from readDict
            lane2 = readDict[item2]['lane']
            tile2 = readDict[item2]['tile']
            coord2 = readDict[item2]['coord']

            if lane1 == lane2 and tile1 == tile2:
                eucDist = math.sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2 )
                if eucDist <= pix:
                    item1Set.add(item2)
        setList.append(item1Set)
    return(setList)


def reduceSet(setList):
    """Step through a list of sets and combine sets that overlap.
       This will create a unique list of sets, where each set contains the headers
       of reads that optical duplicates."""
    setList2 = [setList[0]]
    for item1 in setList[1:]:
        inSetList2 = False
        for item2 in setList2:
            if item1 & item2:       # If there is an intersection
                item2 |= item1      # Combine sets and go to next item in setList
                inSetList2 = True
                break
        if not inSetList2:          # If I could not find any overlaps then append set to setList2
            setList2.append(item1)
    return(setList2)
    

def dupCounts(setList,pdup_cnt,opdup_cnt,opdupset_cnt,flagFQ,dropOp):
    """This function calculates various counts and returns a list of counts."""
    dupset_num = 0
    for item in setList:
        if len(item) > 1:
            opdup_cnt += len(item)
            opdupset_cnt += 1
            dupset_num += 1
            if flagFQ:
                # If the user asks for FASTQ output, store a set of optical
                # duplicate headers. I will keep one optical duplicate from
                # each set, to create a reduced list and not completely remove
                # reads that are optical duplicates.
                myList = list(item)
                dropOp |= set(myList[1:])

    pdup_cnt += len(setList) - dupset_num
    return(pdup_cnt,opdup_cnt,opdupset_cnt,dropOp)


def writeOutput(handle, myList):
    """Function to write output from a list to a csv"""
    handle.write(','.join(str(x) for x in myList) + "\n")


def writeCount(args,myOutHeader,myOut):
    """Write a summary table of counts. If the append flag is added, check if
       the output file exists and try to append to it."""
    if args.append:
        if os.path.exists(args.out):
            try:
                with open(args.out, 'a') as handle:
                    writeOutput(handle,myOut)
            except:
                logging.error("Could not open output file, it must be busy.")
        else:
            try:
                with open(args.out, 'w') as handle:
                    writeOutput(handle,myOutHeader)
                    writeOutput(handle,myOut)
            except:
                logging.error("Could not open output file, do I have write permission?")
    else:
        try:
            with open(args.out, 'w') as handle:
                writeOutput(handle,myOutHeader)
                writeOutput(handle,myOut)
        except:
            logging.error("Could not open output file, do I have write permission?")


def writeTable(table, mySeqDict):
    """Write a table with how many time each sequence is duplicated."""
    myDict = dict()
    for seq in mySeqDict:
        myDict[seq] = len(mySeqDict[seq])

    with open(table,'w') as handle:
        for item in sorted(myDict,key=myDict.get,reverse=True):
            writeOutput(handle, [myDict[item],item])

def writeFQout(fqIn,fqOut,dropOp):
    """Output a FASTQ file that reduces the optical duplicates, so that only a
    single read is left from an optical duplicate set."""
    with open(fqIn, 'r') as FQ:
        with open(fqOut, 'w') as FO:
            for record in SeqIO.parse(FQ, 'fastq'):
                if not {record.name} & dropOp:
                    SeqIO.write(record,FO,"fastq")


def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    ## READ IN DATA ##
    # Read in the FASTQ file and return a dictionary with SEQ as the key and
    # list(HEADERS) as the values. Also return the total number of reads.
    logging.info("Starting to read in FASTQ file and creating dictionary.")
    mySeq, myRead = readFASTQ(args.fq)
    logging.info("Finished reading FASTQ file.")

    # Simple Counts 
    total_read =  len(myRead)    # total number of reads
    uniq_seq = len(mySeq)        # number of uniq sequences

    # Initialize counts
    uniq_read = 0
    pdupCnt = 0 
    opdupCnt = 0 
    opdupSetCnt = 0

    # Initialize a set to store headers of reads I may want to drop from my
    # FASTQ file
    dropOp = set()

    ## Loop Through Data and identify duplicates ##
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
            (pdupCnt, opdupCnt, opdupSetCnt, dropOp) = dupCounts(readSet,pdupCnt, opdupCnt, opdupSetCnt, args.fqOut, dropOp)
        else:
            logging.error("There is something wrong with the fastq dictionary at key %s" % key)
    logging.info("Finished identifying duplicates.")

    ## Write desired output ##
    per_uniq_read = round(float(uniq_read)/total_read * 100, 2)
    myOutHeader= ["File_Name","Total_Reads", "Num_Unique_Reads", "Per_Uniq_Reads", "Num_Unique_Seq", "Num_PCR_Dups", "Num_OP_Dups", "Num_OP_Dup_Sets"]
    myOut = [args.fq,total_read, uniq_read, per_uniq_read, uniq_seq, pdupCnt, opdupCnt, opdupSetCnt]

    # Test if we want to append the file, if so handle the headers correctly.
    logging.info("Writing summary counts table.")
    writeCount(args,myOutHeader,myOut)

    if args.table:
        logging.info("Writing sequence counts table.")
        writeTable(args.table,mySeq)

    if args.fqOut:
        logging.info("Writing reduced fastq file.")
        writeFQout(args.fq, args.fqOut,dropOp)
    logging.info("Script Finished.")


if __name__ == '__main__':
    main()
