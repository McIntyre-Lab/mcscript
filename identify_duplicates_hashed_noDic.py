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
import hashlib
import time
import pdb
import scipy as sp
import numpy as np
import timeit

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

def filterOnFly(f, stride):
    for i, line in enumerate(f):
        if (i % stride == 0) or (i % stride == 1):
            yield line
def readFASTQ2(fname):
    """Read a fastq file and store information into two different dictionaries.
    The first is mySeqDict which has the sequence as the key and a list of
    headers with that have that sequence as the values. The second is
    myReadDict, which has read name as the key and coordinate infromation as
    values."""
    logging.info("Reading the FASTQ file.")

    with open(fname) as f:
        data = np.genfromtxt(filterOnFly(f,4), dtype = 'string')
    mySeqDict = sp.zeros((data.shape[0] / 2,7), dtype = '|S95')
    mySeqDict[:,0]   = data[1::2]
    #TODO: This can be done more efficiently through marix operation (numpy.coredefcahrarray or so)
    mySeqDict[:,1:6] = sp.array([parseHeader(x) for x in data[::2]])
    mySeqDict[:,6]   = data[::2]
    logging.info("Finished reading the FASTQ file.")
    return mySeqDict

def readFASTQ(fname):
    """Read a fastq file and store information into two different dictionaries.
       The first is mySeqDict which has the sequence as the key and a list of
       headers with that have that sequence as the values. The second is
       myReadDict, which has read name as the key and coordinate infromation as
       values."""
    logging.info("Reading the FASTQ file.")
    #mySeqDict = collections.defaultdict(list)
    #myReadDict = dict()
    mySeqDict = []
    #TODO pre-init so that i do not have to re-allocate every time
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            # Create mySeqDict
            tmp = [str(record.seq)]

            # Parse Header and create myReadDict
            match = parseHeader(record.name)
            if len(match) == 4: # if there is no PE read information, append a 1
                match.append(1)
                #myReadDict[record.name] = {'lane':match[0], 'tile':match[1],'coord':(match[2],match[3]),'read':match[4]}
            match.append(record.name)
            tmp.extend(match)
            mySeqDict.append(tmp)
    logging.info("Finished reading the FASTQ file.")
    return sp.array(mySeqDict)


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
    matchIntList = [x for x in matchList]
    if len(matchIntList) ==4:
        matchIntList.append(1)
    return(matchIntList)


def checkOptical(myList, pix):
    """Given a list of FASTQ headers, this function will check if reads are
       within a certain 'args.pix' distance of each other. If reads are within
       this distance, then they will be called optical duplicates."""
    # Create a storage dictionary

    # Create a list of sets, where each set contains headers that are optical duplicates.
    listOfSets = identifySets(myList, pix)

    # reduce the set list so that each set is a group of optical duplicates.
    redSetList = reduceSet(listOfSets)
    return(redSetList)
   

def identifySets(myList, pix):
    """This function steps through a list of headers and creates sets of those
       that fall within the args.pix range. The resulting sets may overlap so
       they need reduced."""
    # Create list to store results
    setList = list()

    # Compare each header and create sets of headers that overlap.
    for i in xrange(myList.shape[0]):
        # Grab coordinate information from readDict
        item1Set = {myList[i,6]}
        lane1 = float(myList[i,1])
        tile1 = float(myList[i,2])
        coord1 = (float(myList[i,3]),float(myList[i,4]))

        lane2 = sp.array(myList[i+1:,1], dtype = 'float')
        tile2 = sp.array(myList[i+1:,2], dtype = 'float')
        iLane = lane1 == lane2
        iTile = tile1 == tile2
        c_x   = sp.array(myList[i+1:,3], dtype = 'float')
        c_y   = sp.array(myList[i+1:,4], dtype = 'float')

        iSel  = sp.sqrt((coord1[0]-c_x)**2 + (coord1[1]-c_y)**2 ) < pix
        if sp.sum(iSel & iLane & iTile) != 0:
            setList.append(item1Set.union(set(myList[i+1:,6][iSel].tolist())))
        else:
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
    start = time.time()
    mySeq = readFASTQ2(args.fq)
    print "Took %f seconds to read" % (time.time() - start)
    logging.info("Finished reading FASTQ file.")
    start = time.time()
    # Simple Counts 
    total_read =  mySeq.shape[0]    # total number of reads
    mySeq = mySeq[sp.argsort(mySeq[:,0]),:]
    uqSeq, uidx, uinv  = sp.unique(mySeq[:,0], return_index = True, return_inverse = True)
    uniq_seq = uqSeq.shape[0]        # number of uniq sequences
    print "Took %f seconds to unique" % (time.time() - start)


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

    start = time.time()
    for i in xrange(uidx.shape[0] - 1):
        # Copy the dictionary value using list() so that I don't modify the
        # dictionary value later in the script.
        #    uqSeq[0] corresponds to mySeq[uinv[uidx[0] == uinv]]
        #myLength = len(myValue)
        if i % 1000 == 0:
            print "Took %f seconds for 1000 entries out of %i" % (time.time() - start, uqSeq.shape[0])
            #pdb.set_trace()
            start = time.time()
            
        myValue  = mySeq[range(uidx[0],uidx[1])]
        myLength = myValue.shape[0]

        if myLength == 1:
            # If there is only one header for a given sequence, then this is a unique read.Is that not the same as uniq_seq?
            uniq_read += 1
        else:
            # If there are more than one header for a given, then these are reads.
            total_dup = myLength
            # Identify if the duplicates are optical dups or something else.
            readSet = checkOptical(myValue, args.pix)
            (pdupCnt, opdupCnt, opdupSetCnt, dropOp) = dupCounts(readSet,pdupCnt, opdupCnt, opdupSetCnt, args.fqOut, dropOp)

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
