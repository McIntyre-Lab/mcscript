#!/usr/bin/env python
import argparse
import logging
from Bio import SeqIO
import itertools
import collections
import re
import numpy
import operator

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool for identifying duplicates and creating various useful output')
    parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A list of fq file [Required]')
    parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
    parser.add_argument('--pixel', dest='pix', action='store', default=100, required=False, help='Number of pixels to consider a read as an optical duplicate [Default:100]')
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
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            myDict[str(record.seq)].append(record.name)
            count += 1
    logging.info("Finished reading the FASTQ file.")
    return(myDict,count)


def checkOptical(myList, pix):
    """Given a list of FASTQ headers, this function will check if reads are
       within a certain 'args.pix' distance of each other. If there are within
       this distance, then they will be called optical duplicates."""

    logging.info("Starting check for optical dups.")

    # Create a dictionary to store optical duplicate information.
    myDict = dict()
    cntPcrDups = 0

    # Create a list of sets, where each set contains headers that are within
    # args.pix of each other.
    listOfSets, cntPcrDups = identifySets(myList, pix, cntPcrDups)

    # reduce the set list so that each set is a group of optical duplicates.
    redSetList = reduceSet(listOfSets)

    logging.info("Finished check for optical dups.")
    return(redSetList,cntPcrDups)


def identifySets(myList, pix, count):
    """This function steps through a list of headers and creates sets of those
       that fall within the args.pix range. These need to be reduced."""

    logging.info("Creating sets of optical dups.")
    # Create list to store results
    setList = list()
    # Compare each header and create sets of headers that overlap.
    for item1 in myList:
        # Parse the first header to grab coordinate information.
        match1 = parseHeader(item1)
        item1Set = {item1}
        for item2 in myList:
            # Don't compare a header to itself.
            if item1 != item2: 
                # Parse the second header to grab coordinate information.
                match2 = parseHeader(item2)

                # Test if headers are on the same lane and tile. If they are
                # then determine if they are within args.pix of each other
                # using numpy's norm function.
                if match1[0] == match2[0] and match1[1] == match2[1]:
                    myCoord = [numpy.array([int(match1[2]),int(match1[3])]),numpy.array([int(match2[2]),int(match2[3])])]       # create a list of numpy arrays [(x1,y1),(x2,y2)]
                    eucDist = numpy.linalg.norm(myCoord[0]-myCoord[1])
                    if eucDist <= pix:
                        item1Set.add(item2)
                    else:
                        count += 1
        setList.append(item1Set)
    logging.info("Finished creating sets of optical dups.")
    return(setList, count)


def parseHeader(fqName):
    """Function to parse the FASTQ header line"""
    match = re.search('.*:?([0-9]):([0-9]+):([0-9]+):([0-9]+).*/?([1-2]*)',fqName)
    matchList = filter(None,match.groups())     # Removes any empty strings, ie if there is no PE information
    return(matchList)


def reduceSet(setList):
    """Function to step through a list of sets and combine sets that overlap.
       This will create a unique list of sets, where each set contains the headers
       of reads that are within args.pix of each other"""

    logging.info("Starting to reduce optical dup sets.")
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

    logging.info("Finished reducing optical dup sets.")
    return(setList2)
    

def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    total_read = 0     # Initialize a counter to keep track of the total number of reads
    uniq_read = 0      # Initialize a counter to keep track of the total number of unique reads
    pcrDup = 0

    # Read in the FASTQ file and return a dictionary with SEQ as the key and
    # list(HEADERS) as the values. Also return the total number of reads.
    myFASTQ, total_read = readFASTQ(args.fq,total_read)

    # Calculate the number of unique SEQUENCES. 
    #NOTE: this is not the same as number of unique READS.
    uniq_seq = len(myFASTQ)

    # Loop through each sequence in the dictionary and examine it for
    # duplicates and optical duplicates.
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
            redSet, cntPcrDup = checkOptical(myValue, args.pix)
            pcrDup += cntPcrDup
        else:
            logging.ERROR("There is something wrong with the fastq dictionary at key %s" % key)

    print redSet, total_read, pcrDup, total_dup, uniq_read

if __name__ == '__main__':
    main()



    """

### Open fq
with open(args.fq,'rb') as input_fq:
# Get the list of all of the sequences
# Sequences are every 4 lines starting with line 1
   seqs=itertools.islice(input_fq,1,None,4)
# Get the count for each sequence
   counts=collections.Counter(seqs)
   total_num=sum(counts.values())
   uniq_num=len(counts)

percent_uniq = float(uniq_num) / (total_num)
with open(args.out,'wb') as dataout:
    dataout.write('total # fq is '+str(total_num)+'\n# unique sequences is '+str(uniq_num)+'\npercent unique is '+str(percent_uniq))

if args.table:
# Sort by count in descending order
        sorted_counts=sorted(counts.iteritems(),key=operator.itemgetter(1),reverse=True)

   	with open(args.table,'wb') as tableOut:
   	    for item in sorted_counts:
                tableOut.write(str(item[1])+' '+str(item[0]).strip()+'\n')
 
"""
