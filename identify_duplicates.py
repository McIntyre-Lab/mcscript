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
    parser.add_argument('-t','--table',dest='table',action='store',required=False, help='Output table of intermediate size information')
    parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
    args = parser.parse_args()
    return(args)


def setLogger(fname,loglevel):
    """Function for handling error logging"""
    logging.basicConfig(filename=fname, level=loglevel, format='%asctime)s - %(levelname)s - %(message)s')


def readFASTQ(fname):
    """Read a fastq file and store information in a dictionary"""
    myDict = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(fname,'r') as FQ:
        for record in SeqIO.parse(FQ, 'fastq'):
            myDict[str(record.seq)]['name'].append(record.name)
            myList = parseHeader(record.name)
            if len(myList) == 5:
                lane, tile, xcord, ycord, pe = myList
            elif len(myList) == 4:
                lane, tile, xcord, ycord = myList
                pe = 0
            else:
                logging.ERROR("malformed FASTQ header %s" % (record.name))

            if pe != 2:
                myDict[str(record.seq)]['lane'].append(lane)
                myDict[str(record.seq)]['tile'].append(tile)
                coords = map(int,(xcord,ycord))
                myDict[str(record.seq)]['coord'].append(numpy.array(coords))
    return(myDict)


def parseHeader(fqName):
    """Function to parse the FASTQ header line"""
    match = re.search('.*:?([0-9]):([0-9]+):([0-9]+):([0-9]+).*/?([1-2]*)',fqName)
    matchList = filter(None,match.groups())     # Removes any empty strings, ie if there is no PE information
    return(matchList)


def calcDist(myCoord)
    """Function to calculate the distance between all pairwise combinations of coordinates"""
    myList = list()
    myIter = itertools.combinations(myCoord,2)     # Create an iterator of all possible pairwise combinations
    for record in myIter:
        myList.append(numpy.linalg.norm(record[0]-record[1]))    # Calculate the distance between each set of coordinates using numpy's norm function
    return(myList)
    

def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
        setLogger(args.log,logging.ERROR)
        setLogger(args.log,logging.WARNING)

    myFASTQ = readFASTQ(args.fq)



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
