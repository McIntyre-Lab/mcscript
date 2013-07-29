#!/usr/bin/env python

import sys,os,csv,argparse,logging,pprint,re,datetime,itertools
from datetime import datetime
csv.field_size_limit(1000000000)
 
def getOptions():
    """Function to pull in arguments from the command line"""

    description="""This script compares the number of reads in a given input fastq file and the combined number of reads of the given output files to determine if the alignment output is the correct size. """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_name", dest="input", action='store', required=True, help="input file")
    parser.add_argument('-a', dest="alignedflag", action='store_true',required=False)
    parser.add_argument("-unaln", "--unaligned_output", dest="unaligned_output", action='store', required=True, nargs = '*', help="unaligned output file list [Required]")
    parser.add_argument("-aln", "--aligned_output", dest="aligned_output", action='store', required=True, nargs = '*', help="aligned output file list [Required]")
    args = parser.parse_args()
    return(args)


def fileLength(filename):
    """Get the number of lines in the given file"""
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getExtension(whole_filename):
    """Get the extension of the given file"""
    fileName,fileExtension=os.path.splitext(whole_filename)
    return fileExtension

def getNumSAMReads(filename,aligned=0):
    """Get the number of reads in a given SAM file. Specify whether we want aligned reads (aligned=1) or unaligned reads (aligned=0)"""
    unaligned=["4","77","141"] #the values for unaligned reads
    unalignedCount=0
    alignedCount=0
    alreadySeenDict={} #To keep track of reads that have been seen, so reads mapping more than once are only counted once
    with open(filename) as myFile:
        myCSV=csv.reader(myFile,delimiter='\t')
        for row in myCSV:
            if row[0][0]=="@":
                continue	#skip any rows that start with @, since that means it's a header
            if row[0] in alreadySeenDict:
                continue #skip duplicate rows
            alreadySeenDict[row[0]]=True
            if row[1] in unaligned: #indicating the row represents an unaligned read
                unalignedCount=unalignedCount+1
            else:
                alignedCount=alignedCount+1
    
    if aligned:
        return alignedCount
    else:
        return unalignedCount

def getNumMAFReads(filename):
    """Get the number of reads in a given MAF file. We assume the maf file consists of paragraphs with a first line that starts with "a ", a second line that starts with "s " and has the name of the reference, and a third line that starts with an "s " and tells the name of the actual read we are mapping to the reference. We remove all duplicate reads by searching in this third line and building/checking against the alreadySeenDict."""

    count=0
    alreadySeenDict={}
    with open(filename) as fileRead:
        fileRead,fileRead2=itertools.tee(fileRead)
        next(fileRead2,None)
        for linePair in itertools.izip(fileRead,fileRead2):
            if linePair[0][0]=='s' and linePair[1][0]=='s':
                line_read = (linePair[1].split(" ")[1])
                if line_read in alreadySeenDict:
                    pass #skip all duplicates
                else:
                    count=count+1 
                    alreadySeenDict[line_read]=True
         
    return count

def getNumFASTQReads(filename):
    return (fileLength(filename))/4


def getNumFastaReads(filename):
    count=0
    with open(filename) as myFile:
        for line in myFile:
            if line[0]=='>':
                count=count+1
    return count

def getNumReads(filename,aligned=0):
    """Get the number of reads in the file"""
    extension=getExtension(filename).lower() #we can ignore any capitalization
    if extension==".fq" or extension == '.txt': #remove this later. txt should not be anything.
        return getNumFASTQReads(filename)
    elif extension==".sam":
        return getNumSAMReads(filename,aligned)
    elif extension==".maf":
        return getNumMAFReads(filename)
    elif extension==".fasta" or extension==".fa":
        return getNumFastaReads(filename)

    else: 
        print("-1") #Denoting the size checker was given a file other than a maf, sam fastq,or fasta file. 

    

def main():
 #   printOutFile=open("/bio/mcintyre/pipeline_dev/scripts/printoutfile.txt","wb")
    args=getOptions()
    input_file=args.input
    if args.alignedflag:
        aligned=1
    else:
        aligned=0
    print("Input:\t\t"+str(datetime.now()))
    input_reads=getNumReads(input_file,aligned)
    print("Input reads: "+str(input_reads))
    output_size_list=[]
    for output_file in args.unaligned_output:
        print("Unaligned:\t"+str(datetime.now()))
        output_size_list.append(getNumReads(output_file,aligned=0))
    for output_file2 in args.aligned_output:
        print("Aligned:\t"+str(datetime.now()))
        output_size_list.append(getNumReads(output_file2,aligned=1))
    print("Done:\t\t"+str(datetime.now()))
    print("Output reads: "+str(output_size_list))
    if sum(output_size_list)==input_reads:
        print("1") # A 1 is printed if the check was successful
    else:
        print("-1") # a 0 is printed if the check found a failure
  #  printOutFile.close()
""" If I am running this script directly then I will run the main function. If
I am importing this script into another script then I will just import all of
my functions """

if __name__ == '__main__':
    main()
