#!/usr/bin/env python
import argparse 
import logging
import collections
import os.path
from Bio import SeqIO

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a fasta file and a 4-column bed file and calculates GC content.")
    parser.add_argument("-f", dest="fname", action='store', required=True, help="Name of reference file in FASTA format [Required]")
    parser.add_argument("-b", dest="bname", action='store', required=False, help="Name of 4-column BED file with regions of interest [Optional]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of output file in csv format [Required]")
    parser.add_argument("-g", dest="log", action='store', required=False, help="Name of log file [Optional]") 
    #args = parser.parse_args(['-f', '/home/jfear/tmp/hsv1_non-redundant.fa', '-b', '/home/jfear/tmp/test.bed', '-o', '/home/jfear/tmp/test_out.csv', '-g', '/home/jfear/tmp/test_out.log'])
    #args = parser.parse_args(['-f', '/home/jfear/tmp/hsv1_non-redundant.fa', '-o', '/home/jfear/tmp/test_out.csv', '-g', '/home/jfear/tmp/test_out.log'])
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')

def getGit():
    """ This function will parse the current Git commit version. This will
    allow recording of exactly which script was used to create a given
    output."""
    import subprocess
    fullname = os.path.abspath(__file__)
    gitdir = os.path.dirname(fullname)
    label = subprocess.check_output(["git", "--git-dir="+gitdir+"/.git", "--work-tree="+gitdir,"rev-parse","HEAD"])
    return(label.rstrip(), gitdir)

def readBed(bname):
    """ Read a 4-column BED file """
    bDict = collections.defaultdict(dict)
    with open(bname, 'r') as BED:
        for row in BED:
            chrom, rstart, rend, rid = row.rstrip().split('\t')
            bDict[chrom][rid] = [rstart, rend]
    return(bDict)

def getRegionGC(record, bDict):
    """ Calculate a regions GC content """
    regionGC = list()
    for rid, location in bedDict[record.id].items():
        rSeq = record.seq[int(location[0]):int(location[1])]
        region_G = rSeq.count("G")
        region_C = rSeq.count("C")
        region_GC = region_G + region_C
        region_len = len(rSeq)
        region_per = region_GC / float(region_len) * 100
        regionGC.append([rid, region_len, region_GC, region_per])
    return(regionGC)

def writeOutput(handle, myOut):
    """ Write output """
    handle.write(','.join([str(x) for x in myOut]) + "\n")

if __name__=='__main__':

    # Get command line options
    args = getOptions()

    # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)
        # Add version information from git to the log
        git_status, gitdir = getGit()
        logging.info("Starting %s", __file__) 
        logging.info("Running script from  %s", gitdir) 
        logging.info("Git commit id: %s", git_status)
    else:
        setLogger(os.devnull,logging.INFO)

    # If a bed file is provided read it to a dictionary. Create correct headers.
    if args.bname:
        logging.info("Bed file provided, building BED dictionary")
        bedDict = readBed(args.bname)
        myHeader = ['chrom', 'total_len', 'total_GC', 'total_per_GC', 'region', 'region_len', 'region_GC', 'region_per_GC']
    else:
        logging.info("No Bed file provided, only looking at genomic regions")
        bedDict = ''
        myHeader = ['chrom', 'total_len', 'total_GC', 'total_per_GC']


    with open(args.oname, 'w') as OUT:
        # Write header
        OUT.write(','.join([str(x) for x in myHeader]) + "\n")

        logging.info("Parsing FASTA file")
        with open(args.fname, 'r') as FA:
            # Loop through FASTA file
            for record in SeqIO.parse(FA, 'fasta'):
                tot_G = record.seq.count("G")
                tot_C = record.seq.count("C")
                tot_GC = tot_G + tot_C
                total = len(record.seq)
                tot_per_gc = tot_GC / float(total) * 100

                totList = [record.id, total, tot_GC, tot_per_gc]

                if bedDict:
                    regions = getRegionGC(record, bedDict)
                    for region in regions:
                        myOut = totList + region
                        writeOutput(OUT,myOut)
                else:
                    writeOutput(OUT,totList)

    logging.info("Script complete.")
