#!/usr/bin/env python

import sys 
import csv 
import argparse 
import collections
import logging
import numpy as np
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
csv.field_size_limit(100000000)

def getOptions():
    """ Function to pull in arguments """

    description="""This script can be used to calculates coverage (RPKM and APN) in two different settings:\n
    (1) Coverage can be calculated across an entire genomic region. To do
        this a 3-column bed file must be provided (Try fasta2bed.py). 
        col[0] = chromosome/fusion name (eg., chr2L or S7_SI)
        col[1] = start position (i.e., '0') 
        col[2] = end position (i.e., length) """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--mpileup", dest="mname", action='store', required=True, help="mpileup file [Required]",metavar="MPILEUP_FILE")
    parser.add_argument("-s", "--sam", dest="sname", action='store', required=True, help="sam alignment file [Required]", metavar="SAM_FILE")
    parser.add_argument("-b", "--bed", dest="bname", action='store', required=True, help="bed file (3 or 4 columns) [Required]", metavar="BED_FILE")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Out File", metavar="OUT_FILE")
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

# BED Functions
def read_bed(bname):
    """ Read BED file and create a Dictionary containing all information """
    bdict = collections.defaultdict(dict)
    with open(bname,'r') as BED:
        reader = csv.reader(BED,delimiter='\t')
        for row in reader:
            if len(row) == 3:
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])

                bdict[chrom]['start'] = start
                bdict[chrom]['end'] = end
                bdict[chrom]['region_length'] = end - start 
                bdict[chrom]['depth'] = 0
                bdict[chrom]['reads_in_region'] = 0
                bdict[chrom]['apn'] = 0
                bdict[chrom]['rpkm'] = 0
                bdict[chrom]['mean'] = 0
                bdict[chrom]['std'] = 0
                bdict[chrom]['cv'] = 0
                bdict[chrom]['counts'] = [0] * (end - start + 1)
            else:
                logging.error("I can only handle 3 column bed files.")
                exit(-1)

    return(bdict)

# SAM Functions
def read_sam(sname):
    """ Read SAM file to get read length and number of mapped reads. Note: if
        you have allowed ambiguous mapping then reads are counted multiple times.  """
    num_mapped_reads = 0
    read_length = 0
    with open(sname,'r') as SAM:
        for row in SAM:
            if not row.startswith('@'):
                record = row.strip().split('\t')
                if record[1] != 4:
                    num_mapped_reads += 1
                read_length = max(read_length,len(record[9]))

    return(num_mapped_reads,read_length)

# MPILEUP Functions
def read_mpileup(mname,bdict):
    """ Read mpileup and store depth and length into a new Dictionary """
    with open(mname, 'r') as MPILEUP:
        reader = csv.reader(MPILEUP, delimiter='\t',quoting=csv.QUOTE_NONE)
        for row in reader:
            mchrom = row[0]
            mpos = int(row[1])
            mdepth = int(row[3])
            if bdict.has_key(mchrom):
                bdict[mchrom]['counts'][mpos-1] = mdepth
                bdict[mchrom]['depth'] += mdepth


# Coverage Functions
def calc_coverage(bdict,num_mapped_reads,read_length):
    """ Calculate different coverage metric: Estimate number of reads in
    region.  Average per nucleotide coverage (APN).  Reads per kilobase per
    million mapped reads (RPKM). """

    for key in bdict:
        depth = bdict[key]['depth']
        counts = bdict[key]['counts']
        region_length = bdict[key]['region_length']
        reads_in_region = float(depth) / float(read_length) # Estimate the number of reads in region based on depth/read_length.
        if not depth == 0:
            bdict[key]['reads_in_region'] = reads_in_region 
            bdict[key]['apn'] = float(depth) / float(region_length) # Calculate average per nucleotide coverage APN (depth in region / length of region).
            bdict[key]['rpkm'] = (1000000000.0 * reads_in_region) / (num_mapped_reads * region_length) # Calculate reads per kilobase per million mapped reads RPKM from Moretzavi et al. 

            if depth != sum(counts):
                logging.error("There was a discrepancy between the depth and the count. Must be a bug!")

            mymean = np.mean(counts)
            mystd = np.std(counts)
            bdict[key]['mean'] = mymean
            bdict[key]['std'] = mystd
            bdict[key]['cv'] = mystd / mymean

# Output Functions
def writeOutput(oname,bdict,num_mapped_reads,read_length):
    """ I tried writing output using the CSV module, but this did not behave
        well with SAS downstream. So I opted for the brute force method. """

    header = ['fusion_id','mapped_reads','read_length','region_length','region_depth','reads_in_region','apn','rpkm','mean','std','cv']
    with open(oname, 'wb') as OUT:
        OUT.write(','.join(header) + '\n')
        for key in bdict:
            OUT.write(','.join(str(x) for x in [key,num_mapped_reads,read_length,bdict[key]['region_length'],bdict[key]['depth'],bdict[key]['reads_in_region'],bdict[key]['apn'],bdict[key]['rpkm'],bdict[key]['mean'],bdict[key]['std'],bdict[key]['cv']]) + '\n')

def main():
    """ MAIN Function to execute everything """
    args = getOptions()

    # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)

    # READ IN BED FILE
    logging.info("Reading the BED File '%s'." % args.bname)
    bdict = read_bed(args.bname)
    logging.info("Finished reading the BED File '%s'." % args.bname)

    # Use SAM file to estimate the number of mapped reads and the max read length
    logging.info("Reading the SAM File '%s'." % args.sname)
    num_mapped_reads, read_length = read_sam(args.sname)
    logging.info("Finished reading the SAM File '%s'." % args.sname)

    # Read Mpileup file
    logging.info("Reading the Mpileup File '%s'." % args.mname)
    read_mpileup(args.mname,bdict)
    logging.info("Finished reading the Mpileup File '%s'." % args.mname)

    logging.info("Performing Calculations")
    calc_coverage(bdict,num_mapped_reads,read_length)    # Use all of this information to estimate the number of reads in genomic region and calculate coverage (APN,RPKM) for the genomic region
    logging.info("Finished Calculations")

    logging.info("Writing Output")
    writeOutput(args.out,bdict,num_mapped_reads,read_length) # Write output to CSV file
    logging.info("Finished writing Output")

if __name__=='__main__':
    main()
    logging.info("Script Complete")
