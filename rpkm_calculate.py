#!/usr/bin/env python

import sys, csv, argparse, logging
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
        col[2] = end position (i.e., length) 

    (2) Coverage can also be calculated by excising specific exons/fusions
        from a genome. For example if you have aligned to the genome, but want
        coverage at the exon level. For this a 4-column bed must be provided.
        col[0] = chromosome name (eg., chr2L)
        col[1] = exon/fusion start position (eg., 2929)
        col[2] = exon/fusion end position (eg., 3090)
        col[3] = exon/fusion name (eg., S7_SI) 
        IMPORTANT: Setting 2 requires a lot of RAM ~10-12GB for calculating
        coverage using fly fusions """

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
def init_fdict(label,start,end,fdict):
    """ Initialize fdict so that everything is printed to output even if there
        is no coverage information (i.e, print lines with all 0's) """

    fdict[label] = {}
    fdict[label]['depth'] = 0
    fdict[label]['region_length'] = end - start + 1
    fdict[label]['reads_in_region'] = 0
    fdict[label]['apn'] = 0
    fdict[label]['rpkm'] = 0

def three_col_bed(bdict,chrom,start,end,fdict):
    """ Read in information from a 3-column BED file and store it in bdict """

    if not fdict.has_key(chrom):        # As you read in BED information go ahead and initialize fdict with all of its 0 values
        init_fdict(chrom,start,end,fdict)

    if not bdict.has_key(chrom):        # Build bdict from BED file
        bdict[chrom] = {}
    bdict[chrom]['start'] = start
    bdict[chrom]['end'] = end

def four_col_bed(bdict,chrom,start,end,label,fdict):
    """ Read in information from a 4-column BED file and store in bdict. Note:
        for the 4-column setting I am expanding bdict to the coordinate level.
        In other words bdict[chrom][postion]... This was the quickest way to
        test if a mpileup position is in region, unfortunately it requires a
        lot of RAM ~10-12GB for fly fusions """

    if not fdict.has_key(label): # As you read in BED information go ahead and initialize fdict with all of its 0 values
        init_fdict(label,start,end,fdict)

    if not bdict.has_key(chrom):           # Build bdict from BED file
        bdict[chrom] = {}
    for pos in xrange(start,end+1):
        if not bdict[chrom].has_key(pos):
            bdict[chrom][pos] = {}
        bdict[chrom][pos]['start'] = start
        bdict[chrom][pos]['end'] = end
        bdict[chrom][pos]['label'] = label

def read_bed(args,fdict):
    """ Read BED file and create a Dictionary containing all information """
    if args.log:
        logging.info("Reading the BED File '%s'." % args.bname)

    bdict = {}
    with open(args.bname,'r') as BED:
        reader = csv.reader(BED,delimiter='\t')
        for row in reader:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            if len(row) == 4:                    # If BED file has 4 columns
                label = row[3]
                four_col_bed(bdict,chrom,start,end,label,fdict)
                flag_4 = 1
            elif len(row) == 3:
                three_col_bed(bdict,chrom,start,end,fdict)
                flag_4 = 0
            else:
                logging.error("I can only handle 3 or 4 column bed files. See Help for descriptions")
                exit(-1)

    return(bdict,flag_4)

# SAM Functions
def read_sam(args):
    """ Read SAM file to get read length and number of mapped reads. Note: if
        you have allowed ambiguous mapping then reads are counted multiple times.  """
    logging.info("Reading the SAM File '%s'." % args.sname)

    num_mapped_reads = 0
    read_length = 0
    with open(args.sname,'r') as SAM:
        for row in SAM:
            if not row.startswith('@'):
                record = row.strip().split('\t')
                if record[1] != 4:
                    num_mapped_reads += 1
                read_length = max(read_length,len(record[9]))

    return(num_mapped_reads,read_length)

# MPILEUP Functions
def flag4_mpileup(fdict,bdict,chrom,pos,depth):
    """ Sum depth over parts of a genomic region """

    if bdict[chrom].has_key(pos):
        start = bdict[chrom][pos]['start']
        end = bdict[chrom][pos]['end']
        label = bdict[chrom][pos]['label']
        if fdict.has_key(label):
            fdict[label]['depth'] += depth
        else:
            logging.error("fdict did not have a key found in bdict")
            exit(-1)


def noflag4_mpileup(fdict,bdict,chrom,depth):
    """ Sum depth over entire genomic region """

    if fdict.has_key(chrom):
        fdict[chrom]['depth'] += depth
    else:
        logging.error("fdict did not have a key found in bdict")
        exit(-1)

def read_mpileup(args,flag_4,bdict,fdict):
    """ Read mpileup and store depth and length into a new Dictionary """
    logging.info("Reading the Mpileup File '%s'." % args.mname)

    with open(args.mname, 'r') as MPILEUP:
        reader = csv.reader(MPILEUP, delimiter='\t',quoting=csv.QUOTE_NONE)
        for row in reader:
            mchrom = row[0]
            mpos = int(row[1])
            mdepth = int(row[3])
            if bdict.has_key(mchrom):
                if flag_4:
                    flag4_mpileup(fdict,bdict,mchrom,mpos,mdepth)
                else:
                    noflag4_mpileup(fdict,bdict,mchrom,mdepth)

# Coverage Functions
def calc_coverage(fdict,num_mapped_reads,read_length):
    """ Calculate different coverage metric:
        Estimate number of reads in region.
        Average per nucleotide coverage (APN).
        Reads per kilobase per million mapped reads (RPKM). """
    logging.info("Calculating Coverage Counts")

    for key in fdict:
        if not fdict[key]['depth'] == 0:
            fdict[key]['reads_in_region'] = (fdict[key]['depth']) / (read_length*1.0)                                              # Estimate the number of reads in region based on depth/read_length. Multiplying by 1.0 to tell python to use decimals.
            fdict[key]['apn'] = (fdict[key]['depth']) / (fdict[key]['region_length'] * 1.0)                                        # Calculate average per nucleotide coverage APN (depth in region / length of region). Multiplying by 1.0 to tell python to use decimals.
            fdict[key]['rpkm'] = (1000000000.0 * fdict[key]['reads_in_region']) / (num_mapped_reads * fdict[key]['region_length']) # Calculate reads per kilobase per million mapped reads RPKM from Moretzavi et al. 

# Output Functions

def writeOutput(args,fdict,num_mapped_reads,read_length):
    """ I tried writing output using the CSV module, but this did not behave
        well with SAS downstream. So I opted for the brute force method. """
    logging.info("Writing Output")

    header = ['fusion_id','mapped_reads','read_length','region_length','region_depth','reads_in_region','apn','rpkm']
    with open(args.out, 'wb') as OUT:
        OUT.write(','.join(header) + '\n')
        for key in fdict:
            OUT.write(','.join(str(x) for x in [key,num_mapped_reads,read_length,fdict[key]['region_length'],fdict[key]['depth'],fdict[key]['reads_in_region'],fdict[key]['apn'],fdict[key]['rpkm']]) + '\n')

def main():
    """ MAIN Function to execute everything """

    args = getOptions()                                  # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)

    fdict = {}                                           # Initialize the dictionary that will hold all of my information for output
    bdict, flag_4 = read_bed(args,fdict)                 # Read through BED file and create dictionary with BED information. Note: A 4 column BED will require a large amount of RAM ~10-12GB for FLY
    read_mpileup(args,flag_4,bdict,fdict)                # Read Mpileup file and add depth to fdict
    del bdict                                            # Remove bed information because I don't need it any more
    num_mapped_reads, read_length = read_sam(args)       # Use SAM file to estimate the number of mapped reads and the max read length
    calc_coverage(fdict,num_mapped_reads,read_length)    # Use all of this information to estimate the number of reads in genomic region and calculate coverage (APN,RPKM) for the genomic region
    writeOutput(args,fdict,num_mapped_reads,read_length) # Write output to CSV file


if __name__=='__main__':
    main()
