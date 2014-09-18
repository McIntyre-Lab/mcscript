#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys

# Add-on packages
import numpy as np
import gffutils

# McLab Packages
import mclib
from mclib import gff as mcgff

def getOptions():
    """ Function to pull in arguments """
    description = """ This script constructs wiggle plots and gene models.
    Wiggle plots are created from sorted BAM files. Gene models require a GFF
    file at this time.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--gff", dest="gffName", action='store', required=True, help="Name with PATH, to a GFF file. Note if the GFFutils database has not been created it will be, this could take around 30min [Required]")
    parser.add_argument("--sd", dest='strand', action='store_false', required=False, help="Flag if you want to make strand dependent fusions [Optional]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output PNG. [Required]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    #args = parser.parse_args()
    args = parser.parse_args(['--gff', '/home/jfear/Desktop/dmel-all-no-analysis-r5.51.gff', '-o', '/home/jfear/tmp/inr.bed'])
    return(args)

if __name__ == '__main__':
    args = getOptions()

    # Turn on logging
    mclib.logger.set_logger(args.log)
    # Output git commit version to log, if user has access
    mclib.git.git_to_log(__file__)

    ################################################################################
    # Initilize Variables and Files
    ################################################################################

    # Connect to the database
    flyGff = mcgff.FlyGff(args.gffName)

    # Connect to output file
    OUT = open(args.oname, 'w')

    # Initilize Counter
    cnt = 1

    # Initilize Suffix
    if args.strand:
        sfx = '_SI'
    else:
        sfx = '_SD'

    ################################################################################
    # Build Fusions
    ################################################################################

    # Get list of all genes in the genome
    genes = flyGff.get_genes()

    # Iterate over genes in the genome
    for gene in genes:
        # What strand is the current gene on 
        strand = gene.strand

        # Pull a list of exons for the current gene
        exons = flyGff.get_exons(gene)

        # Create a list of fusions by merging overlapping exons
        fusions = list(flyGff.db.merge(exons, ignore_strand = args.strand))

        # I need to figure out fusions are from a single exon 'singleton' and
        # which are multiple overlapping exons 'fusion'. Create a list of
        # coordinates for both exons and fusions. Fusion coordinate is in the
        # exons list then it is a singleton.
        exonTup = [(i.chrom, i.start, i.end) for i in exons]   # Make a list of (chrom, start, ends) for exons
        fusTup = [(i.chrom, i.start, i.end) for i in fusions]  # Make a list of (chrom, start, ends) for fusion

        for fusion in fusTup:
            if fusion in exonTup:
                # singleton
                name = "S{0}{1}".format(cnt, sfx)
            else:
                # fusion
                name = "F{0}{1}".format(cnt, sfx)

            # Write output in a bed format 
            chrom = str(fusion[0])
            start = str(fusion[1]) - 1      #remember bed is a 0-based format and gff is 1-based
            end = str(fusion[2])
            myOut = '\t'.join([chrom, start, end, name, '.', strand]) + "\n"
            OUT.write(myOut)

            # increment fusion counter
            cnt +=1

    OUT.close()
