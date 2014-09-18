#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys

# Add-on packages
import numpy as np

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
    parser.add_argument("--sd", dest='strand', action='store_true', required=False, help="Flag if you want to make strand dependent fusions. The default is to make strand independent fusions [Optional]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output PNG. [Required]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    args = parser.parse_args()
    #args = parser.parse_args(['--gff', '/home/jfear/Desktop/dmel-all-no-analysis-r5.51.gff', '-o', '/home/jfear/Desktop/fb551.bed'])
    return(args)

def writeOutput(chrom, fusions, strand, sfx, OUT):
    """ This functions names the fusions and writes their output. If a fusion
    is made from a single exon then it is a singleton and will be prefixed with
    the letter 'S'. If a fusion is made up of overlapping fusions then it is
    prefixed with a 'F'.

    Arguments:
    chrom (str) = the current chromosome id
    exonList (list) = a list of exons to be merged
    fusionList (list) = a list of merged exons
    strand (str) = the strand to output; {'-', '+'} for SD and {'.'} for SI fusions
    sfx (str) = the suffix to append on to the fusion id {'_SI', '_SD'}
    OUT (obj) = File output object
    """

    # Attach the global counter
    global cnt

    for fusion in fusions:
        # Name the Fusion based on if it is a singleton or fusion.
        if fusion['merged']:
            # singleton
            name = "S{0}{1}".format(cnt, sfx)
        else:
            # fusion
            name = "F{0}{1}".format(cnt, sfx)

        # Write output in a bed format 
        start = str(fusion['start'] - 1)      #remember bed is a 0-based format and gff is 1-based
        end = str(fusion['end'])
        myOut = '\t'.join([chrom, start, end, name, '.', strand]) + "\n"
        OUT.write(myOut)

        # increment fusion counter
        cnt +=1

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
    logging.info('connecting to the gff file')
    flyGff = mcgff.FlyGff(args.gffName)

    # Connect to output file
    OUT = open(args.oname, 'w')

    # Initialize Counter
    cnt = 1

    ################################################################################
    # Build Fusions
    ################################################################################

    # Get list of chromosomes
    chromosomes = flyGff.get_chrom()

    for chrom in chromosomes:
        if args.strand:
            # Create Strand dependent fusions
            for currStrand in ('-', '+'):
                # Get list of all genes in the genome
                logging.info('Pulling list of exons for chromosome: {0} on strand {1}'.format(chrom.id,currStrand))
                exons = list(flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end),strand = currStrand))

                # Create a list of fusions by merging overlapping exons
                logging.info('Merging overlapping exons on strand '+currStrand)
                fusions = list(flyGff.merge(exons, ignore_strand=False))
                writeOutput(chrom.id, fusions, currStrand, '_SD', OUT)
        else:
            # Create Strand independent fusions
            # Get list of all genes in the genome
            logging.info('Pulling list of exons for chromosome: '+chrom.id)
            exons = list(flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end)))

            # Create a list of fusions by merging overlapping exons
            logging.info('Merging overlapping exons')
            fusions = list(flyGff.merge(exons, ignore_strand=True))
            writeOutput(chrom.id, fusions, '.' , '_SI', OUT)

    OUT.close()
