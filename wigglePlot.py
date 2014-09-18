#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys

# Add-on packages
import numpy as np
import matplotlib 
matplotlib.use('Agg')

# McLab Packages
import mclib
from mclib import bam as mcbam
from mclib import gff as mcgff
from mclib import wiggle as mcwiggle

# TODO: Add option to take coordinates instead of just a gene
# TODO: Make it so you can give a BED file instead of just a GFF
# TODO: Add ability to make overlapping wiggle when up to XX groups are given
# TODO: Add option to provide a FudgeFactor for correction
# Minor TODO: For flies, add the ability to give it a fusion file and plot the fusions with labels

def getOptions():
    """ Function to pull in arguments """
    description = """ This script constructs wiggle plots and gene models.
    Wiggle plots are created from sorted BAM files. Gene models require a GFF
    file at this time.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--gff", dest="gffName", action='store', required=True, help="Name with PATH, to a GFF file. Note if the GFFutils database has not been created it will be, this could take around 30min [Required]")
    parser.add_argument("--bam", nargs='*', dest="bamList", action='store', required=True, help="Name with PATH, to sorted BAM files that will be averaged. [Required]")
    parser.add_argument("-g", dest="geneName", action='store', required=True, help="Name of the gene you want to make wiggles of. [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output PNG. [Required]")
    parser.add_argument("--fudge", dest="ff", action='store_true', required=False, help="Use the built in fudge factor. [Optional]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    args = parser.parse_args()
    #args = parser.parse_args(['--gff', '/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff', '--bam','/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V1.sorted.bam','/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V2.sorted.bam', '-g', 'InR', '-o', '/home/jfear/tmp/inr.png'])
    return(args)

if __name__ == '__main__':
    args = getOptions()

    # Turn on logging
    mclib.logger.set_logger(args.log)
    # Output git commit version to log, if user has access
    mclib.git.git_to_log(__file__)

    ################################################################################
    # GENE ANNOTATION 
    ################################################################################

    ## Import GFF database
    myGffDb = mcgff.FlyGff(args.gffName)

    ## Pull gene from database
    myGene = mcgff.FlyGene(args.geneName, myGffDb)

    ## Create Gene Model
    myModel = mcwiggle.GeneModel(myGene)

    ################################################################################
    # GENE COVERAGE
    ################################################################################

    # Pull in bam file and make gene pileup
    pileups = []
    for bam in args.bamList:
        currBam = mcbam.Bam(bam)
        pileups.append(currBam.get_pileup(myGene.chrom, myGene.start, myGene.end))

    # Average Pileups together
    avgPileup = mcbam.avg_pileups(pileups, fudgeFactor=args.ff)
    
    ################################################################################
    # MAKE WIGGLES
    ################################################################################

    mcwiggle.plot_wiggle(avgPileup, myModel, args.oname)

