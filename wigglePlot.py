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

def getOptions():
    """ Function to pull in arguments """
    description = """"This script starts with a begining pathway in the format
    of a 'path' file. It then builds SAS code for SEMs after adding a new gene
    to all possible locations within a network.

        (a) Downstream of endogenous genes
        (b) Downstream of exogenous genes
        (c) Upstream of endogenous genes
        (d) Upstream of exogenous genes
        (e) In between exogenous genes and their endogenous targets
        (e) In betwee endogenous genes and their endogenous targets
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", dest="pname", action='store', required=True, help="Name of the 'path' file [Required]")
    parser.add_argument("-l", dest="lname", action='store', required=True, help="Path to sas library that has the SAS dataset that will be analyzed [Required]")
    parser.add_argument("-m", dest="mname", action='store', required=True, help="Name of SAS dataset that will be analyzed [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name with PATH, of the output file ending with '.sas'; NOTE: model number will be appended to the filename [Required]")
    parser.add_argument("-n", dest="newGene", nargs='+', action='store', required=True, help="List of new gene isoforms separted by spaces [Required]")
    parser.add_argument("-g", dest="gname", action='store', required=True, help="The gene name of the new gene/isoforms being added [Required]")
    parser.add_argument("-t", dest="template", action='store', required=False, help="Name of the PROC CALIS template file [Optional]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 

    args = parser.parse_args()
    return(args)

if __name__ == '__main__':
    #args = getOptions()

    # Turn on logging
    mclib.logger.set_logger(args.log)
    try:
        # Output git commit version to log, if user has access
        mclib.git.git_to_log()
    except:
        pass

    # Grab gene annotation 
    ## Import GFF database
    gffName = '/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff'
    myGffDb = mcgff.FlyGff(gffName)

    ## Pull gene from database
    myGene = mcgff.FlyGene('InR', myGffDb)

    ################################################################################

    # Pull in pileup 
    bamName1 = '/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V1.sorted.bam'
    bamName2 = '/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V2.sorted.bam'

    b1 = mcbam.Bam(bamName1, myGene.chrom, myGene.start, myGene.end, fudgeFactor=True)
    b2 = mcbam.Bam(bamName2, myGene.chrom, myGene.start, myGene.end, fudgeFactor=True)

    b = mcbam.avg_pileups([b1,b2], fudgeFactor=True)
    

print sum((Counter(dict(x.pileups)) for x in [b1, b2]), Counter())



