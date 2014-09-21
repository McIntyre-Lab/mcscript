#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys
import re

# Add-on packages
import numpy as np
import matplotlib 
matplotlib.use('Agg')

# McLab Packages
import mclib
from mclib import bam as mcbam
from mclib import gff as mcgff
from mclib import vcf2 as mcvcf
from mclib import wiggle as mcwiggle

# TODO: Add ability to make overlapping wiggle when up to XX groups are given

def getOptions():
    """ Function to pull in arguments """

    description = """ This script constructs wiggle plots and gene models.
    Wiggle plots are created from sorted BAM files. Gene models require a GFF
    file at this time.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    group1 = parser.add_argument_group('Data Types', 'The following types of data are capable of being handled.')
    group1.add_argument("--bam", nargs='*', dest="bamList", action='store', required=True, help="Name with PATH to sorted BAM file. If multiple files are given they will be averaged. [Required]")
    group1.add_argument("--gff", dest="gffName", action='store', required=False, help="Name with PATH to a GFF file. Note if a GFFutils database has not been created it will be, this could take around 30min [Optional]")
    #group1.add_argument("--bed", dest="bedName", action='store', required=False, help="Name with PATH to a bed file [Optional].") # TODO: Add bed functionality
    group1.add_argument("--vcf", dest="vcfName", action='store', required=False, help="Name with PATH to a vcf file zipped using bgzip. [Optional]")
    #group1.add_argument("--fusions", dest="fusName", action='store', required=False, help="Name with PATH to a bed file contianing the genome coordinates to fusion name [Optional].") # TODO: Add fusion functionality

    group2 = parser.add_argument_group('Region of Interest', 'Select which region you want to focus on, one of the following is Required:')
    group2a = group2.add_mutually_exclusive_group(required=True)
    group2a.add_argument("--gene", dest="geneName", action='store', help="Name of the gene you want to make wiggles of.")
    group2a.add_argument("--region", dest="region", action='store', help="Name of the region of interset in the format 'chrom:start-end'")

    group3 = parser.add_argument_group('Factors', 'Various options to tweak output')
    group3.add_argument("--sample", dest="sample", action='store', required=False, help="Name of the sample you are working on. Note if you are using VCF file this name needs to match. [Optional]")
    group3.add_argument("--fudge", dest="ff", action='store', type=int, default=0, required=False, help="Use a fudge factor, if 999 then use the built in fudge factor calculation. [Optional]")
    group3.add_argument("--debug", dest="debug", action='store_true', required=False, help="Trun on debug output. [Optional]")

    group4 = parser.add_argument_group('Output')
    group4.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output PNG. [Required]")
    group4.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 

    args = parser.parse_args()
    #args = parser.parse_args(['--gff', '/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff', '--bam','/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V1.sorted.bam','/mnt/storage/cegs_aln/bam_fb551_genome_nodup/r101_V2.sorted.bam', '-g', 'InR', '-o', '/home/jfear/tmp/inr.png'])
    return(args)

if __name__ == '__main__':
    args = getOptions()

    # Turn on logging
    if args.debug:
        mclib.logger.set_logger(args.log, 'debug')
    else:
        mclib.logger.set_logger(args.log)

    # Output git commit version to log, if user has access
    mclib.git.git_to_log(__file__)

    ################################################################################
    # GENE ANNOTATION 
    ################################################################################

    if args.gffName and args.geneName:
        logging.info('Getting gene annotation from GFF file')
        ## Import GFF database
        myGffDb = mcgff.FlyGff(args.gffName)

        ## Pull gene from database
        myGene = mcgff.FlyGene(args.geneName, myGffDb)

        chrom = myGene.chrom
        start = myGene.start
        end = myGene.end

        ## Create Gene Model
        geneModel = mcwiggle.GeneModel(myGene)

    elif args.region:
        location = re.split(':|-',args.region)
        chrom = location[0]
        start = int(location[1])
        end = int(location[2])

        # TODO: would be nice to plot all gene models in a region
        geneModel = None

    else:
        logging.error('You need to specify either a gene name along with a GFF or a region to plot')
        raise ValueError

    ################################################################################
    # GENE COVERAGE
    ################################################################################

    logging.info('Creating pileups')

    # Pull in bam file and make gene pileup
    pileups = []
    for bam in args.bamList:
        currBam = mcbam.Bam(bam)
        pileups.append(currBam.get_pileup(chrom, start, end))

    # Average Pileups together
    avgPileup = mcbam.avg_pileups(pileups, fudgeFactor=args.ff)

    ################################################################################
    # Pull Variants if requested
    ################################################################################
    if args.vcfName:
        logging.info('Processing VCF file')

        # Attach vcf file
        myVcf = mcvcf.Vcf(args.vcfName)

        # Pull the region of interest
        region = myVcf.pull_vcf_region(chrom, start, end)

        # Grab inidividuals that are homozygous for an alternate base
        homz = myVcf.pull_homz(region=region, snp=True)

        # Create list of positions that have a variants. If a certain sample is
        # given only output variants for that sample.
        variantPos = list()
        for pos in homz:
            if args.sample:
                for line in homz[pos]:
                    if line == args.sample:
                        variantPos.append(pos)
            else:
                variantPos.append(pos)
    else:
        logging.debug('Setting VCF to None')
        variantPos = None

    ################################################################################
    # Make fusion models if requested
    ################################################################################

    # TODO: add fusion model 
    fusionModel=None

    ################################################################################
    # MAKE WIGGLES
    ################################################################################

    mcwiggle.plot_wiggle(avgPileup, args.oname, chrom, start, end, geneModel=geneModel, fusionModel=fusionModel, variantPos=variantPos)

