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
    description = """ This script combines overlapping exons into exonic
    regions (i.e., fusions). The script takes a GFF file, it then takes each
    chromosome and identifies overlapping exons. The script outputs a BED FILE
    with genomic coordinates and the associated fusion_id.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--gff", dest="gffName", action='store', required=True, help="Name with PATH, to a GFF file. Note if the GFFutils database has not been created it will be, this could take around 30min [Required]")
    parser.add_argument("--sd", dest='strand', action='store_true', required=False, help="Flag if you want to make strand dependent fusions. The default is to make strand independent fusions [Optional]")
    parser.add_argument("--obed", dest="obed", action='store', required=True, help="Name of the output BED. [Required]")
    parser.add_argument("--otable", dest="otable", action='store', required=True, help="Name of the output Table Relating Fusions to Exons. [Required]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    #args = parser.parse_args()
    args = parser.parse_args(['--gff', '/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.57.gff', '--obed', '/home/jfear/Desktop/fb557_si.bed', '--otable', '/home/jfear/Desktop/fb551_si.tsv','--debug'])
    return(args)

def writeOutput(chrom, fusions, strand, sfx, OUTbed, OUTtable):
    """ This functions names the fusions and writes their output. If a fusion
    is made from a single exon then it is a singleton and will be prefixed with
    the letter 'S'. If a fusion is made up of overlapping fusions then it is
    prefixed with a 'F'.

    Arguments:
    ----------
    chrom (str) = the current chromosome id
    fusions (generator) = is a generator of merged exons
    strand (str) = the strand to output; {'-', '+'} for SD and {'.'} for SI fusions
    sfx (str) = the suffix to append on to the fusion id {'_SI', '_SD'}
    OUTbed (obj) = File output object for the BED file
    OUTtable (obj) = File output object for the Table file

    """

    # Attach the global counter
    global cnt

    for fusion in fusions:
        # Get list of exons in a fusion
        exons = fusion['exonId']

        # Name the Fusion based on if it is a singleton or fusion.
        if fusion['merged']:
            # fusion
            name = "F{0}{1}".format(cnt, sfx)

            # Are there multiple genes in this fusion
            ## Get a list of genes
            genes = set([x.split(':')[0] for x in exons])

            ## If there are more than 1 gene set flag_multigene
            if len(genes) == 1:
                flag_multigene = '0'
            else:
                flag_multigene = '1'
        else:
            # singleton
            name = "S{0}{1}".format(cnt, sfx)

            # singletons by definition don't have multiple genes 
            flag_multigene = '0'

        # Write output in a bed format 
        start = str(fusion['start'] - 1) # remember bed is a 0-based format and gff is 1-based
        end = str(fusion['end'])
        myOut = '\t'.join([chrom, start, end, name, '.', strand]) + "\n"
        OUTbed.write(myOut)

        # Write output table linking fusion_id to exon_ID
        for exon in exons:
            gene = exon.split(':')[0]
            myout = '\t'.join([name, exon, gene, flag_multigene]) + "\n"
            OUTtable.write(myout)

        # increment fusion counter
        cnt +=1

if __name__ == '__main__':
    # Turn on Logging if option -g was given
    args = getOptions()

    # Turn on logging
    logger = logging.getLogger()
    if args.debug:
        mclib.logger.setLogger(logger, args.log, 'debug')
    else:
        mclib.logger.setLogger(logger, args.log)

    # Output git commit version to log, if user has access
    mclib.git.git_to_log(__file__)

    ################################################################################
    # Initialize Variables and Files
    ################################################################################

    # Connect to the database
    logger.info('connecting to the gff file')
    flyGff = mcgff.FlyGff(args.gffName)

    # Open output BED file
    OUTbed = open(args.obed, 'w')

    # Open output Table file and add header
    OUTtable = open(args.otable, 'w')
    myout = '\t'.join(['fusion_id', 'exon_id', 'primary_FBgn', 'flag_multigene']) + "\n"
    OUTtable.write(myout)

    # Initialize FUSION_ID Counter
    cnt = 1

    ################################################################################
    # Build Fusions
    ################################################################################

    # Get list of chromosomes
    chromosomes = flyGff.get_chrom()

    # Take each chromosome and identify overlapping exons
    for chrom in chromosomes:
        if args.strand:
            # Create Strand dependent fusions
            for currStrand in ('-', '+'):
                # Get list of all genes in the genome
                logger.info('Pulling list of exons for chromosome: {0} on strand {1}'.format(chrom.id,currStrand))
                exons = flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end),strand = currStrand)

                # Create a list of fusions by merging overlapping exons
                logger.info('Merging overlapping exons on strand '+currStrand)
                fusions = flyGff.merge(exons, ignore_strand=False)
                writeOutput(chrom.id, fusions, currStrand, '_SD', OUTbed, OUTtable)
        else:
            # Create Strand independent fusions
            # Get list of all genes in the genome
            logger.info('Pulling list of exons for chromosome: '+chrom.id)
            exons = flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end))

            # Create a list of fusions by merging overlapping exons
            logger.info('Merging overlapping exons')
            fusions = flyGff.merge(exons, ignore_strand=True)
            writeOutput(chrom.id, fusions, '.' , '_SI', OUTbed, OUTtable)

    OUTbed.close()
    OUTtable.close()
