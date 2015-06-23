#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging

# Add-on packages

# McLab Packages
import mclib_Python as mclib
from mclib_Python import gff as mcgff

# Initialize global FUSION_ID Counter
COUNT = 1


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
    args = parser.parse_args()
    return(args)


def writeOutput(chrom, fusions, sfx, OUTbed, OUTtable):
    """ This functions names the fusions and writes their output. If a fusion
    is made from a single exon then it is a singleton and will be prefixed with
    the letter 'S'. If a fusion is made up of overlapping fusions then it is
    prefixed with a 'F'.

    Arguments:
    ----------
    chrom (str) = the current chromosome id
    fusions (generator) = is a generator of merged exons
    sfx (str) = the suffix to append on to the fusion id {'_SI', '_SD'}
    OUTbed (obj) = File output object for the BED file
    OUTtable (obj) = File output object for the Table file
    """

    # Attach the global counter
    global COUNT

    for fusion in fusions:
        # Get list of exons in a fusion
        exons = fusion['exonId']

        # Name the Fusion based on if it is a singleton or fusion.
        if fusion['merged']:
            # fusion
            name = "F{0}{1}".format(COUNT, sfx)

            # Are there multiple genes in this fusion
            ## Get a list of genes
            genes = set([':'.join(x.split(':')[:-1]) for x in exons])

            ## If there are more than 1 gene set flag_multigene
            if len(genes) == 1:
                flag_multigene = '0'
            else:
                flag_multigene = '1'
        else:
            # singleton
            name = "S{0}{1}".format(COUNT, sfx)

            # singletons by definition don't have multiple genes
            flag_multigene = '0'

        # Write output in a bed format
        start = str(fusion['start'] - 1)  # remember bed is a 0-based format and gff is 1-based
        end = str(fusion['end'])
        strand = str(fusion['strand'])
        myOut = '\t'.join([chrom, start, end, name, '.', strand]) + "\n"
        OUTbed.write(myOut)

        # Write output table linking fusion_id to exon_ID
        for exon in exons:
            gene = ':'.join(exon.split(':')[:-1])
            myout = '\t'.join([name, exon, gene, flag_multigene]) + "\n"
            OUTtable.write(myout)

        # increment fusion counter
        COUNT += 1

    logger.debug('%d fusions' % COUNT)


def main(args):
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
                logger.info('Pulling list of exons for chromosome: {0} on strand {1}'.format(chrom.id, currStrand))
                exons = flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end), strand = currStrand)

                # Create a list of fusions by merging overlapping exons
                logger.info('Merging overlapping exons on strand ' + currStrand)
                fusions = flyGff.merge(exons, ignore_strand=False)
                writeOutput(chrom.id, fusions, '_SD', OUTbed, OUTtable)
        else:
            # Create Strand independent fusions
            # Get list of all genes in the genome
            logger.info('Pulling list of exons for chromosome: ' + chrom.id)
            exons = flyGff.db.features_of_type('exon', limit=(chrom.id, chrom.start, chrom.end))

            # Create a list of fusions by merging overlapping exons
            logger.info('Merging overlapping exons')
            fusions = flyGff.merge(exons, ignore_strand=True)
            writeOutput(chrom.id, fusions, '_SI', OUTbed, OUTtable)

    logger.debug('%d total fusions' % COUNT)
    OUTbed.close()
    OUTtable.close()

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

    # Run Main part of the script
    main(args)
    logger.info("Script complete.")
