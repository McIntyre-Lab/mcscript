#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import csv

csv.excel_tab

# Add-on packages

# McLab packages
import mclib_Python as mclib
from mclib_Python import vcf2 as mcvcf
from mclib_Python import bed as mcbed


def getOptions():
    """ Function to pull in arguments """
    description = """ Make a table of the number of snps present in a bed region. """
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    group1 = parser.add_argument_group(description="Input Files")
    group1.add_argument("--vcf", dest="vcfName", action='store', required=True, help="Name of VCF file. If not zipped using bgzip will try to create zip. [Required]")
    group1.add_argument("--bed", dest="bedName", action='store', required=True, help="Name of the 4-column bed file that contains the genomic regions you want to count SNPs in. [Required]")

    group2 = parser.add_argument_group(description="Output Files")
    group2.add_argument("--out", dest="oname", action='store', required=True, help="Name of output table in CSV format. [Required]")

    group2.add_argument("--log", dest="log", action='store', required=False, help="Name of the LOG file [Optional]")
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.")
    args = parser.parse_args()
#     args = parser.parse_args(['--vcf', '/home/jfear/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/r228_w11182r228_UPD.vcf.gz',
#                               '--bed', '/home/jfear/mclab/useful_dmel_data/flybase551/si_fusions/fb551_si_fusions.bed',
#                               '--out', '/home/jfear/tmp/test.csv',
#                               '--debug'])
    return(args)


def main(args):
    """ Main Script """
    # Import Files
    vcf = mcvcf.Vcf(args.vcfName)
    bed = mcbed.Bed(args.bedName)

    # Open output
    OUT = open(args.oname, 'w')
    OUT.write('regionID,num_snps,num_indels,total\n')

    # Iterate over bed regions and pull SNPs
    for chrom, start, end, region in bed:
        try:
            poly = vcf.pull_vcf_region(chrom, start, end)
            snps = 0
            indels = 0
            for record in poly:
                if record.is_snp:
                    snps += 1
                elif record.is_indel:
                    indels += 1

            total = snps + indels
            OUT.write(','.join(str(x) for x in [region, snps, indels, total]) + '\n')
        except:
            OUT.write('{},0,0,0\n'.format(region))

    OUT.close()


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
