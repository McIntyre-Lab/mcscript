#!/usr/bin/env python

# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys
import re
from collections import defaultdict

# Add-on packages
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import numpy as np

# McLab Packages
import mclib
from mclib import vcf2 as mcvcf
from mclib import bed as mcbed

def getOptions():
    """ Function to pull in arguments """

    description = """ This script takes a VCF file and updates a FASTA file to incorporate variant calls. """
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    group1 = parser.add_argument_group(description="Input Files")
    group1.add_argument("--vcf", dest="vcfName", action='store', required=True, help="Name of VCF file. If not zipped using bgzip will try to create zip. [Required]")
    group1.add_argument("--fasta", dest="fastaName", action='store', required=True, help="Name of fasta file. [Required]")
    group1.add_argument("--bed", dest="bed", action='store', required=False, help="Name of a 4-column BED file, locations in BED file will be sliced out of FASTA file. [Optional]")

    group2 = parser.add_argument_group(description="Output Files")
    group2.add_argument("-o", dest="oname", action='store', required=True, help="Name of output FASTA. [Required]")
    group2.add_argument("--log", dest="log", action='store', required=False, help="Name of the LOG file [Optional]") 

    group3 = parser.add_argument_group(description="SNP options")
    group3.add_argument("--mask", dest="mask", action='store_true', required=False, help="Mask SNPs with 'N' instead of updating reference. [Optional]")
    group3.add_argument("--snps_only", dest="snpsOnly", action='store_true', required=False, help="Only update SNPs, ignore indels. Indels take significantly longer to run. [Optional]")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 

    args = parser.parse_args()
    return(args)

def buildCoordIndex(seqRecord):
    """ Create a numpy array where the index is the original coordinate and the
    value will be the updated coordinate. 
    
    Arguments:
    ----------
    seqRecord (Bio.SeqIO.Seq) = Bio-python SeqIO seq object containing the sequence
                                for the current chromosome

    Returns:
    --------
    coordIndex (numpy array) = array where the index is the original coordinate
                               and the value will be updated to the new coordinate.
    """
    bases = len(seqRecord.seq)
    coordIndex = np.array(xrange(0,bases))
    return coordIndex

def buildVariantDict(myVcf, snpsOnly):
    """ Build a dictionary where the key is chromosome and the value is a list
    of variants for that chromosomes.

    Arguments:
    ----------
    myVcf (mclib.vcf2.Vcf object) = an object created using the mclib vcf2 library
    
    Returns:
    --------
    variants (dict) = dictionary of variants where key is chromosome and value
                      is a list of variants organized in a list [position,
                      reference allele, alternative allele, difference between ref - alt]
    """
    # Iterate through each record and add to storage dictionary (variants)
    variants = defaultdict(list)
    for record in myVcf.vcf_reader:
        # Make sure that this VCF only contains a single sample by checking for commas
        if ',' in record.REF or ',' in record.ALT:
            logging.error('VCF file needs to be split by sample for this script to work')
            raise ValueError

        start = record.POS - 1      # VCF files are 1-based, convert to 0-based
        refbase = record.REF
        altbase = str(record.ALT[0])
        diff = len(altbase) - len(refbase)
        end = start + diff          # End is not really necessary until we start slicing bed files
        if snpsOnly:
            # If snpsOnly is on, only add SNPs, ie when diff = 0
            if diff == 0:
                variants[record.CHROM].append([start, end, refbase, altbase, diff])
        else:
            # Else add all variants
            variants[record.CHROM].append([start, end, refbase, altbase, diff])

    return variants

def adjustCoords(varList, coordList):
    """ Adjust the variant coordinates from a given location due to indels
    prior to the current variant. 

    Arguments:
    ----------
    varList (list) = List of variants with (position, reference allele,
                     alternative allele, difference between ref-alt) built from
                     buildVariantDict function

    coordIndex (numpy array) = array where the index is the original coordinate
                               and the value will be updated to the new coordinate.
    Returns:
    --------
    Updates the position in-place.
    """
    for record in varList:
        start = record[0]
        delta = record[4]
        if delta != 0:
            # Don't adjust coordinates for SNPs
            coordList[start+1:] = coordList[start+1:] + delta

def updateSeq(Seq, varList, coordList, chrom):
    """ Update the genomic sequence given a list of variants
    Arguments:
    ----------
    Seq (Bio.SeqIO.Seq) = Bio-python SeqIO seq object containing the sequence
                          for the current chromosome

    varList (list) = List of variants with (position, reference allele,
                     alternative allele, difference between ref-alt) built from
                     buildVariantDict function. Coordinates have been updated
                     by the adjustVarCoords function.

    coordIndex (numpy array) = array where the index is the original coordinate
                               and the value will be updated to the new coordinate.

    chrom (str) = Current chromosome, only used to debug output

    Returns:
    --------
    Updates the Seq object in-place, by adding variants to the sequence.
    """
    # Create a mutable sequence
    mut = Seq.seq.tomutable()

    # Iterate through variants and update
    for var in varList:
        origStart, origEnd, ref, alt, diff = var
        newStart = coordList[origStart]
        newEnd = coordList[origEnd]

        if diff == 0:
            # If a SNP
            if mut[newStart] == ref:
                if args.mask:
                    mut[newStart] = 'N'
                else:
                    mut[newStart] = alt
            else:
                logging.error('coordinates appear to be off for a SNP')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart], ref, chrom, origStart, newStart))
                raise ValueError
        elif diff > 0:
            # If a Insertion
            if mut[newStart] == ref[0]:
                cnt = newStart + 1
                for base in alt[1:]:
                    mut.insert(cnt, base)
                    cnt += 1
            else:
                logging.error('coordinates appear to be off for a Insertion')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart], ref, chrom, origStart, newStart))
                raise ValueError
        elif diff < 0:
            # If a Deletion
            if mut[newStart:newStart + len(ref)] == ref:
                mut[newStart:newStart + len(ref)] = alt
            else:
                logging.error('coordinates appear to be off for a deletion')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart:newStart + len(ref)], ref, chrom, origStart, newStart))
                raise ValueError
    Seq.seq = mut.toseq()

def sliceAndDiceSeq(bedRow, seqRecord):
    """ Slice out fusions from the genome
    Arguments:
    ----------
    bedRow (mcbed.Bed.BedRow) = An updated row from a bed file

    seqRecord (Bio.SeqIO.Seq) = Bio-python SeqIO seq object containing the sequence
                                for the current chromosome

    Returns:
    --------
    fusRecord (Bio.SeqIO.Seq) = A new Bio-python SeqIO seq object containing the sequence
                                for the current fusion
    """
    fusID = bedRow['name']
    fusSeq = seqRecord[bedRow['chromStart']:bedRow['chromEnd']].seq
    fusRecord = SeqRecord(fusSeq, id=fusID, description='')
    return fusID, fusRecord

def updateBed(coordIndex, chrom, mySeq, myBed, fusions):
    """ Take a bed file and update it coordinate and then slice the genomic
    region.

    Arguments:
    ----------
    coordIndex (numpy array) = array where the index is the original coordinate
                               and the value will be updated to the new coordinate.

    chrom (str) = Current chromosome, only used to debug output

    mySeq (Bio.SeqIO) = A dictionary where key is chromosome and the value is a
                        SeqRecord for that chromosome

    myBed (mclab.bed.Bed) = A reader for a Bed file

    fusions (dict) = A dictionary where key is a fusion id and the value is a
                     fusions SeqRecord
    
    Returns:
    --------
    Updates fusions in place.
    """
    try:
        for row in myBed.get_rows(name=chrom):
            start = row['chromStart']
            end = row['chromEnd']
            newStart = coordIndex[start]
            newEnd = coordIndex[end]

            row['chromStart'] = newStart
            row['chromEnd'] = newEnd
            fusID, fusRecord = sliceAndDiceSeq(row, mySeq[chrom])
            fusions[fusID] = fusRecord
    except:
        logging.warn('The chromosome: {0} did not have any fusions associated with it.'.format(chrom))

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
    # Import VCF information
    ################################################################################
    logging.info('Importing VCF information: %s' % args.vcfName)
    myVcf = mcvcf.Vcf(args.vcfName)
    variants = buildVariantDict(myVcf, args.snpsOnly)

    if args.snpsOnly:
        logging.INFO('You are running in SNPONLY mode, remove --snps_only flag to include indels')

    ################################################################################
    # Import FASTA
    ################################################################################
    logging.info('Importing FASTA: %s' % args.fastaName)
    mySeq = SeqIO.to_dict(SeqIO.parse(open(args.fastaName, 'r'), 'fasta'))

    ################################################################################
    # Import Bed File
    ################################################################################
    fusions = dict()
    if args.bed:
        logging.info('Importing BED: %s' % args.bed)
        myBed = mcbed.Bed(args.bed)
    else:
        pass

    ################################################################################
    # Iterate through the chromosomes and update the genome
    ################################################################################
    logging.info('Identifying variants and updating genome')
    for chrom in mySeq:
        logging.info('{0}: Building coordinate Index'.format(chrom))
        coordIndex = buildCoordIndex(mySeq[chrom])

        if not args.snpsOnly:
            logging.info('{0}: Adjusting coordinates'.format(chrom))
            adjustCoords(variants[chrom], coordIndex)

        logging.info('{0}: Updating sequences'.format(chrom))
        updateSeq(mySeq[chrom], variants[chrom], coordIndex, chrom)

        if args.bed:
            # If a BED file was provided slice out the coordinates from updated
            # sequence
            logging.info('{0}: Updating BED coordinates'.format(chrom))
            updateBed(coordIndex, chrom, mySeq, myBed, fusions)
    
    ################################################################################
    # Output Updated FASTA file
    ################################################################################
    if fusions:
        # If there are fusions (i.e. if a bed file) then output fusions.
        logging.info('Outputing updated fusions')
        myOut = fusions
    else:
        logging.info('Outputing updated genome')
        myOut = mySeq

    with open(args.oname, 'w') as OUT:
        for record in myOut.values():
            SeqIO.write(record, OUT, "fasta")
