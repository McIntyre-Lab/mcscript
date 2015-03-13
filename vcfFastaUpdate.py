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
    group3.add_argument("--force-mask", dest="fmask", action='store', required=False, help="Given a BED file of regions, these regions will be masked with 'N' before updating the reference. [Optional]")
    group3.add_argument("--snps-only", dest="snpsOnly", action='store_true', required=False, help="Only update SNPs, ignore indels. Indels take significantly longer to run. [Optional]")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    parser.add_argument("--debug-chrom", dest="dchrom", action='store', required=False, help="When debugging, select one chromosome to parse.") 
    parser.add_argument("--debug-fusion", dest="dfus", action='store', required=False, help="When debugging, output information related to this fusion.") 

    args = parser.parse_args()
    return(args)

def buildVariantDict(myVcf, snpsOnly):
    """ Build a dictionary where the key is chromosome and the value is a list
    of variants for that chromosomes.

    Args:
        myVcf (mclib.vcf2.Vcf object): an object created using the mclib vcf2
            library
    
    Returns:
        variants (dict): dictionary of variants where: 
            key is chromosome, values are a list of lists with: 
            [start, reference allele, alternative allele, 
            difference between alt - ref, length of the reference base]. 
            
            For example:
                {'2L': [[12, 'AAATTA', 'A', -5, 6], [32, 'A', 'T', 0, 1]}

    """
    # Iterate through each record and add to storage dictionary (variants)
    variants = defaultdict(list)
    for record in myVcf.vcf_reader:
        # Make sure that this VCF only contains a single sample by checking for commas
        if ',' in record.REF or ',' in record.ALT:
            logger.error('VCF file needs to be split by sample for this script to work')
            logger.debug(record)
            raise ValueError

        start = record.POS - 1      # VCF files are 1-based, convert to 0-based
        refbase = record.REF
        lref = len(refbase)
        altbase = str(record.ALT[0])
        lalt = len(altbase)
        diff = lalt - lref

        if snpsOnly:    # If snpsOnly is on, only add SNPs, ie when diff = 0
            if diff == 0:
                variants[record.CHROM].append([start, refbase, altbase, diff, lref])
        else:   # Else add all variants
            variants[record.CHROM].append([start, refbase, altbase, diff, lref])

    return variants

def debugVariants(dFus, variants):
    """ Print the variants located within a fusion 
    Args:
        dFus (obj): Is False when no fusion is passed to --debug-fusion.
            Otherwise, contains coordinate information about a given fusion.

        variants (dict): Dictionary of variants created by buildVariantDict function.

    """
    if dFus:
        positions = set(range(dFus.start, dFus.end))
        for variant in variants[dFus.chrom]:
            pos, ref, alt, diff, lref = variant

            if diff == 0:
                locs  = set([pos,])
            else:
                locs = set(range(pos, pos + lref))

            if positions & locs:
                logger.debug('Printing variant information for fusion {0}:{1}--{2}.\n{3}\n'.format(dFus.id, dFus.start, dFus.end,variant))

def force_masking(Seq, chrom, maskBed):
    """ Masks positions with N instead of changing them.
    Args:
        Seq (Bio.SeqIO.Seq): Bio-python SeqIO seq object containing the
            sequence for the current chromosome

        chrom (str): Current chromosome, only used to debug output

        maskBed (str): File name for a BED file containing the coordinates to mask

    Returns:
        Masks the Seq object in-place, using coordinates provided by maskBed.

    """
    # Create a mutable sequence
    mut = Seq.seq.tomutable()

    # Import the BED file for masking
    myBed = mcbed.Bed(maskBed)

    # Iterate over BED file and mask regsions with 'N'
    for row in myBed.get_all_rows(name=chrom):
        start = row['chromStart']
        end = row['chromEnd']
        mut = mut[start:end] = 'N'
    Seq.seq = mut.toseq()

def buildCoordIndex(seqRecord):
    """ Create chromosome coordinate array.
    
    Create a numpy array where the index is the original coordinate and the
    value will be the updated coordinate. 
    
    Args:
        seqRecord (Bio.SeqIO.Seq): Bio-python SeqIO seq object containing the
            sequence for the current chromosome

    Returns:
        coordIndex (numpy array): array where the index is the original
            coordinate and the value will be updated to the new coordinate.

        delMask (numpy array): array with a 0 for each position in the genome.
            0 will be updated to 1 if there was a deletion at that position.

    """
    bases = len(seqRecord.seq)
    coordIndex = np.arange(0,bases)  
    delMask = np.zeros(bases)
    return coordIndex, delMask

def debugCoords(dFus, coords, description):
    """ Output additional debugging information for coordinates array.

    When running debugging I want to be able to output additional information
    about a specific fusion.

    Args:
        dFus (obj): Is False when no fusion is passed to --debug-fusion.
            Otherwise, contains coordinate information about a given fusion.

        coords (list): A list of coordinates, where the index is the original
            coordinate and the value is the updated coordinate. The original
            list is from the buildCoordIndex function.

        description (str): A string describing what step we are on.

    Returns:
        Outputs current value and debug message to STDOUT.

    """
    if dFus:
        logger.debug(description + ' for fusion {0}:{1}--{2}.\n{3}\n'.format(dFus.id, dFus.start, dFus.end,[coords[dFus.start], coords[dFus.end]]))

def adjustCoords(variants, coordList, delMask):
    """ Adjust the variant coordinates. 
    
    To correct for coordinate changes due to indels, cycle through each
    variant and update coordinate locations for all downstream variants. 

    Args:
        variants (list): List of variants with (position, reference allele,
            alternative allele, difference between ref-alt) built from
            buildVariantDict function

        coordList (numpy array): array where the index is the original
            coordinate and the value will be updated to the new coordinate.

    Returns:
        Updates position in coordList (numpy array).

    """
    for variant in variants:
        start, ref, alt, diff, lref = variant
        if diff != 0:      # Don't adjust coordinates for SNPs
            coordList[start + lref:] = coordList[start + lref:] + diff
            if diff < 0:    # Mask deletions
                delMask[start + len(alt):start + lref] = 1

def debugMask(dFus, delMask):
    """ Output additional debugging information for mask array.

    When running debugging I want to be able to deletion masking from a
    specific fusion.

    Args:
        dFus (obj): Is False when no fusion is passed to --debug-fusion.
            Otherwise, contains coordinate information about a given fusion.

        delMask (list): A list of 0|1 for the entire length of the chromosome.
            A 1 represents that base was deleted.

    Returns:
        List of positions, where 1 indicates a deletion has taken place at that
            position. 

    """
    if dFus:
        logger.debug('Masked Array for fusion {0}:{1}--{2}.\n{3}\n'.format(dFus.id, dFus.start, dFus.end, delMask[dFus.start:dFus.end]))

def updateSeq(Seq, variants, coordList, chrom):
    """ Update the genomic sequence given a list of variants

    Args:
        Seq (Bio.SeqIO.Seq): Bio-python SeqIO seq object containing the
            sequence for the current chromosome

        variants (list): List of variants with (position, reference allele,
            alternative allele, difference between ref-alt) built from
            buildVariantDict function. Coordinates have been updated by the
            adjustVarCoords function.

        coordIndex (numpy array): array where the index is the original
            coordinate and the value will be updated to the new coordinate.

        chrom (str): Current chromosome, only used to debug output

    Returns:
        Updates the Seq object in-place, by adding variants to the sequence.

    """
    # Create a mutable sequence
    mut = Seq.seq.tomutable()

    # Iterate through variants and update
    for variant in variants:
        origStart, ref, alt, diff, lref = variant
        newStart = coordList[origStart]

        if diff == 0:       # If a SNP
            if mut[newStart] == ref:
                if args.mask:
                    mut[newStart] = 'N'
                else:
                    mut[newStart] = alt
            else:
                logger.error('coordinates appear to be off for a SNP')
                logger.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart], ref, chrom, origStart+1, newStart))
                raise ValueError
        elif diff > 0:      # If a Insertion
            if mut[newStart:newStart + lref] == ref:
                mut[newStart:newStart + lref] = alt
            else:
                logger.error('coordinates appear to be off for a Insertion')
                logger.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart], ref, chrom, origStart+1, newStart))
                raise ValueError
        elif diff < 0:      # If a Deletion
            if mut[newStart:newStart + lref] == ref:
                mut[newStart:newStart + lref] = alt
            else:
                logger.error('coordinates appear to be off for a deletion')
                logger.debug('Seq ref: {0}, VCF ref: {1}, Chrom: {2} Original Pos: {3} New Pos {4}'.format(mut[newStart:newStart + len(ref)], ref, chrom, origStart, newStart))
                raise ValueError

    Seq.seq = mut.toseq()

def sliceAndDiceSeq(bedRow, seqRecord):
    """ Slice out fusions from the genome

    Args:
        bedRow (mcbed.Bed.BedRow): An updated row from a bed file

        seqRecord (Bio.SeqIO.Seq): Bio-python SeqIO seq object containing the
            sequence for the current chromosome

    Returns:
        fusRecord (Bio.SeqIO.Seq): A new Bio-python SeqIO seq object
            containing the sequence for the current fusion

    """
    fusID = bedRow['name']
    fusSeq = seqRecord[bedRow['chromStart']:bedRow['chromEnd']].seq
    fusRecord = SeqRecord(fusSeq, id=fusID, description='')
    return fusID, fusRecord

def updateBed(coordIndex, delMask, chrom, mySeq, myBed, fusions):
    """ Take a bed file and update it coordinate and then slice the genomic
    region.

    Args:
        coordIndex (numpy array): array where the index is the original
            coordinate and the value will be updated to the new coordinate.

        delMask (list): A list of 0|1 for the entire length of the chromosome.
            A 1 represents that base was deleted.

        chrom (str): Current chromosome, only used to debug output

        mySeq (Bio.SeqIO): A dictionary where key is chromosome and the value
            is a SeqRecord for that chromosome

        myBed (mclab.bed.Bed): A reader for a Bed file

        fusions (dict): A dictionary where key is a fusion id and the value is
            a fusions SeqRecord
    
    Returns:
        Updates fusions in place.

    """
    try:
        for row in myBed.get_rows(name=chrom):
            start = row['chromStart']
            end = row['chromEnd']

            # Move fusion start position up if it was deleted
            upStart = start
            flagDel = True
            while flagDel:
                if delMask[upStart] == 0:
                    flagDel = False
                else:
                    upStart += 1

            if upStart != start:
                logger.debug('Incremented start position {0} to {1} because it was deleted.'.format(start, upStart))

            newStart = coordIndex[upStart]
            newEnd = coordIndex[end]

            row['chromStart'] = newStart
            row['chromEnd'] = newEnd
            fusID, fusRecord = sliceAndDiceSeq(row, mySeq[chrom])
            if len(fusRecord) > 0:
                fusions[fusID] = fusRecord
            else:
                logger.warn("The exonic region: {0} had a length of 0 [{1}-{2}] and is being ignored.".format(fusID,newStart, newEnd))
    except:
        logger.warn('The chromosome: {0} did not have any fusions associated with it.'.format(chrom))

def main(args):
    ################################################################################
    # Import Bed File
    ################################################################################
    fusions = dict()
    if args.bed:
        logger.info('Importing BED: %s' % args.bed)
        myBed = mcbed.Bed(args.bed)
    else:
        pass

    # When debugging with a test fusion, grab fusions information from BED file
    if args.debug and args.dfus and args.bed:
        dFus = myBed.get_rows(args.dfus)[0]
        dFus.chrom = dFus[0]
        dFus.start = dFus[1]
        dFus.end = dFus[2]
        dFus.id = dFus[3]
        args.dchrom = dFus.chrom
        logger.debug('Printing additional information for fusion: {0}'.format((dFus.id, dFus.chrom, dFus.start, dFus.end)))
    else:
        dFus = False

    ################################################################################
    # Import VCF information
    ################################################################################
    logger.info('Importing VCF information: %s' % args.vcfName)
    myVcf = mcvcf.Vcf(args.vcfName)
    variants = buildVariantDict(myVcf, args.snpsOnly)
    debugVariants(dFus, variants)

    if args.snpsOnly:
        logger.info('You are running in SNPONLY mode, remove --snps-only flag to include indels')

    ################################################################################
    # Import FASTA
    ################################################################################
    logger.info('Importing FASTA: %s' % args.fastaName)
    mySeq = SeqIO.to_dict(SeqIO.parse(open(args.fastaName, 'r'), 'fasta'))
    if args.fmask:
        logger.info('Force Masking Provided Regions: %s' % args.fmask)
        for chrom in mySeq:
            force_masking(mySeq[chrom], args.fmask)

    ################################################################################
    # Iterate through the chromosomes and update the genome
    ################################################################################
    logger.info('Identifying variants and updating genome')
    for chrom in mySeq:
        # If debugging and dchrom is given then only process the given chromosome
        if args.debug and args.dchrom and chrom != args.dchrom:
            logger.debug('Skipping chromosome {0}.'.format(chrom))
            continue

        logger.info('{0}: Building coordinate Index'.format(chrom))
        coordIndex, delMask = buildCoordIndex(mySeq[chrom])
        debugCoords(dFus, coordIndex, 'Initial coordinates')

        if not args.snpsOnly:
            logger.info('{0}: Adjusting coordinates'.format(chrom))
            adjustCoords(variants[chrom], coordIndex, delMask)
            debugCoords(dFus, coordIndex, 'Adjusted coordinates')
            debugMask(dFus, delMask)

        logger.info('{0}: Updating sequences'.format(chrom))
        updateSeq(mySeq[chrom], variants[chrom], coordIndex, chrom)

        # If a BED file was provided slice out the coordinates from updated
        # sequence
        if args.bed:
            logger.info('{0}: Updating BED coordinates'.format(chrom))
            updateBed(coordIndex, delMask, chrom, mySeq, myBed, fusions)
    
    ################################################################################
    # Output Updated FASTA file
    ################################################################################
    if fusions:
        # If there are fusions (i.e. if a bed file) then output fusions.
        logger.info('Outputing updated fusions')
        myOut = fusions
    else:
        logger.info('Outputing updated genome')
        myOut = mySeq

    with open(args.oname, 'w') as OUT:
        for record in myOut.values():
            SeqIO.write(record, OUT, "fasta")

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
