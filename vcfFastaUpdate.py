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

# McLab Packages
import mclib
from mclib import vcf2 as mcvcf
from mclib import bed as mcbed

def getOptions():
    """ Function to pull in arguments """
    description = """ This script updates FASTA files. """
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", dest="vcfName", action='store', required=True, help="Name with PATH to a vcf file zipped using bgzip. [Required]")
    parser.add_argument("--fasta", dest="fastaName", action='store', required=True, help="Name with PATH to a fasta file. [Required]")
    parser.add_argument("--mask", dest="mask", action='store_true', required=False, help="Select if you just want to mask output with an 'N'. [Optional]")
    parser.add_argument("--bed", dest="bed", action='store', required=False, help="Name and PATH to a bed file containing fusions locations. [Optional]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name with PATH of the output FASTA. [Required]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    #args = parser.parse_args()
    args = parser.parse_args(['--vcf', '/home/jfear/tmp/indels/r101.vcf', '--fasta','/home/jfear/storage/useful_dmel_data/dmel-all-chromosome-r5.57.fasta', '-o', '/home/jfear/tmp/indel.fasta', '--bed', '/home/jfear/storage/useful_dmel_data/dmel_si_fusions_r5.57.bed', '--debug'])
    #args = parser.parse_args(['--vcf', '/home/jfear/tmp/test.vcf', '--fasta','/home/jfear/tmp/test.fasta', '-o', '/home/jfear/tmp/indel.fasta', '--debug'])
    return(args)

def buildVariantDict(myVcf):
    """ Build a dictionary of variants organized by chromosome 
    Arguments:
    ----------
    myVcf (mclib.vcf2.Vcf object) = an object created using the mclib vcf library
    
    Returns:
    --------
    variants (dict) = dictionary of variants where key is chromosome and value
                      is a list of variants orgznized in a list [position,
                      reference allele, alternative allele, difference between ref - alt]
    """
    variants = defaultdict(list)
    varIndex = defaultdict(dict)
    for record in myVcf.vcf_reader:
        # Make sure that this is a split VCF by checking for commas
        if ',' in record.REF or ',' in record.ALT:
            logging.error('VCF file needs to be split by sample for this scirpt to work')
            raise ValueError

        start = record.POS - 1
        refbase = record.REF
        altbase = str(record.ALT[0])
        diff = len(altbase) - len(refbase)
        end = start + diff
        variants[record.CHROM].append([start, end, refbase, altbase, diff])

    for chrom in variants:
        for index, value in enumerate(variants[chrom]):
            varIndex[chrom][value[0]] = index
    return variants, varIndex

def adjustVarCoords(varList):
    """ Function to adjust the variant coordinates from a given location 
    Arguments:
    ----------
    Var (list of variants organized in a list) = from buildVariantDict function
    index (int) = current index of the variant list
    adjust (int) = How to adjust the positions 
                    {0: SNP no adjustment, 
                    positve: Insertion -- number of position to add to positions, 
                    negative: Deletion -- number of positions to subract
    """
    for index, record in enumerate(varList):
        delta = record[4]
        if index + 1 < len(varList):
            for row in varList[index + 1:]:
                row[0] = row[0] + delta
                row[1] = row[1] + delta

def updateSeq(Seq, varList):
    """ Update the sequence given a list of variants """
    # Create a mutable sequence
    mut = Seq.seq.tomutable()
    for var in varList:
        start, end, ref, alt, diff = var
        if diff == 0:
            if mut[start] == ref:
                if args.mask:
                    mut[start] = 'N'
                else:
                    mut[start] = alt
            else:
                logging.error('coordinates appear to be off for a SNP')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Pos: {2}'.format(mut[start], ref, start))
                raise ValueError
        elif diff > 0:
            if mut[start] == ref:
                cnt = start + 1
                for base in alt[1:]:
                    mut.insert(cnt, base)
                    cnt += 1
            else:
                logging.error('coordinates appear to be off for a Insertion')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Pos: {2}'.format(mut[start], ref, start))
                raise ValueError
        elif diff < 0:
            if mut[start:start + len(ref)] == ref:
                mut[start:start + len(ref)] = alt
            else:
                logging.error('coordinates appear to be off for a deletion')
                logging.debug('Seq ref: {0}, VCF ref: {1}, Pos: {2}'.format(mut[start:start + len(ref)], ref, start))
                raise ValueError
    Seq.seq = mut.toseq()

def getVariants(chrom, start, end, variants, varIndex):
    origVarIndex = []
    for pos in xrange(start, end+1):
        try:
            origVarIndex.append(varIndex[chrom][pos])
        except:
            pass
    if origVarIndex:
        newPos = []
        for index in origVarIndex:
            start = variants[chrom][index][0]
            end = variants[chrom][index][1]
            newPos + [start, end]
        newStart = min(newPos)
        newEnd = max(newPos)
    else:
        newStart = start
        newEnd = end
    return newStart, newEnd

def sliceAndDiceSeq(bedRow, seqRecord):
    fusID = bedRow['name']
    fusSeq = seqRecord[bedRow['chromStart']:bedRow['chromEnd']].seq
    fusRecord = SeqRecord(fusSeq, id=fusID, description='')
    return fusID, fusRecord

def updateBedAndSlice(Bed, variants, varIndex, mySeq):
    bedReader = mcbed.Bed(Bed)
    updateSeq = dict()
    for row in bedReader.get_all_rows():
        chrom = row['chrom']
        start = row['chromStart']
        end = row['chromEnd']
        newStart, newEnd = getVariants(chrom, start, end, variants, varIndex)
        if start != newStart:
            logging.debug('Original Coords: {0}, New Coords: {1}'.format((start, end), (newStart, newEnd)))
        row['chromStart'] = newStart
        row['chromEnd'] = newEnd
        fusID, fusRecord = sliceAndDiceSeq(row, mySeq[chrom])
        updateSeq[fusID] = fusRecord
    return updateSeq

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
    logging.info('Importing VCF information')
    myVcf = mcvcf.Vcf(args.vcfName)
    variants, varIndex = buildVariantDict(myVcf)

    ################################################################################
    # Import FASTA
    ################################################################################
    logging.info('Importing FASTA')
    mySeq = SeqIO.to_dict(SeqIO.parse(open(args.fastaName, 'r'), 'fasta'))

    ################################################################################
    # Iterte through the chromosomes and update the genome
    ################################################################################
    logging.info('Identifying variants and update geneome')
    for chrom in mySeq:
        adjustVarCoords(variants[chrom])
        updateSeq(mySeq[chrom], variants[chrom])

    ################################################################################
    # Update bed coordinates
    ################################################################################
    if args.bed:
        fusions = updateBedAndSlice(args.bed, variants, varIndex, mySeq)
        myOut = fusions
    else:
        myOut = mySeq

    with open(args.oname, 'w') as OUT:
        for record in myOut.values():
            SeqIO.write(record, OUT, "fasta")
