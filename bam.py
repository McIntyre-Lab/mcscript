import os.path
from collections import defaultdict

import numpy as np
import pysam

class Bam(object):
    def __init__(self, filename):
        """ Basic wrapper for pysam module, to allow the quick pulling of pileups
        for a certain region 

        Arguments:
        ----------
        filename = The full path to a sorted BAM file

        Attributes:
        -----------
        bam = a pysam Samfile object
        """
        # If there is not BAM index make it
        if not os.path.exists(filename + '.bai'):
            pysam.index(filename)
        
        # Import BAM and pull pileups for a gene
        self.bam = pysam.Samfile(filename, 'rb')

    def get_pile_dict(self, chrom, start, end):
        """ Create a dictionary where key is pos and value is count.

        Arguments:
        ----------
        chrom (str) = chromosome name.
        start (int) = start location
        end (int) = end location

        Returns:
        --------
        A dictionary where the key is the base position and the value is the
        base count.
        """
        pileups = self.bam.pileup(chrom, start, end)
        pileDict = {}
        for base in pileups:
            pileDict[base.pos]= base.n
        return pileDict

    def get_gene_read_count(self, chrom, start, end):
        """ Get the number of reads that aligned to a particular regions.
        Arguments:
        ----------
        chrom (str) = chromosome name.
        start (int) = start location
        end (int) = end location

        Returns:
        --------
        An (int) with the number of reads overlapping the region.
        """
        count = 0
        for alnRead in self.bam.fetch(chrom, start, end):
            count += 1
        return count

    def calc_base_level_fudgeFactor(self, normFactor=10**9):
        """ When comparing wiggle accross treatments, it is necessary to
        normalize reads. If you look at 
        
        RPKM = (# reads * 10^9) / (# total mapped reads * exon length)

        In contrast to RPKM, I want to calculate a fudge factor to apply to
        each base so I will the following:
            (1) Use number of mapped bases instead of reads
            (2) Remove the exon length

        Arguments:
        ----------
        normFactor = this is related to the 10^9 factor from RPKM

        Returns:
        --------
        A float with the global fudge factor 10^9 / sum(per base coverage)
        """
        base_count = 0
        for read in self.bam.fetch():
            base_count += read.inferred_length

        # Calculate the fudge factor similar to RPKM, but at the base level
        # instead of read level
        ff = float(normFactor) / float(base_count)
        return ff
    
def sum_pileups(pileups):
    """ Sum a list of pileup dictionaries

    Arguments:
    ----------
    pileups (list) = List of pileup dictionaries

    Returns:
    --------
    A new dcitionary containing the sum for each base
    """
    combinedDict = defaultdict(int)
    for pileup in pileups:
        for pos in pileup:
            combinedDict[pos] += pileup[pos]
    return combinedDict

def avg_pileups(pileups, fudgeFactor=False):
    """ Average a list of pileup dictionaries

    Arguments:
    ----------
    pileups (list) = List of pileup dictionaries
    fudgeFactor (bool) = do you want normalization or not

    Returns:
    --------
    A new dcitionary containing the average coverage for each base
    """
    num = len(pileups)

    if fudgeFactor:
        # Calculate the average fudge factor
        fList = [pileup.calc_base_level_fudgeFactor() for pileup in pileups]
        ff = np.mean(np.array(fList))

    combinedDict = sum_pileups(bamList)

    for pos in combinedDict:
        try:
            if fudgeFactor:
                combinedDict[pos] = combinedDict[pos] / num * ff
            else:
                combinedDict[pos] = combinedDict[pos] / num
        except:
            if combinedDict[pos] == 0:
                pass
            else:
                print 'There is something wrong with your pileup'
                raise ValueError
    return combinedDict
