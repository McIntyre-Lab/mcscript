import os.path
import collections

import numpy as np
import pysam

class Bam(object):
    """ Basic wrapper for pysam module, to allow the quick pulling of pileups
    for a certain region 
    """
    def __init__(self, filename, chrom, start, end, fudgeFactor=False):
        """ Arguments:
        filename = The full path to a sorted BAM file
        chrom = The gene of interests chromosome
        start = The gene of interests start location
        end = The gene of interests end location

        Attributes:
        bam = a pysam Samfile object
        pileups = dictionary of pileups where key is position and value is count
        fudgeFactor = a fudgeFactor related to RPKM that can be used to normalize base counts
        """

        # If there is not BAM index make it
        if not os.path.exists(filename + '.bai'):
            pysam.index(filename)
        
        # Import BAM and pull pileups for a gene
        self.bam = pysam.Samfile(filename, 'rb')
        self.pileups = self._get_pile_dict(chrom, start, end)

        # If needed calculate the normaliztion fudge factor
        if fudgeFactor:
            self.fudgeFactor = self._calc_base_level_fudgeFactor()

    def _get_pile_dict(self, chrom, start, end):
        """ Create a dictionary where key is pos and value is count """
        pileups = self.bam.pileup(chrom, start, end)
        pileDict = {}
        for base in pileups:
            pileDict[base.pos]= base.n
        return pileDict

    def _calc_base_level_fudgeFactor(self, normFactor=10**9):
        """ When comparing wiggle accross treatments, it is necessary to
        normalize reads. If you look at 
        
        RPKM = (# reads * 10^9) / (# total mapped reads * exon length)

        In contrast to RPKM, I want to calculate a fudge factor to apply to
        each base so I will the following:
            (1) Use number of mapped bases instead of reads
            (2) Remove the exon length

        Arguments:
        normFactor = this is related to the 10^9 factor from RPKM

        Returns a global fudge factor
        """
        base_count = 0
        for read in self.bam.fetch():
            base_count += read.inferred_length

        # Calculate the fudge factor similar to RPKM, but at the base level
        # instead of read level
        ff = float(normFactor) / float(base_count)
        return ff
    
    def get_pile(self, chrom, start, end):
        """ Create two lists first with positions and second with counts """
        pileups = self.bam.pileup(chrom, start, end)
        count = []
        for base in pileups:
            pos.append(base.pos)
            count.append(base.n)
        return pos, count
    
    def get_gene_read_count(self, chrom, start, end):
        """ Get the number of reads that aligned to a particular regions """
        count = 0
        for alnRead in self.bam.fetch(chrom, start, end):
            count += 1
        return count

def sum_pileups(bamList):
    """ Given a list of pileup dictionaries, create a new dictionary that is
    the sum accross positions 
    """
    combinedDict = collections.Counter()
    for bam in bamList:
        combinedDict.update(bam.pileups)
    return combinedDict

def avg_pileups(bamList, fudgeFactor=False):
    """ Given a list of pileup dictionaries, create a new dictionary that is
    the sum accross positions 
    Arguments:
    pileupList (list) = a list containing mclib.bam.Bam objects
    fudgeFactor (bool) = do you want normalization or not
    """
    num = len(bamList)

    if fudgeFactor:
        # Calculate the average fudge factor
        fList = [bam.fudgeFactor for bam in bamList]
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
