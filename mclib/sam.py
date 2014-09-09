import pysam
import os.path

class Sam(object):
    def __init__(self, filename):

        if not os.path.exists(filename + '.bai'):
            pysam.index(filename)
        
        self.sam = pysam.Samfile(filename, 'rb')

    def get_pile(self, chrom, start, end):
        pileups = self.sam.pileup(chrom, start, end)
        pos = []
        count = []
        for base in pileups:
            pos.append(base.pos)
            count.append(base.n)
        return pos, count
    
    def get_read_count(self, chrom, start, end):
        count = 0
        for alnRead in self.sam.fetch(chrom, start, end):
            count += 1
        return count
