import os.path
import subprocess
import logging
import collections

import vcf as pyvcf

class Vcf(object):
    def __init__(self, filename, build=True):

        # Determine if the input file was compressed by bgzip and has a tabix
        # index
        ftest = subprocess.check_output(["file", filename])
        if "gzip compressed data, extra field" not in ftest:
            if os.path.isfile(filename + '.gz'):
                filename = filename + '.gz'
            else:
                logger.warn('Input file needs to be compressed by bgzip')
                logger.info('I will try to run bgzip now')

                bgz = filename+'.gz'
                with open(bgz, 'w') as FH:
                    p = subprocess.call(['bgzip', '-c', filename], stdout=FH) 

                if p == 0:
                    logger.info('Sucessfully made bzipped file')
                    filename = bgz
                else:
                    logger.error('I could not compress the file using bgzip, please run "bgzip {0}"'.format(filename))
                    raise IOError

        if not os.path.isfile(filename + '.tbi'):
            logger.warn('A Tabix index file was not found.')
            logger.info('Trying to create tabix index file')
            try:
                subprocess.check_call(['tabix', filename, '-p', 'vcf'])
                logger.info('Sucessfully built tabix index')
            except:
                logger.error('I could not create the index file, please run "tabix {0} -p vcf"'.format(filename))

        # Create pyvcf object
        self.vcf_reader = pyvcf.Reader(filename=filename, compressed=True)

    def pull_vcf_region(self, chrom, start, end):
        """ PUll all records that belong to a spcific region 

        Arguments:
        ----------
        chrom (str) = chromosome
        start (int) = start location
        end (int) = end location

        Returns:
        --------
        A pyVCF object that contains all of the records for the given region.
        """
        return self.vcf_reader.fetch(chrom, start, end)

    def pull_homz(self, region=False, snp=True, indel=False, sv=False):
        """ Pull homozygous individuals that have a differ from ref base. This
        difference can be due to a SNP, InDel or structural variant.

        Arguments:
        ----------
        region (pyVCF) = pyVCF object that contains a list of records. If not
        object is passed will use the full pyVCF file stored by the Vcf class
        snp (bool) = True if you want snps to be included
        indel (bool) = True if you want indels to be included
        sv (bool) = True if you want structural variants to be included

        """
        if not region:
            region = self.vcf_reader

        variants = collections.defaultdict(list)
        for record in region:
            if snp and record.is_snp:
                for homz in record.get_hom_alts():
                    variants[record.POS].append(homz.sample)
            if indel and record.is_indel:
                for homz in record.get_hom_alts():
                    variants[record.POS].append(homz.sample)
            if sv and record.is_sv:
                for homz in record.get_hom_alts():
                    variants[record.POS].append(homz.sample)
        return variants
