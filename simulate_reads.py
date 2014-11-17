#!/usr/bin/env python

# Built-in packages
import argparse
import logging

# Add-on packages
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

# McLab Packages
import mclib

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool to simulate SNPs into a FASTA Reference')
    parser.add_argument('-i','--input',dest='fname', action='store', required=True, help='FASTA file [Required]')
    parser.add_argument('-o','--out', dest='oname', action='store', required=True, help='FASTQ file with simulated reads from FASTA [Required]')
    parser.add_argument('-n','--num', dest='n', action='store', required=True, help='Read length for simulation [Required]')
    parser.add_argument('--prefix', dest='prefix', action='store', required=False, help='A prefix to add to the read id for distinguishing where reads came from [Optional]')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    args = parser.parse_args()
    #args = parser.parse_args(['-i', '/home/jfear/tmp/fa/fb551_si_fusions.fa', '-o', '/home/jfear/tmp/slide_test.fq', '-n', '95'])
    return(args)

def sliding_window(seq, n):
    """ Returns a generator with a sliding window """
    from itertools import islice
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) <= n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def main(args):
    num = int(args.n)
    with open(args.fname, 'r') as FA:
        with open(args.oname, 'w') as OUT:
            for record in SeqIO.parse(FA, 'fasta'):
                winGen = sliding_window(record.seq, num)
                start = 1
                end = num
                for winSeq in winGen:
                    currseq = ''.join(winSeq)
                    if len(currseq) >= 25:
                        if args.prefix:
                            currid = "%s_%s:%s-%s" % (args.prefix,record.id,start,end)
                        else:
                            currid = "%s:%s-%s" % (record.id,start,end)
                        currqual = ''.join(['g']*len(currseq))
                        curr = "@%s\n%s\n+\n%s\n" % (currid, currseq, currqual)
                        OUT.write(curr)
                    start = start + 1
                    end = end + 1

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
