#!/usr/bin/env python
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool to simulate SNPs into a FASTA Reference')
    parser.add_argument('-i','--input',dest='fname', action='store', required=True, help='FASTA file [Required]')
    parser.add_argument('-o','--out', dest='oname', action='store', required=True, help='FASTQ file with simulated reads from FASTA [Required]')
    parser.add_argument('-n','--num', dest='n', action='store', required=True, help='Read length for simulation [Required]')
    parser.add_argument('-g','--log',dest='log',action='store',required=False, help='Create an error log')
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


def main():
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

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
                        currid = "%s:%s-%s" % (record.id,start,end)
                        currqual = ''.join(['g']*len(currseq))
                        curr = "@%s\n%s\n+\n%s\n" % (currid, currseq, currqual)
                        OUT.write(curr)
                    start = start + 1
                    end = end + 1

if __name__ == '__main__':
    main()
    logging.info("Script complete")
