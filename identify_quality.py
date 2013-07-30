#!/usr/bin/env python
# This script is a simplification of code found in the fastqident script by
# "Ryan C. Thompson". Its goal is to return what kind of quality score a file
# is using.
import argparse
from Bio import SeqIO
import itertools


def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool to identify a FASTQ files quality score.')
    parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A FASTQ file [Required]')
    args = parser.parse_args()
    return(args)

def get_qual(filename):
    # Set minium quality scores for the various types: http://en.wikipedia.org/wiki/Fastq
    possible_encodings = set(('phred33', 'phred64'))

    phred33_min = 33
    phred64_min = 59
    phred64_offset = phred64_min - phred33_min # 26
    maxQuality = 40

    # Set some place holders
    minSeen = 999
    maxSeen = 0

    # I am importing everything as a sanger so I have to off set the solexa and
    # illumina minumus, hence the threshold above.
    seqio = SeqIO.parse(filename, "fastq-sanger") 
    seqSlice = itertools.islice(seqio, 0, 5000, 100)

    for record in seqSlice:
        quals = record.letter_annotations["phred_quality"]
        minSeen = min(minSeen, min(quals))
        maxSeen = max(maxSeen, max(quals))

        if 'phred33' in possible_encodings and maxSeen > maxQuality:
            possible_encodings.remove('phred33')
        if 'phred64' in possible_encodings and minSeen < phred64_offset:
            return 'phred33'

        if len(possible_encodings) == 1:
            return possible_encodings.pop()
        elif len(possible_encodings) == 0:
            raise ValueError("Could not identify FASTQ file %s: eliminated all possible encodigns." % (filename,))


def main():
    args = getOptions()
    quality = get_qual(args.fq)
    print(quality)

if __name__ == '__main__':
    main()
