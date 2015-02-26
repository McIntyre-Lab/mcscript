#!/usr/bin/env python

#
#    DESCRIPTION: This script takes a fASTQ file and subsets it using a list of Fastq IDs.
#    It creates a new FASTQ file containing only the reads in the subset list.
#
#    AUTHOR: Alison Morse
#

import os,csv,sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Parse command line arguments
parser = argparse.ArgumentParser(description='Subset FQ file using a list of FQ IDs')
parser.add_argument("-i", "--input_file", dest="input_file", required=True, help="fq file to subset from")
parser.add_argument("-l", "--id_file", dest="id_file", required=True, help="list if IDs want to subset")
parser.add_argument("-o", "--output_file", dest="output_file", required=True, help="output file for the subsetted fq reads")
args = parser.parse_args()

wanted = set()
with open(args.id_file, "r") as id:
    for line in id:
	line = line.strip()
	if line !="":
            wanted.add(line)

fq_sequences = SeqIO.parse(open(args.input_file, "r"), "fastq")
with open(args.output_file, "w") as out:
    for seq in fq_sequences:
        seq.id = seq.id[:-2]
	if seq.id in wanted:
            SeqIO.write([seq], out, "fastq")

