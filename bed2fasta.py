#!/usr/bin/env python

#    
#    DESCRIPTION: This script uses a junctions BED file and reads the coordinates to slice out the sequences
#    out of the reference FASTA file. It creates a new FASTA file with all of these junctions' sequences.
#
#    AUTHOR: Chelsea Tymms
#



import os,csv,sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


# Parse command line arguments
parser = argparse.ArgumentParser(description='Using Cooridnates from a BED file, slice sequences out of a FASTA file')
parser.add_argument('-f','--fasta', dest='fa', action='store', required=True, help='A fasta file [Required]')
parser.add_argument('-b','--bed', dest='bed', action='store', required=True, help='A bed file for fusions [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for fusions in FASTA format [Required]')
args = parser.parse_args()

# Open FASTA file and read it into a list
with open(args.fa,"r") as FASTA:
    records = list(SeqIO.parse(FASTA, "fasta"))

# Open a BED file and Read it into a list
with open(args.bed) as BED:
    bedReader=csv.reader(BED,delimiter="\t")
    bedrecord=list(bedReader)

with open(args.out,'w') as FH:
    for row in bedrecord:
        chrom = row[0]
        start = int(row[1]) 
        end = int(row[2])
        name=row[3]        
        #Get the actual sequences and write the junction to the output FASTA
        for ch in records:
            if ch.id == chrom:
                sequ=Seq("")
                sequ=sequ+ch.seq[start:end]
                #current = SeqIO.SeqRecord(seq = sequ, id=name, description = '{0}:{1}--{2}'.format(chrom, start, end)) 
                current = SeqIO.SeqRecord(seq = sequ, id=name, description = '') 
                SeqIO.write(current,FH,"fasta")
