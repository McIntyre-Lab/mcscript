#!/usr/bin/env python
import os
import logging
import argparse 
from Bio.Blast import NCBIXML

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Parse BLASTs XML format to a basic table. This was originally done for the Bloom project, the script may need addjusted to output more information.")
    parser.add_argument("-i", "--input", dest="fname", action='store', required=True, help="Input XML file [Required]")
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Output CSV file [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    """ MAIN Function to execute everything """

    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    with open(args.fname, 'r') as xml_handle:
        with open(args.out, 'w') as OUT:
            OUT.write(','.join(['query_name','query_seq','target_name','target_id','percent_identity','score','evalue'])+"\n")
            blast_records = NCBIXML.parse(xml_handle)
            for record in blast_records:
                query_name = record.query
                query_length = record.query_length
                for aln in record.alignments:
                    hit_def = aln.hit_def
                    hit_id = aln.hit_id
                    for hsp in aln.hsps:
                        query_seq = hsp.query
                        aln_len = hsp.align_length
                        per_ident = float(aln_len) / float(query_length) * 100
                        score = hsp.score
                        evalue = hsp.expect
                        OUT.write(','.join([str(x) for x in [query_name,query_seq,hit_def,hit_id,per_ident,score,evalue]]) + "\n")


if __name__=='__main__':
    main()
    logging.info("Script complete.")
