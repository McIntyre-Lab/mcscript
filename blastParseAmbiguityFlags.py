#!/usr/bin/env python
import argparse 
import logging
import os
from collections import defaultdict as dd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a single-end (SE) or paired-end (PE) file and splits out unique and duplicate reads.")
    parser.add_argument("-i", dest="fname", action='store', required=True, help="Input BLAST file in outfmt 6 [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Output file for counts in csv format [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    #args = parser.parse_args(['-i', '/home/jfear/tmp/blast/ambiguity_blast_fb551_non-redundant_fusions.tsv', '-o', '/home/jfear/tmp/blast/test.csv', '-g', '/home/jfear/tmp/blast/test.log'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    mydict = dd(dict)

    with open(args.oname, 'w') as OUT:
        header = ['fusion_id','flag_self','flag_ambig','cnt_ambig','ambig_fusion_cat']
        with open(args.fname,'r') as IN:
            logging.info('Reading Blast File and Building Dictionary.')
            for row in IN:
                query, hit = row.split('\t')[:2]
                if query == hit:
                    mydict[query]['self'] = 1
                else:
                    try:
                        mydict[query]['ambig'].append(hit) 
                    except:
                        mydict[query]['ambig'] = [hit] 

        logging.info('Wirting Output.')
        for key in mydict:
            myout = [key,'0','0','0','NA']
            if mydict[key]['self']:
                myout[1] = '1'
                try:
                    ambig = mydict[key]['ambig']
                    cnt = len(ambig)
                    myout[2] = '1'
                    myout[3] = str(cnt)
                    myout[4] = '|'.join([str(x) for x in ambig])
                except:
                    pass
            else:
                logging.error("Exonic region %s did not map to itself." % key)
            OUT.write(','.join(myout) + "\n")

if __name__=='__main__':
    main()
    logging.info("Script complete.")
