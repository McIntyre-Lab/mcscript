#!/usr/bin/env python
import argparse 
import logging
import os.path
from Bio import SeqIO

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a fasta file and a 3-column bed file and calculates GC content.")
    parser.add_argument("-f", dest="fname", action='store', required=True, help="Name of reference file in FASTA format [Required]")
    parser.add_argument("-b", dest="bname", action='store', required=False, help="Name of 3-column BED file with regions of interest [Optional]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of output file in csv format [Required]")
    parser.add_argument("-g", dest="log", action='store', required=False, help="Name of log file [Optional]") 
    args = parser.parse_args()
    #args = parser.parse_args(['-r1', '/home/jfear/tmp/fq/r1.fq', '-r2', '/home/jfear/tmp/fq/r2.fq', '--outdir', '/home/jfear/tmp/files', '-o', '/home/jfear/tmp/files/counts.csv', '-t', '/home/jfear/tmp/files/cnts_table.tsv', '-g', '/home/jfear/tmp/files/test.log'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')

def getGit():
    """ This function will parse the current Git commit version. This will
    allow recording of exactly which script was used to create a given
    output."""
    import subprocess
    # get full path to script
    fullname = os.path.abspath(__file__)
    gitdir = os.path.dirname(fullname)
    label = subprocess.check_output(["git", "--git-dir="+gitdir+"/.git", "--work-tree="+gitdir,"rev-parse","HEAD"])
    return(label.rstrip(), gitdir)



if __name__=='__main__':

    # Get command line options
    args = getOptions()

    # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)

        # Add version information from git to the log so we always know exactly
        # what version was used
        git_status, gitdir = getGit()
        logging.info("Starting %s", __file__) 
        logging.info("Running script from  %s", gitdir) 
        logging.info("Git commit id: %s", git_status)
    else:
        setLogger(os.devnull,logging.INFO)








    logging.info("Script complete.")
