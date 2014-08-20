#!/usr/bin/env python
import argparse 
import os
import sys
import re
import logging

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a set of tabular files (e.g. CSV, TSV) and concatenates together.")
    parser.add_argument("-f", dest="fname", action='store', nargs='*', required=True, help="List of files to concatneate [Required]")
    parser.add_argument("--odir", dest="odir", action='store', required=False, help="Output directory, if not provided will create output in the original file's folder [Optional]")
    parser.add_argument("--oname", dest="outname", action='store', required=False, help="Name of output file, if not provided I will try to guess [Optional]")
    parser.add_argument("--header", dest="header", action='store_true', help="Indicate if the file in question has a header [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file [Optional]") 
    args = parser.parse_args()
    #args = parser.parse_args(['-f', '/home/jfear/tmp/*.txt', '-g', '/home/jfear/tmp/test.log', '--header'])
    return(args)

def setLogger(fname, loglevel, stream):
    """ Function to handle error logging """
    if stream:
        logging.basicConfig(stream=sys.stdout, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
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

def removeOldOutput(args):
    """ check if the old output exists and remove the old version """
    try:
        os.remove(args.oname)
        logging.info("Removed old version of output file: {0}".format(args.oname))
    except:
        pass

def createOutput(args):
    """ Use input file information to create an output file """
    if args.outname:
        args.oname = os.path.join(args.odir, args.outname)
    else:
        bname = os.path.splitext(os.path.basename(args.fname[0]))
        try:
            # Files names incremented: foo_1.csv, foo_2.csv
            sumName = re.sub('_\d+$', '_summary', bname[0])
        except:
            try:
                # Files name are incremented by fusion id: foo_S1000_SI.csv, foo_F1002_SI.csv
                sumName = re.sub('[SF]\d+_SI$', '_summary', bname[0])
            except:
                # Everything else just append summary
                sumName = bname[0] + '_summary'
        args.oname = os.path.join(args.odir, sumName + bname[1])
        logging.info("File will be output to: {0}".format(args.oname))

    # remove old output if it is there
    removeOldOutput(args)

def getHeader(args):
    """ Pull the header from the first file """
    with open(args.fname[0], 'r') as INPUT:
        args.header = INPUT.readline()

def main(args):
    # If an output directory was not given then put output in the same
    # directory as input
    if not args.odir:
        args.odir = os.path.dirname(args.fname[0])

    # Create an output file
    createOutput(args)

    # get the list of files to process
    logging.info("A total of {0} files will be combined".format(len(args.fname)))

    with open(args.oname, 'w') as OUT:
        # If header was passed pull the header line
        if args.header:
            logging.info("Pulling header line")
            getHeader(args)
            OUT.write(args.header)

        for currFile in args.fname:
            with open(currFile, 'r') as INPUT:
                INPUT.next()
                for row in INPUT:
                    OUT.write(row)

if __name__=='__main__':
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO, False)
    else:
        setLogger(os.devnull,logging.INFO, True)

    # Get Git information and start log
    git_status, gitdir = getGit()
    logging.info("Starting %s", __file__) 
    logging.info("Running script from  %s", gitdir) 
    logging.info("Git commit id: %s", git_status)

    # Run Main part of the script
    main(args)
    logging.info("Script complete.")
