#!/usr/bin/env python

# Built-in packages
import argparse 
import os
import sys
import re
import logging

# McLab Packages
import mclib_Python as mclib

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a set of tabular files (e.g. CSV, TSV) and concatenates together.")
    parser.add_argument("-f", dest="fname", action='store', nargs='*', required=True, help="List of files to concatneate [Required]")
    parser.add_argument("--odir", dest="odir", action='store', required=False, help="Output directory, if not provided will create output in the original file's folder [Optional]")
    parser.add_argument("--oname", dest="outname", action='store', required=False, help="Name of output file, if not provided I will try to guess [Optional]")
    parser.add_argument("--header", dest="header", action='store_true', help="Indicate if the file in question has a header [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file [Optional]") 
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    args = parser.parse_args()
    #args = parser.parse_args(['-f', '/home/jfear/tmp/*.txt', '-g', '/home/jfear/tmp/test.log', '--header'])
    return(args)

def removeOldOutput(args):
    """ check if the old output exists and remove the old version """
    try:
        os.remove(args.oname)
        logger.info("Removed old version of output file: {0}".format(args.oname))
    except:
        pass

def createOutput(args):
    """ Use input file information to create an output file """
    if args.outname:
        args.oname = os.path.join(args.odir, args.outname)
    else:
        bname = os.path.splitext(os.path.basename(args.fname[0]))

        # try to guess how files are named
        # REGEX patterns
        incremental = re.compile('_\d+$')   ## Files names are incremented: foo_1.csv, foo_2.csv
        fusion = re.compile('_[FS]\d+_SI$') ## Files name contain fusion id: foo_S1000_SI.csv, foo_F1002_SI.csv

        if re.search(incremental, bname[0]):
            sumName = re.sub(incremental, '_summary', bname[0])
        elif re.search(fusion, bname[0]):
            sumName = re.sub(fusion, '_summary', bname[0])
        else:
            # Everything else just append summary
            sumName = bname[0] + '_summary'

        args.oname = os.path.join(args.odir, sumName + bname[1])
        logger.info("File will be output to: {0}".format(args.oname))

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
    logger.info("A total of {0} files will be combined".format(len(args.fname)))

    with open(args.oname, 'w') as OUT:
        # If header was passed pull the header line
        if args.header:
            logger.info("Pulling header line")
            getHeader(args)
            OUT.write(args.header)

        for currFile in args.fname:
            with open(currFile, 'r') as INPUT:
                if args.header:
                    INPUT.next()
                for row in INPUT:
                    OUT.write(row)

if __name__=='__main__':
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
