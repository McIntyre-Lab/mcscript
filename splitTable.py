#!/usr/bin/env python

# Built-in packages
import argparse 
import os.path
import logging

# McLab Packages
import mclib

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a tabular file (e.g. CSV, TSV) and splits its rows up into multiple files.")
    parser.add_argument("-f", dest="fname", action='store', required=True, help="Path to tabular file [Required]")
    parser.add_argument("--prefix", dest="prefix", action='store', required=False, help="Output file prefix, if prefix not provided input file name will be used [Optional]")
    parser.add_argument("-o", dest="odir", action='store', required=False, help="Output directory, if not provided will create output in the original file's folder [Optional]")
    parser.add_argument("--header", dest="header", action='store_true', help="Indicate if the file in question has a header [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file [Optional]") 

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--nfiles", dest="nfiles", type=int, action='store', help="The number of files that you want to create")
    group.add_argument("--nlines", dest="nlines", type=int, action='store', help="The number of lines that you want in each output file")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 

    args = parser.parse_args()
    #args = parser.parse_args(['-f', '/home/jfear/tmp/test.csv', '--prefix', 'bob', '-g', '/home/jfear/tmp/test.log','--nfiles', '100', '--header'])
    return(args)

def file_len(fname):
    """ Count the number of lines in a file """
    with open(fname) as fn:
        for i, l in enumerate(fn):
            pass
    return i + 1

def nfiles(args, numRows):
    """ Split a tabular file into N files """
    with open(args.fname, 'r') as fname:
        if args.header:
            # If there is a header pull it and remove one from row count
            header = fname.next()
            numRows -= 1

        # How many blocks do we need to output entire file
        if numRows % args.nfiles == 0:
            # The number of rows divides evenly by the number of files
            numBlocks = numRows / args.nfiles

        elif numRows <= args.nfiles :
            # The number of rows is less than or equal to the number of files
            numBlocks = 1
            logger.error("You are trying to split a file with {0} rows into {1} files.".format(numRows, args.nfiles))
            exit(1)

        else:
            # The number of rows are not evenly divisible
            numBlocks = numRows / args.nfiles + 1

        for fileNum in xrange(args.nfiles):
            with open(os.path.join(args.odir,"{0}_{1}.csv".format(args.prefix, fileNum)), 'w') as OUT:
                OUT.write(header)
                for J in xrange(numBlocks):
                    try:
                        row = fname.next()
                        OUT.write(row)
                    except:
                        break

def nlines(args):
    """ Split a tabular file such that all output files have N rows """
    header = ''
    writeFlag = 0
    fileNum = 1
    with open(args.fname, 'r') as fname:
        for index, row in enumerate(fname):
            if index == 0 and args.header:
                header = row
            else:
                if writeFlag > 0 and writeFlag < args.nlines:
                    OUT.write(row)
                    writeFlag += 1
                else:
                    writeFlag = 0
                    OUT = open(os.path.join(args.odir,"{0}_{1}.csv".format(args.prefix, fileNum)), 'w')
                    fileNum += 1
                    OUT.write(header)
                    OUT.write(row)
                    writeFlag += 1
    logger.info("Created {0} files.".format(fileNum))

def main(args):
    if not args.odir:
        args.odir = os.path.dirname(args.fname)

    if not args.prefix:
        args.prefix = os.path.splitext(os.path.basename(args.fname))[0]

    if args.nfiles:
        logger.info("Splitting into {0} files.".format(args.nfiles))
        numRows = file_len(args.fname)
        nfiles(args, numRows)
    else:
        logger.info("Creating split files eah with {0} rows.".format(args.nlines))
        nlines(args)

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
