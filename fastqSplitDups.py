#!/usr/bin/env python
import argparse 
import collections
import re
import os.path
import logging
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a single-end (SE) or paired-end (PE) file and splits out unique and duplicate reads.")
    parser.add_argument("-r1", dest="r1", action='store', required=True, help="Name of read 1 or SE FASTQ file [Required]")
    parser.add_argument("-r2", dest="r2", action='store', required=False, help="Name of read 2 FASTQ file. If SE leave blank.")
    parser.add_argument("--outdir", dest="odir", action='store', required=True, help="Directory to store output files [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def getName(odir, fname):
    bname = os.path.basename(fname)
    name = os.path.splitext(bname)[0]
    uname = os.path.join(odir, name + '_uniq.fq')
    dname = os.path.join(odir, name + '_duplicate.fq')
    upuniq = os.path.join(odir, name + 'unpaired_uniq.fq')
    updup = os.path.join(odir, name + 'unpaired_duplicate.fq')
    return(uname, dname, upuniq, updup)

def readFq(fname,hdict):
    with open(fname, 'r') as FQ:
        for header, seq, qual in FastqGeneralIterator(FQ):
            header = re.sub('/[1-2]','',header)
            hdict[header].append((seq,qual))

def idDups(hdict):
    sdict = collections.defaultdict(list)
    for key in hdict:
        try:
            seq = hdict[key][0][0] + hdict[key][1][0]
        except:
            seq = hdict[key][0][0]
        sdict[seq].append(key)
    return(sdict)

def buildOutSE(hdict,vlist):
    myout = ''
    for value in vlist:
        myout += '\n'.join(str(x) for x in ['@' + value, hdict[value][0][0],'+', hdict[value][0][1]]) + '\n'
    return(myout)

def writeOutputSE(hdict,sdict,uname,dname):
    UNIQ = open(uname[0], 'w')
    DUPS = open(dname[0], 'w')

    for value in sdict.values():
        if len(value) == 1:
            myout = buildOutSE(hdict,value)
            UNIQ.write(myout)
        else:
            myout = buildOutSE(hdict,value)
            DUPS.write(myout)

def buildOutPE(hdict,vlist):
    myout1 = ''
    myout2 = ''
    upout = ''
    for value in vlist:
        try:
            myout2 += '\n'.join(str(x) for x in ['@' + value + '/2', hdict[value][1][0],'+', hdict[value][1][1]]) + '\n'
            myout1 += '\n'.join(str(x) for x in ['@' + value + '/1', hdict[value][0][0],'+', hdict[value][0][1]]) + '\n'
        except:
            upout += '\n'.join(str(x) for x in ['@' + value , hdict[value][0][0],'+', hdict[value][0][1]]) + '\n'

    return(myout1, myout2, upout)

def writeOutputPE(hdict,sdict,uname,dname,upuniq,updup):
    UNIQ1 = open(uname[0], 'w')
    DUPS1 = open(dname[0], 'w')
    UNIQ2 = open(uname[1], 'w')
    DUPS2 = open(dname[1], 'w')
    UPU = open(upuniq, 'w')
    UPD = open(updup, 'w')

    for value in sdict.values():
        if len(value) == 1:
            myout1, myout2, upout = buildOutPE(hdict,value)
            UNIQ1.write(myout1)
            UNIQ2.write(myout2)
            UPU.write(upout)
        else:
            myout1, myout2, upout = buildOutPE(hdict,value)
            DUPS1.write(myout1)
            DUPS2.write(myout2)
            UPD.write(upout)

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    odir = os.path.abspath(args.odir)

    # Construct output files for read 1 or SE reads
    uname1, dname1, upuniq, updup = getName(odir, args.r1)
    uname = [uname1]
    dname = [dname1]

    # Read in first fastq file.
    hdict = collections.defaultdict(list)
    logging.info("Reading '%s'" % (args.r1))
    readFq(args.r1,hdict)
    logging.info("Finished reading '%s'" % (args.r1))

    if args.r2:
        # Construct output files for read 2
        uname2, dname2, upuniq, updup = getName(odir, args.r2)
        uname.append(uname2)
        dname.append(dname2)

        # Read in second fastq file
        logging.info("Reading '%s'" % (args.r2))
        readFq(args.r2,hdict)
        logging.info("Finished reading '%s'" % (args.r2))

        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputPE(hdict,sdict,uname,dname,upuniq,updup)
    else:
        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputSE(hdict,sdict,uname,dname)


if __name__=='__main__':
    main()
    logging.info("Script complete.")
