#!/usr/bin/env python
import argparse 
import collections
import itertools
import operator
import re
import os.path
import logging
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a single-end (SE) or paired-end (PE) file and splits out unique and duplicate reads.")
    parser.add_argument("-r1", dest="r1", action='store', required=True, help="Name of read 1 or SE FASTQ file [Required]")
    parser.add_argument("-r2", dest="r2", action='store', required=False, help="Name of read 2 FASTQ file. If SE leave blank.")
    parser.add_argument("--outdir", dest="odir", action='store', required=True, help="Directory to store FASTQ output files [Required]")
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Output file for counts in csv format [Required]")
    parser.add_argument("-t", "--table", dest="tname", action='store', required=True, help="Output table of intermediate size information [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    #args = parser.parse_args(['-r1', '/home/jfear/tmp/fq/r1.fq', '-r2', '/home/jfear/tmp/fq/r2.fq', '--outdir', '/home/jfear/tmp/files', '-o', '/home/jfear/tmp/counts.csv', '-t', '/home/jfear/tmp/cnts_table.tsv'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def getName(odir, fname):
    """ This function creates unique and duplicate dataset names. """
    bname = os.path.basename(fname)
    name = os.path.splitext(bname)[0]
    uname = os.path.join(odir, name + '_uniq.fq')
    dname = os.path.join(odir, name + '_duplicate.fq')
    return(uname, dname)

def getNameUnpaired(odir, fname1, fname2):
    """ This function creates the unpaired dataset names. """
    name1 = os.path.splitext(os.path.basename(fname1))[0]
    name2 = os.path.splitext(os.path.basename(fname2))[0]
    upuniq = os.path.join(odir, name1 + '_' + name2 + '_unpaired_uniq.fq')
    updup = os.path.join(odir, name1 + '_' + name2 + '_unpaired_duplicate.fq')
    return(upuniq, updup)

def readFq(fname,hdict):
    """ Function to read a FQ file, I used the FastqGeneralIterator for its
    speeds since I don't need a full SeqIO object. """
    with open(fname, 'r') as FQ:
        for header, seq, qual in FastqGeneralIterator(FQ):
            if header.count(' '):
                # Header was created using CASAVA 1.8+ 
                (mhead, suphead) = header.split(' ')
                hdict[mhead].append((seq,qual))
            else:
                # Header was created using older versions of CASAVA
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

def writeOutputSE(hdict,sdict,uname,dname,oname,tname):
    UNIQ = open(uname[0], 'w')
    DUPS = open(dname[0], 'w')

    total_num = 0
    num_uniq_reads = 0
    uniq_num = 0

    for value in sdict.values():
        if len(value) == 1:
            myout = buildOutSE(hdict,value)
            UNIQ.write(myout)
            total_num += 1
            num_uniq_reads += 1
            uniq_num += 1
        else:
            myout = buildOutSE(hdict,value)
            DUPS.write(myout)
            total_num += len(value)
            uniq_num += 1

    UNIQ.close()
    DUPS.close()

    percent_uniq_seq = float(uniq_num) / float(total_num) * 100
    percent_uniq_reads = float(num_uniq_reads) / float(total_num) * 100

    with open(oname, 'w') as ON:
        ON.write("file_name,total_reads,num_unique_seq,per_uniq_seq,num_unique_reads,per_uniq_reads\n")
        myout = [os.path.basename(oname),total_num,uniq_num,percent_uniq_seq,num_uniq_reads,percent_uniq_reads]
        ON.write(','.join([str(x) for x in myout])+"\n")

    counts = dict()
    for item in sdict:
        counts[item] = len(sdict[item])
    with open(tname, 'w') as OT:
        sorted_counts = sorted(counts.iteritems(), key=operator.itemgetter(1), reverse=True)
        for item in sorted_counts:
            OT.write(str(item[1]) + ' ' + str(item[0]).strip() + "\n")

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

def writeOutputPE(hdict,sdict,uname,dname,upuniq,updup,oname,tname):
    UNIQ1 = open(uname[0], 'w')
    DUPS1 = open(dname[0], 'w')
    UNIQ2 = open(uname[1], 'w')
    DUPS2 = open(dname[1], 'w')
    UPU = open(upuniq, 'w')
    UPD = open(updup, 'w')

    total_num = 0
    num_uniq_reads = 0
    uniq_num = 0

    for value in sdict.values():
        if len(value) == 1:
            myout1, myout2, upout = buildOutPE(hdict,value)
            UNIQ1.write(myout1)
            UNIQ2.write(myout2)
            UPU.write(upout)
            num_uniq_reads += 1
            total_num += 1
            uniq_num += 1
        else:
            myout1, myout2, upout = buildOutPE(hdict,value)
            DUPS1.write(myout1)
            DUPS2.write(myout2)
            UPD.write(upout)
            total_num += len(value)
            uniq_num += 1

    UNIQ1.close()
    DUPS1.close()
    UNIQ2.close()
    DUPS2.close()
    UPU.close()
    UPD.close()

    percent_uniq_seq = float(uniq_num) / float(total_num) * 100
    percent_uniq_reads = float(num_uniq_reads) / float(total_num) * 100
    with open(oname, 'w') as ON:
        ON.write("file_name,total_reads,num_unique_seq,per_uniq_seq,num_unique_reads,per_uniq_reads\n")
        myout = [os.path.basename(oname),total_num,uniq_num,percent_uniq_seq,num_uniq_reads,percent_uniq_reads]
        ON.write(','.join([str(x) for x in myout])+"\n")

    counts = dict()
    for item in sdict:
        counts[item] = len(sdict[item])
    with open(tname, 'w') as OT:
        sorted_counts = sorted(counts.iteritems(), key=operator.itemgetter(1), reverse=True)
        for item in sorted_counts:
            OT.write(str(item[1]) + ' ' + str(item[0]).strip() + "\n")


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
    uname1, dname1 = getName(odir, args.r1)
    uname = [uname1]
    dname = [dname1]

    # Read in first fastq file.
    hdict = collections.defaultdict(list)
    logging.info("Reading '%s'" % (args.r1))
    readFq(args.r1,hdict)
    logging.info("Finished reading '%s'" % (args.r1))

    if args.r2:
        # Construct output files for read 2
        uname2, dname2 = getName(odir, args.r2)
        upuniq, updup = getNameUnpaired(odir, args.r1, args.r2)
        uname.append(uname2)
        dname.append(dname2)

        # Read in second fastq file
        logging.info("Reading '%s'" % (args.r2))
        readFq(args.r2,hdict)
        logging.info("Finished reading '%s'" % (args.r2))

        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputPE(hdict,sdict,uname,dname,upuniq,updup,args.oname,args.tname)
    else:
        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputSE(hdict,sdict,uname,dname,args.oname,args.tname)


if __name__=='__main__':
    main()
    logging.info("Script complete.")
