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
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Path and name of the summary table of read counts in csv format [Required]")
    parser.add_argument("-t", "--table", dest="tname", action='store', required=True, help="Path and name of the table showing each sequence and the number of reads with that sequence in tsv format [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file") 
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


def getName(odir, fname):
    """ This function creates unique and duplicate dataset names. """
    bname = os.path.basename(fname)
    name = os.path.splitext(bname)[0]
    uname = os.path.join(odir, name + '_uniq.fq')
    duname = os.path.join(odir, name + '_duplicate.fq')
    diname = os.path.join(odir, name + '_distinct.fq')
    return(uname, duname, diname)

def getNameUnpaired(odir, fname1, fname2):
    """ This function creates the unpaired dataset names. """
    name1 = os.path.splitext(os.path.basename(fname1))[0]
    name2 = os.path.splitext(os.path.basename(fname2))[0]
    upuniq = os.path.join(odir, name1 + '_' + name2 + '_unpaired_uniq.fq')
    updup = os.path.join(odir, name1 + '_' + name2 + '_unpaired_duplicate.fq')
    updist = os.path.join(odir, name1 + '_' + name2 + '_unpaired_distinct.fq')
    return(upuniq, updup, updist)

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

def silentRemove(filename):
    """ Remove a file if it exists """
    try:
        os.remove(filename)
    except:
        pass

def buildOutSE(hdict,vlist):
    myout = ''
    for value in vlist:
        myout += '\n'.join(str(x) for x in ['@' + value, hdict[value][0][0],'+', hdict[value][0][1]]) + '\n'
    return(myout)

def buildDistinctOutSE(hdict,vlist):
    myout = ''
    for value in vlist:
        myout += '\n'.join(str(x) for x in ['@' + value, hdict[value][0][0],'+', hdict[value][0][1]]) + '\n'
        return(myout)

def writeOutputSE(hdict,sdict,uname,duname,diname,oname,tname):
    UNIQ = open(uname[0], 'w')
    DUPS = open(duname[0], 'w')
    DIST = open(diname[0], 'w')

    total_num = 0
    num_uniq_reads = 0
    uniq_num = 0
    dup_num = 0
    dist_num = 0

    for value in sdict.values():
        if len(value) == 1:
            # Unique Reads Output
            myout = buildOutSE(hdict,value)
            UNIQ.write(myout)
            total_num += 1
            num_uniq_reads += 1
            uniq_num += 1
        else:
            # Duplicated Reads Output
            myout = buildOutSE(hdict,value)
            DUPS.write(myout)
            total_num += len(value)
            dup_num += len(value)
            uniq_num += 1

        # Distinct Reads Output
        myout = buildDistinctOutSE(hdict,value)
        DIST.write(myout)
        dist_num += 1

    # Close Output Files
    UNIQ.close()
    DUPS.close()
    DIST.close()

    # Create Output Count Table. This is the same output as Alison's original duplicate count script.
    percent_uniq_seq = float(uniq_num) / float(total_num) * 100
    percent_uniq_reads = float(num_uniq_reads) / float(total_num) * 100
    percent_dist_reads = float(dist_num) / float(total_num) * 100
    percent_dup_reads = float(dup_num) / float(total_num) * 100

    with open(oname, 'w') as ON:
        ON.write("file_name,total_reads,num_unique_reads,per_unique_reads,num_distinct_reads,per_distinct_reads,num_duplicate_reads,per_duplicate_reads\n")
        myout = [os.path.basename(oname),total_num,num_uniq_reads,percent_uniq_reads,dist_num,percent_dist_reads,dup_num,percent_dup_reads]
        ON.write(','.join([str(x) for x in myout])+"\n")

    counts = dict()
    for item in sdict:
        counts[item] = len(sdict[item])
    with open(tname, 'w') as OT:
        sorted_counts = sorted(counts.iteritems(), key=operator.itemgetter(1), reverse=True)
        for item in sorted_counts:
            OT.write(str(item[1]) + '\t' + str(item[0]).strip() + "\n")

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

def buildDistinctOutPE(hdict,vlist):
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

def writeOutputPE(hdict,sdict,uname,duname,diname,upuniq,updup,updist,oname,tname):
    UNIQ1 = open(uname[0], 'w')
    DUPS1 = open(duname[0], 'w')
    DIST1 = open(diname[0], 'w')
    UNIQ2 = open(uname[1], 'w')
    DUPS2 = open(duname[1], 'w')
    DIST2 = open(diname[1], 'w')

    # Remove unpaired files if they already exists since I am doing append steps
    silentRemove(upuniq)
    silentRemove(updup)
    silentRemove(updist)

    total_num = 0
    num_uniq_reads = 0
    uniq_num = 0
    dup_num = 0
    dist_num = 0

    # Raise a flag if there are reads that are missing thier mate-pair
    flag_unpaired = 0

    for value in sdict.values():
        if len(value) == 1:
            # Unique Reads Output
            myout1, myout2, upout = buildOutPE(hdict,value)
            UNIQ1.write(myout1)
            UNIQ2.write(myout2)
            if upout:
                flag_unpaired = 1
                with open(upuniq, 'a') as UPU:
                    UPU.write(upout)
            num_uniq_reads += 1
            total_num += 1
            uniq_num += 1
        else:
            # Duplicated Reads Output
            myout1, myout2, upout = buildOutPE(hdict,value)
            DUPS1.write(myout1)
            DUPS2.write(myout2)
            if upout:
                flag_unpaired = 1
                with open(updup, 'a') as UPD:
                    UPD.write(upout)
            total_num += len(value)
            dup_num += len(value)
            uniq_num += 1

        # Distinct Reads Output
        myout1, myout2, upout = buildDistinctOutPE(hdict,value)
        DIST1.write(myout1)
        DIST2.write(myout2)
        if upout:
            flag_unpaired = 1
            with open(updist, 'a') as UPDIST:
                UPDIST.write(upout)
        dist_num += 1

    if flag_unpaired == 1:
        logging.warn("Some of your reads are missing their mate pair!")

    # Close Output Files
    UNIQ1.close()
    DUPS1.close()
    DIST1.close()
    UNIQ2.close()
    DUPS2.close()
    DIST2.close()

    # Create Output Count Table. This is the same output as Alison's original duplicate count script.
    percent_uniq_seq = float(uniq_num) / float(total_num) * 100
    percent_uniq_reads = float(num_uniq_reads) / float(total_num) * 100
    percent_dist_reads = float(dist_num) / float(total_num) * 100
    percent_dup_reads = float(dup_num) / float(total_num) * 100

    with open(oname, 'w') as ON:
        ON.write("file_name,total_reads,num_unique_reads,per_unique_reads,num_distinct_reads,per_distinct_reads,num_duplicate_reads,per_duplicate_reads\n")
        myout = [os.path.basename(oname),total_num,num_uniq_reads,percent_uniq_reads,dist_num,percent_dist_reads,dup_num,percent_dup_reads]
        ON.write(','.join([str(x) for x in myout])+"\n")

    counts = dict()
    for item in sdict:
        counts[item] = len(sdict[item])
    with open(tname, 'w') as OT:
        sorted_counts = sorted(counts.iteritems(), key=operator.itemgetter(1), reverse=True)
        for item in sorted_counts:
            OT.write(str(item[1]) + '\t' + str(item[0]).strip() + "\n")


def main(args):
    """ MAIN Function to execute everything """
    odir = os.path.abspath(args.odir)

    # Construct output files for read 1 or SE reads
    uname1, duname1, diname1 = getName(odir, args.r1)
    uname = [uname1]
    duname = [duname1]
    diname = [diname1]

    # Read in first fastq file.
    hdict = collections.defaultdict(list)
    logging.info("Reading '%s'" % (args.r1))
    readFq(args.r1,hdict)
    logging.info("Finished reading '%s'" % (args.r1))

    if args.r2:
        # Construct output files for read 2
        uname2, duname2, diname2 = getName(odir, args.r2)
        upuniq, updup, updist = getNameUnpaired(odir, args.r1, args.r2)
        uname.append(uname2)
        duname.append(duname2)
        diname.append(diname2)

        # Read in second fastq file
        logging.info("Reading '%s'" % (args.r2))
        readFq(args.r2,hdict)
        logging.info("Finished reading '%s'" % (args.r2))

        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputPE(hdict,sdict,uname,duname,diname,upuniq,updup,updist,args.oname,args.tname)
    else:
        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        sdict = idDups(hdict)
        writeOutputSE(hdict,sdict,uname,duname,diname,args.oname,args.tname)


if __name__=='__main__':
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    # Get Git information and start log
    git_status, gitdir = getGit()
    logging.info("Starting %s", __file__) 
    logging.info("Running script from  %s", gitdir) 
    logging.info("Git commit id: %s", git_status)

    # Run Main part of the script
    main(args)

    logging.info("Script complete.")
