#!/usr/bin/env python
# Standard Libraries
import argparse 
from collections import defaultdict
import itertools
import operator
import re
import os.path
import sys
import logging

# AddOn Libraries
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# McLab Libraries
import mclib_Python as mclib

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a single-end (SE) or paired-end (PE) file and splits out unique and duplicate reads.")
    parser.add_argument("-r1", dest="r1", action='store', required=True, help="Name of read 1 or SE FASTQ file [Required]")
    parser.add_argument("-r2", dest="r2", action='store', required=False, help="Name of read 2 FASTQ file. If SE leave blank.")
    parser.add_argument("--outdir", dest="odir", action='store', required=True, help="Directory to store FASTQ output files [Required]")
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Path and name of the summary table of read counts in csv format [Required]")
    parser.add_argument("-t", "--table", dest="tname", action='store', required=False, help="Path and name of the table showing each sequence and the number of reads with that sequence in tsv format [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file") 
    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 
    args = parser.parse_args()
    #args = parser.parse_args(['-r1', '/home/jfear/tmp/fq/r1.fq', '-r2', '/home/jfear/tmp/fq/r2.fq', '--outdir', '/home/jfear/tmp/files', '-o', '/home/jfear/tmp/files/counts.csv', '-t', '/home/jfear/tmp/files/cnts_table.tsv', '-g', '/home/jfear/tmp/files/test.log'])
    return(args)

def getName(odir, fname, fileNames):
    """ This function creates unique and duplicate dataset names. """
    bname = os.path.basename(fname)
    name1 = os.path.splitext(bname)[0]
    name2 = re.sub('1$','2', name1)
    # Read 1 output file names
    fileNames['uname'].append(os.path.join(odir, name1 + '_uniq.fq'))
    fileNames['duname'].append(os.path.join(odir, name1 + '_duplicate.fq'))
    fileNames['diname'].append(os.path.join(odir, name1 + '_distinct.fq'))
    fileNames['sduname'].append(os.path.join(odir, name1 + '_single_duplicate.fq'))

    # Read 2 output file names
    fileNames['uname'].append(os.path.join(odir, name2 + '_uniq.fq'))
    fileNames['duname'].append(os.path.join(odir, name2 + '_duplicate.fq'))
    fileNames['diname'].append(os.path.join(odir, name2 + '_distinct.fq'))
    fileNames['sduname'].append(os.path.join(odir, name2 + '_single_duplicate.fq'))

    # Unpaired output file names
    fileNames['upuniq'] = os.path.join(odir, name1 + '_' + name2 + '_unpaired_uniq.fq')
    fileNames['updup'] = os.path.join(odir, name1 + '_' + name2 + '_unpaired_duplicate.fq')
    fileNames['updist'] = os.path.join(odir, name1 + '_' + name2 + '_unpaired_distinct.fq')

def readFq(fname, readDict):
    """ Function to read a FQ file, I used the FastqGeneralIterator for its
    speeds since I don't need a full SeqIO object. """
    with open(fname, 'r') as FQ:
        for header, seq, qual in FastqGeneralIterator(FQ):
            if header.count(' '):
                # Header was created using CASAVA 1.8+ 
                (mhead, suphead) = header.split(' ')
                readDict[mhead].append((seq,qual))
            else:
                # Header was created using older versions of CASAVA
                header = re.sub('/[1-2]','',header)
                readDict[header].append((seq,qual))

def idDups(readDict):
    seqDict = defaultdict(list)
    for key in readDict:
        try:
            seq = readDict[key][0][0] + readDict[key][1][0]
        except:
            seq = readDict[key][0][0]
        seqDict[seq].append(key)
    return(seqDict)

def silentRemove(filename):
    """ Remove a file if it exists """
    try:
        os.remove(filename)
    except:
        pass

def buildOutSE(readDict, vlist):
    myout = ''
    for value in vlist:
        myout += '\n'.join(str(x) for x in ['@' + value, readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
    return(myout)

def buildOutPE(readDict, vlist):
    myout1 = ''
    myout2 = ''
    upout = ''
    for value in vlist:
        try:
            myout2 += '\n'.join(str(x) for x in ['@' + value + '/2', readDict[value][1][0],'+', readDict[value][1][1]]) + '\n'
            myout1 += '\n'.join(str(x) for x in ['@' + value + '/1', readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
        except:
            upout += '\n'.join(str(x) for x in ['@' + value , readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
    return(myout1, myout2, upout)

def buildDistinctOutSE(readDict, vlist):
    myout = ''
    value = vlist[0]
    myout += '\n'.join(str(x) for x in ['@' + value, readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
    return(myout)

def buildDistinctOutPE(readDict, vlist):
    myout1 = ''
    myout2 = ''
    upout = ''
    value = vlist[0]
    try:
        myout2 += '\n'.join(str(x) for x in ['@' + value + '/2', readDict[value][1][0],'+', readDict[value][1][1]]) + '\n'
        myout1 += '\n'.join(str(x) for x in ['@' + value + '/1', readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
    except:
        upout += '\n'.join(str(x) for x in ['@' + value , readDict[value][0][0],'+', readDict[value][0][1]]) + '\n'
    return(myout1, myout2, upout)

def writeTables(seqDict, stats, fileNames):
    # Create Output Count Table. This is the same output as Alison's original duplicate count script.
    percent_uniq_seq = float(stats['uniq_num']) / float(stats['total_num']) * 100
    percent_uniq_reads = float(stats['num_uniq_reads']) / float(stats['total_num']) * 100
    percent_dist_reads = float(stats['dist_num']) / float(stats['total_num']) * 100
    percent_dup_reads = float(stats['dup_num']) / float(stats['total_num']) * 100

    with open(fileNames['oname'], 'w') as ON:
        ON.write("file_name,total_reads,num_unique_reads,per_unique_reads,num_distinct_reads,per_distinct_reads,num_duplicate_reads,per_duplicate_reads\n")
        myout = [os.path.basename(fileNames['oname']),stats['total_num'],stats['num_uniq_reads'],percent_uniq_reads,stats['dist_num'],percent_dist_reads,stats['dup_num'],percent_dup_reads]
        ON.write(','.join([str(x) for x in myout])+"\n")


    if fileNames['tname']:
        counts = dict()
        for item in seqDict:
            counts[item] = len(seqDict[item])

        with open(fileNames['tname'], 'w') as OT:
            sorted_counts = sorted(counts.iteritems(), key=operator.itemgetter(1), reverse=True)
            for item in sorted_counts:
                OT.write(str(item[1]) + '\t' + str(item[0]).strip() + "\n")

def writeOutputSE(readDict, seqDict, fileNames):
    UNIQ = open(fileNames['uname'][0], 'w')
    DUPS = open(fileNames['duname'][0], 'w')
    DIST = open(fileNames['diname'][0], 'w')
    SDUPS = open(fileNames['sduname'][0], 'w')

    stats = dict()
    stats['total_num'] = 0
    stats['num_uniq_reads'] = 0
    stats['uniq_num'] = 0
    stats['dup_num'] = 0
    stats['dist_num'] = 0

    for value in seqDict.values():
        if len(value) == 1:
            # Unique Reads Output
            myout = buildOutSE(readDict,value)
            UNIQ.write(myout)
            stats['total_num'] += 1
            stats['num_uniq_reads'] += 1
            stats['uniq_num'] += 1
        else:
            # Duplicated Reads Output
            myout = buildOutSE(readDict,value)
            DUPS.write(myout)
            stats['total_num'] += len(value)
            stats['dup_num'] += len(value)
            stats['uniq_num'] += 1
            myoutSdup = buildDistinctOutSE(readDict,value)
            SDUPS.write(myoutSdup)

        # Distinct Reads Output
        myout = buildDistinctOutSE(readDict,value)
        DIST.write(myout)
        stats['dist_num'] += 1

    # Close Output Files
    UNIQ.close()
    DUPS.close()
    DIST.close()
    SDUPS.close()

    # Write out summary and counts tables
    writeTables(seqDict, stats, fileNames)

def writeOutputPE(readDict, seqDict, fileNames):
    UNIQ1 = open(fileNames['uname'][0], 'w')
    DUPS1 = open(fileNames['duname'][0], 'w')
    DIST1 = open(fileNames['diname'][0], 'w')
    SDUPS1 = open(fileNames['sduname'][0], 'w')
    UNIQ2 = open(fileNames['uname'][1], 'w')
    DUPS2 = open(fileNames['duname'][1], 'w')
    DIST2 = open(fileNames['diname'][1], 'w')
    SDUPS2 = open(fileNames['sduname'][1], 'w')

    # Remove unpaired files if they already exists since I am doing append steps
    silentRemove(fileNames['upuniq'])
    silentRemove(fileNames['updup'])
    silentRemove(fileNames['updist'])

    stats = dict()
    stats['total_num'] = 0
    stats['num_uniq_reads'] = 0
    stats['uniq_num'] = 0
    stats['dup_num'] = 0
    stats['dist_num'] = 0

    # Raise a flag if there are reads that are missing thier mate-pair
    flag_unpaired = 0

    for value in seqDict.values():
        if len(value) == 1:
            # Unique Reads Output
            myout1, myout2, upout = buildOutPE(readDict,value)
            UNIQ1.write(myout1)
            UNIQ2.write(myout2)
            if upout:
                flag_unpaired = 1
                with open(fileNames['upuniq'], 'a') as UPU:
                    UPU.write(upout)
            stats['num_uniq_reads'] += 1
            stats['total_num'] += 1
            stats['uniq_num'] += 1
        else:
            # Duplicated Reads Output
            myout1, myout2, upout = buildOutPE(readDict,value)
            DUPS1.write(myout1)
            DUPS2.write(myout2)
            if upout:
                flag_unpaired = 1
                with open(fileNames['updup'], 'a') as UPD:
                    UPD.write(upout)
            stats['total_num'] += len(value)
            stats['dup_num'] += len(value)
            stats['uniq_num'] += 1
            myout1Sdups, myout2Sdups, upoutSdups = buildDistinctOutPE(readDict,value)
            SDUPS1.write(myout1Sdups)
            SDUPS2.write(myout2Sdups)

        # Distinct Reads Output
        myout1, myout2, upout = buildDistinctOutPE(readDict,value)
        DIST1.write(myout1)
        DIST2.write(myout2)
        if upout:
            flag_unpaired = 1
            with open(fileNames['updist'], 'a') as UPDIST:
                UPDIST.write(upout)
        stats['dist_num'] += 1

    if flag_unpaired == 1:
        logger.warn("Some of your reads are missing their mate pair!")

    # Close Output Files
    UNIQ1.close()
    DUPS1.close()
    DIST1.close()
    UNIQ2.close()
    DUPS2.close()
    DIST2.close()
    SDUPS1.close()
    SDUPS2.close()

    # Write out summary and counts tables
    writeTables(seqDict, stats, fileNames)

def main(args):
    """ MAIN Function to execute everything """
    odir = os.path.abspath(args.odir)

    # Construct output file names
    fileNames = defaultdict(list)
    fileNames['oname'] = args.oname
    fileNames['tname'] = args.tname
    getName(odir, args.r1, fileNames)

    # Read in first fastq file and create dictionary of read IDs
    readDict = defaultdict(list)
    logger.info("Reading '%s'" % (args.r1))
    readFq(args.r1,readDict)
    logger.info("Finished reading '%s'" % (args.r1))

    if args.r2:
        # Read in second fastq file and create dictionary of read IDs
        logger.info("Reading '%s'" % (args.r2))
        readFq(args.r2,readDict)
        logger.info("Finished reading '%s'" % (args.r2))

        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        seqDict = idDups(readDict)
        writeOutputPE(readDict, seqDict, fileNames)
    else:
        # Create a dictionary where the key is unique sequences and value is a list
        # of headers
        seqDict = idDups(readDict)
        writeOutputSE(readDict, seqDict, fileNames)

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
