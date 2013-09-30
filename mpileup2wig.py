#!/usr/bin/env python

import logging
import argparse 

def getOptions():
    """ Function to pull in arguments """

    parser = argparse.ArgumentParser(description="Converts a mpileup file into a Wiggle file. Can also do coordinate conversion if the appropriate BED file is included.")
    parser.add_argument("-i", "--input", dest="fname", action='store', required=True, help="Input MPILEUP file [Required]")
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Output WIG file [Required]")
    parser.add_argument("-b", "--bed", dest="bed", action='store', required=False, help="A 4+ column BED file for coordinate conversion or splitting junctions. [Optional]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    #args = parser.parse_args(['-i','/home/jfear/tmp/mpileup/small_junc.mpileup', '-o', '/home/jfear/tmp/mpileup/small_junc.wig','-b', '/home/jfear/tmp/mpileup/fb551_canonical_200bpJunctions.bed'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def rawMpileup2Wig(fname,oname):
    with open(oname, 'w') as OUT:
        with open(fname, 'r') as IN:
            currChrom = ''
            for row in IN:
                chrom, pos, refBase, count, readBase, readQual = row.split('\t')
                if chrom == currChrom:
                    OUT.write(str(pos) + '\t' + str(count) + "\n")
                else:
                    OUT.write("variableStep chrom=" + str(chrom) + "\n")
                    OUT.write(str(pos) + '\t' + str(count) + "\n")
                    currChrom = chrom
def makeCoordDict(bedName,flagJunc):
    bDict = dict()
    with open(bedName, 'r') as BED:
        for row in BED:
            rowList = row.split('\t')
            if len(rowList) > 4:
                chrom, cstart, cend, fusion, score, strand, tstart, tend, color, blockNum, sizes, starts = rowList
                sizeList = sizes.split(',')
                startList = starts.split(',')
                globalS1 = int(cstart) + 1 # convert from 0-based to 1-based start
                globalS2 = int(cstart) + int(startList[1]) + 1 # convert from 0-based to 1-based start
                bDict[fusion] = (int(sizeList[0]), chrom, globalS1, globalS2)
                if not flagJunc:
                    flagJunc = 1
            elif len(rowList) == 4:
                chrom, cstart, cend, fusion = rowList
                globalS1 = int(cstart) + 1 # convert from 0-based to 1-based start
                bDict[fusion.rstrip()] = (chrom, globalS1)
            else:
                raise ValueError("BED file must have at least 4 columns")
    return(bDict,flagJunc)

def convMpileup2Wig(fname,oname,bedDict,flagJunc):
    with open(oname, 'w') as OUT:
        pileList = list()
        with open(fname, 'r') as IN:
            if flagJunc:
                for row in IN:
                    fusion, pos, refBase, count, readBase, readQual = row.split('\t')
                    splitFus = fusion.split('|')
                    joinFus = '|'.join(splitFus[0:2])
                    currFusion = bedDict[joinFus]
                    juncSize = int(currFusion[0])
                    chrom = currFusion[1]
                    globalS1 = int(currFusion[2])
                    globalS2 = int(currFusion[3])
                    if int(pos) < juncSize + 1:
                        currPos = globalS1 + int(pos)
                    else:
                        localPos = int(pos) - juncSize
                        currPos = globalS2 + localPos
                    pileList.append((str(chrom), currPos, int(count)))
            else:
                for row in IN:
                    fusion, pos, refBase, count, readBase, readQual = row.split('\t')
                    currFusion = bedDict[fusion]
                    chrom = currFusion[0]
                    currPos = currFusion[1] + int(pos)
                    pileList.append((str(chrom), currPos, int(count)))

        pileList.sort(key=lambda tup: (tup[0],tup[1]))

        currChrom = ''
        for row in pileList:
            chrom, pos, count = row
            if chrom == currChrom:
                OUT.write(str(pos) + '\t' + str(count) + "\n")
            else:
                OUT.write("variableStep chrom=" + str(chrom) + "\n")
                OUT.write(str(pos) + '\t' + str(count) + "\n")
                currChrom = chrom

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)

    if not args.bed:
        logging.info("Converting '%s' to '%s'" % (args.fname,args.oname))
        rawMpileup2Wig(args.fname,args.oname)
        logging.info("Finished converting '%s' to '%s'" % (args.fname,args.oname))
    else:
        flagJunc = 0
        bDict,flagJunc = makeCoordDict(args.bed,flagJunc)
        logging.info("Converting '%s' to '%s'" % (args.fname,args.oname))
        convMpileup2Wig(args.fname,args.oname,bDict,flagJunc)
        logging.info("Finished converting '%s' to '%s'" % (args.fname,args.oname))


if __name__=='__main__':
    main()
    logging.info("Script complete.")
