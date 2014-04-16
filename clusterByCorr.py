#!/usr/bin/env python
import argparse 
import os.path
import logging
import numpy as np
import scipy as sc
import pandas as pd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a dataset and collapses columns based on correlation coefficients.")
    parser.add_argument("-i", dest="fname", action='store', required=True, help="Input CSV file [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output CSV file [Required]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file [Optional]") 
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

def getGeneName(df):
    name = df.columns.values.to_list()[0].split('_')
    return(name)

def getMaskList(mask):
    myList = []
    for index, row in mask.iterrows():
        myList.append(np.where(row)[0])
    return(myList)

def merge(lsts):
    sets = [set(lst.tolist()) for lst in lsts]
    merged = 1
    while merged:
        merged = 0
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = 1
                    common |= x
            results.append(common)
        sets = results
    return(sets)

def alphaList():
    import string
    return(list(string.ascii_lowercase))

def createNewIso(geneName, mySet, myAlpha):
    myDict = dict()
    for ind in xrange(len(mySet)):
        name = geneName + '_' + myAlpha[ind]
        value = df.iloc[:, mySet[ind]].apply(np.mean, axis=1)
        myDict[name] = value.round(decimals=5)
    return(myDict)

def main(args):
    # Get file name and import into a data frame
    fname = '/home/jfear/tmp/b52_test.csv';
    df = pd.read_csv(fname, index_col=['patRIL', 'matRIL'])
    gname = getGeneName(df)

    # Calculate the correlation and figure out where all of the correlation coefficients is > 0.8
    myCorr = df.corr()
    mask = myCorr > 0.8

    # Get a list of locations that the mask was true
    myList = getMaskList(mask)

    # Merge the list of mask locations to get sets of isoforms gruped together
    mySet = merge(myList)

    # Now create new "isoform" groups with an alpha numeric index append it to the gene name
    myAlpha = alphaList()
    myDict = createNewIso(gname, mySet, myAlpha)

    # Combine these into a data frame and write out
    newDf = pd.DataFrame(myDict)
    newDf.to_csv('/home/jfear/tmp/b52_merged.csv')

if __name__ == '__main__':
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

