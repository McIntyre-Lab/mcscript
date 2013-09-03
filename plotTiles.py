#!/usr/bin/env python
import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes coordinate table created by fastqDumpCoords.py and plots the flowcell location for the top N duplicated reads.")
    parser.add_argument("-i", "--input", dest="input", action='store', required=True, help="fastqDumpCoords.py coordinate table [Required]")
    parser.add_argument("-o", "--outdir", dest="out", action='store', required=True, help="Directory to output PNGs [Required]")
    parser.add_argument("-N",  dest="num", action='store', default=10, required=False, help="Number of sequences to plot [Default 10]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File") 
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def incZ(x):
    one = x['plane'] - 1
    two = x['tileNum'] - 1
    three = x['swath'] - 1
    zarray[one,two,three] += 1


def plotFlow(zarray, pp, maxTileNum, seq):
    column_labels = list('123')
    row_labels = list(range(1,maxTileNum+1))
    fig, (ax1, ax2) = plt.subplots(1,2,sharey=True)
    heatmap1 = ax1.pcolor(zarray[0], cmap=plt.cm.Blues)
    heatmap2 = ax2.pcolor(zarray[1], cmap=plt.cm.Blues)

    # put the major ticks at the middle of each cell
    ax1.set_xticks(np.arange(zarray[0].shape[1])+0.5, minor=False)
    ax1.set_yticks(np.arange(zarray[0].shape[0])+0.5, minor=False)
    ax2.set_xticks(np.arange(zarray[1].shape[1])+0.5, minor=False)
    ax2.set_yticks(np.arange(zarray[1].shape[0])+0.5, minor=False)

    # Change axis labels to look better
    ax1.invert_yaxis()

    ax1.set_xticklabels(column_labels, minor=False)
    ax1.set_yticklabels(row_labels, minor=False)
    ax2.set_xticklabels(column_labels, minor=False)
    ax2.set_yticklabels(row_labels, minor=False)

    # add titles
    ax1.set_title("top")
    ax2.set_title("bottom")

    # Add gradient bar to figure
    fig.colorbar(heatmap1)
    plt.suptitle(seq)
    pp.savefig()

def main():
    """ MAIN Function to execute everything """
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    fname = os.path.basename(args.input)
    myname = os.path.splitext(fname)[0]
    pname = os.path.join(args.out, myname + '.pdf')
    pp = PdfPages(pname)

    logging.info("Importing coordinate table.")
    df = pd.read_csv(args.input)
    logging.info("Finished importing coordinate table.")

    logging.info("Summarizing sequence counts.")
    counts = pd.value_counts(df['sequence'])
    logging.info("Finished summarizing sequence counts.")

    for i in xrange(int(args.num)):
        seq = counts.index[i]
        subset = df[df['sequence'] == seq]
        maxTileNum = max(subset['tileNum'])
        global zarray
        zarray = np.zeros((2,maxTileNum,3))
        subset.apply(incZ, axis=1)
        plotFlow(zarray, pp, maxTileNum, seq)

    pp.close()


if __name__=='__main__':
    main()
    logging.info("Script complete.")
