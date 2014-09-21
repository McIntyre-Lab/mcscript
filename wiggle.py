import os
import logging
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib import gridspec

class GeneModel(object):
    """ Class to construct a gene model """
    def __init__(self, geneObj, height=2):
        """ Arguments:
        geneObj (obj) = a gene object created from a subclass of _Gene in mclib.gff 
        height (int) = the height of the exon model

        Attributes:
        xLoc (tuple) = Gene start and end coordinates.
        yLoc (list) = List of y-coordinates for plotting each transcript on a different row
        patches (list) = List of gene models as matplotlib patches

        """
        self.xLoc = (geneObj.start, geneObj.end)

        # Each transcript needs to be plotted on a different row. Create a list of
        # y-coordinates for plotting.
        self.yLoc = self._get_y(1, 5, geneObj.transCnt)

        self.patches = []

        # Make list of transcript IDs and sort them depending on which strand the gene
        # is on. Sort by start if on + strand and by end if on - strand.
        tsList = geneObj.transcript.keys()
        if geneObj.strand == '-':
            tsList.sort(key=lambda x: geneObj.transcript[x]['tsEnd'])
        else:
            tsList.sort(key=lambda x: geneObj.transcript[x]['tsStart'])

        # Draw the gene model as matplotlib patches. Iterate through the list
        # of transcripts and draw UTR and CDS if they are available, otherwise
        # draw exons.
        for index, ts in enumerate(tsList):
            try:
                # Add 3' UTR
                self._build_patch(geneObj.transcript[ts]['utr'][0], index, height=height, color='grey')

                # Add CDS
                self._build_patch(geneObj.transcript[ts]['cds'], index, height=height)

                # Add 5' UTR
                self._build_patch(geneObj.transcript[ts]['utr'][1], index, height=height, color='grey')
            except:
                # Just plot exons
                self._build_patch(geneObj.transcript[ts]['exons'], index, height=height)

            # Build the introns
            self._build_patch(geneObj.transcript[ts]['introns'], index, height=height, intron=True)

    def _build_patch(self, annoList, index, height=2, color='black', intron=False):
        """ Create a rectangle patch to represent a exon 
        Arguments:
        annoList (list) = list of start end coordinates for gene annotations
        index (int) = number of which transcript we are on, so I can plot different transcripts on different lines
        height (int) = height of the exons in the gen model
        color (str) = color for exon blocks
        intron (bool) = indicate if annoList are introns
        """
        for coord in annoList:
            start = coord[0] - 1    # coordinate convert from 1-based annotation to 0-based wiggle
            end = coord[1]
            width = end - start
            if intron:
                # Put intron in the middle of the exon
                yloc = self.yLoc[index] + height / 2 - 0.025 
                self.patches.append(Rectangle((start, yloc), width, 0.05, fill=True, color=color))
            else:
                self.patches.append(Rectangle((start, self.yLoc[index]), width, height, fill=True, color=color))

    def _get_y(self, start, step, count):
        """ Return a list of numbers equally spaced by a step size
        Arguments:
        start (int) = starting location, typically use 1
        step (int) = step size to space the number
        count (int) = number of equally spaced number to produce
        """
        yList = []
        for i in range(0,count):
            yList.append(start)
            start += step
        return yList

def plot_wiggle(pileDict, outName, chrom, start, end, geneModel=None, fusionModel=None, variantPos=None, title=None):
    """ Function to construct a wiggle plot with gene models
    Arguments:
    pileDict (dict) = where keys are genome coordinates and values are counts at that position
    outName (str) = name of the output png
    geneModel (obj) = a mclib/wiggle.py GeneModel object
    fusionModel (obj) = 
    variantPos (List) = a list of positions that have a variant of interest
    title (str) = a title for the plot
    """

    # Set up the figure
    fig = plt.figure(figsize=(20, 10))

    # Create multiple rows for each subplot
    if geneModel and variantPos and fusionModel:
        gs = gridspec.GridSpec(4, 1, height_ratios=[1, .08, 1, .5])
        ax1 = plt.subplot(gs[0])
        vax = plt.subplot(gs[1], sharex=ax1)
        gax = plt.subplot(gs[2], sharex=ax1)
        fax = plt.subplot(gs[3], sharex=ax1)

    elif geneModel and variantPos:
        gs = gridspec.GridSpec(3, 1, height_ratios=[1, .08, 1])
        ax1 = plt.subplot(gs[0])
        vax = plt.subplot(gs[1], sharex=ax1)
        gax = plt.subplot(gs[2], sharex=ax1)

    elif geneModel and fusionModel:
        gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, .5])
        ax1 = plt.subplot(gs[0])
        gax = plt.subplot(gs[1], sharex=ax1)
        fax = plt.subplot(gs[2], sharex=ax1)

    elif variantPos and fusionModel:
        gs = gridspec.GridSpec(3, 1, height_ratios=[1, .08,.5])
        ax1 = plt.subplot(gs[0])
        vax = plt.subplot(gs[1], sharex=ax1)
        fax = plt.subplot(gs[2], sharex=ax1)

    elif geneModel:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])
        gax = plt.subplot(gs[1], sharex=ax1)

    elif variantPos:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, .08])
        ax1 = plt.subplot(gs[0])
        vax = plt.subplot(gs[1], sharex=ax1)

    elif fusionModel:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,.5])
        ax1 = plt.subplot(gs[0])
        fax = plt.subplot(gs[1], sharex=ax1)

    else:
        ax1 = plt.subplot(111)

    # Make wiggle plot
    ## ax1 is the wiggle
    ax1.set_ylim(0, max(pileDict.values())+150)
    ax1.set_xlim(start-300, end+300)
    n, bins, patches = ax1.hist(pileDict.keys(), bins=len(pileDict.keys()), weights=pileDict.values())

    if variantPos:
        # Plot variants
        logging.debug('Creating variant plot')
        ## vax will be the variant subplot
        vax.scatter(variantPos, [0.5]*len(variantPos), marker="^")
        vax.axis('off')  # Hide y-axis on gene model plot

    if geneModel:
        # Plot gene models
        logging.debug('Creating geneModel plot')
        ## gax will be the gene model subplot
        gax.set_ylim(max(geneModel.yLoc)+3, min(geneModel.yLoc)-3)
        gax.axis('off')  # Hide y-axis on gene model plot

        # Put the gene models together and create their plots
        p = PatchCollection(geneModel.patches, match_original=True)
        gax.add_collection(p)

    if fusionModel:
        # Plot gene models
        logging.debug('Creating fusion plot')
        ## fax will be the fusion model subplot
        fax.set_ylim(min(geneModel.yLoc), max(geneModel.yLoc)+5)
        gax.axis('off')  # Hide y-axis on gene model plot

        # Put the fusion models together and create their plots
        p = PatchCollection(fusionModel.patches, match_original=True)
        fax.add_collection(p)

    # Save output
    fig.suptitle(title, fontsize=20, fontweight='bold')
    plt.show()
    fig.savefig(outName)
