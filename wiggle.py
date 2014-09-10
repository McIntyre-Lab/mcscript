import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection

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












def drawExon(exonTupleList, exonDict, exoncolors, introns):
	""" Draws an exon at the proper location. """
	# Rectangle (xpos (left), ypos (bottom), width, height, kwargs)
	patches = [Rectangle((0, -1), 1, 2, fill=False)] # initialize patches with a blank rectangle to avoid errors
	loc = 0 
	exonColor = 1

	for exon in exonTupleList:
		# width assumes exclusive nature for exons: that is, exon (20,22) is [20,22) and a width 1 rectangle will be drawn
		if introns: # if drawing introns
			start = exonDict[exon[0]] # make the start the location given by the exonDict
			width = exonDict[exon[1]-1] - exonDict[exon[0]] # make the width equal to the width between the locations
		else: # if drawing only exons
			start = loc # make the start at our current location
			width = exon[1] - exon[0] - 1 # make the width our distance between pairs
		# for exon colors
		if exonColor % 2 == 1: col = exoncolors[0] # if odd, use exoncolor1
		else: col = exoncolors[1] # if even, use exoncolor2
			
		patches.append(Rectangle((start, -1), width, 2, color=col)) # draw the exon
		loc += width + 1 # move the current location
		exonColor += 1
	return PatchCollection(patches, match_original=True)

def SetupPlot(plotSize, dimensions, title, chrom, codons):
	'''Set up the parameters for the graph.  Receives dimensions [xmin, xmax, ymin, ymax], alleleColor [color1, color2], title for the graph, and chromosome number.'''
	fig = plt.figure(figsize=(plotSize[0], plotSize[1]))# set window size to width, height 
	# add axes in rectangle left position, bottom position, width, height
	ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
	ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], sharex=ax1)
	ax1.axis(dimensions) # set up axis ranges
	ax2.set_ylim(-5, 5)
	# set up titles, labels, and ticks
	ax1.set_title(title)
	if codons: ax2.set_xlabel("Chromosome %s (codons from start of gene)"%chrom)
	else: ax2.set_xlabel("Chromosome %s (bp from start of gene)"%chrom)
	ax1.grid(True)		
	ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int((abs(x))))))# set ticks to absolute value
	# if we choose to graph codons, then divide all ticks by 3
	if codons: ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: round(x/3, 1)))
	ax2.set_yticks([]) # eliminate yticks from the 2nd axis
	for label in ax1.xaxis.get_ticklabels():
		label.set_fontsize(0)# eliminate the labels from our top plot without eliminating them from the bottom
	horizontalLine = Line2D([-10000000000, 10000000000], [0, 0], linewidth=1, color='black')# draw horizontal x-axis
	ax1.add_line(horizontalLine)

	return ax1, ax2, fig

if __name__ == '__main__':

    from mclib import gff as mcgff
    from mclib import wiggle as mcwiggle

    fname = '/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff'
    myGff = mcgff.FlyGff(fname)
    myGene = mcgff.FlyGene('InR', myGff)
    model = mcwiggle.GeneModel(myGene)












    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(211)
    ax1.set_ylim(0, 200)    # NOTE: set ylim equal to 0, max
    ax1.set_xlim(model.xLoc[0]-100, model.xLoc[1]+100) # NOTE: set xlim equal to start-100, end-100
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.set_ylim(min(model.yLoc), max(model.yLoc)+5)
    ax2.get_yaxis().set_visible(False)  # Hide y-axis on gene model plot
    ax2.get_xaxis().set_visible(False)  # Hide y-axis on gene model plot
    p = PatchCollection(model.patches, match_original=True)
    ax2.add_collection(p)
    #ax2.add_line(Line2D((10, 20), (2,2),linewidth=2,color='black'))
    plt.show()
    #fig.savefig('/home/jfear/Desktop/test.png')


