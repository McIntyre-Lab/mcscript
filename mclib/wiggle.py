import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection


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


fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax2.set_ylim(-5, 20)
ax2.set_xlim(-5, 100)
ax2.set_axis_off()
patches = []
patches.append(Rectangle((0, 0), 10, 4, fill=True, color='black'))
patches.append(Rectangle((20, 0), 10, 4, fill=True, color='black'))
p = PatchCollection(patches, match_original=True)
ax2.add_collection(p)
ax2.add_line(Line2D((10, 20), (2,2),linewidth=2,color='black'))
plt.show()
fig.savefig('/home/jfear/Desktop/test.png')






