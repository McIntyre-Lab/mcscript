#!/usr/bin/env python

#	DESCRIPTION: This program uses the transcript FASTA file and extracts the junctions made by each pair of exons. 
#	The output is a bed file with all of these junctions. 
#
#	AUTHOR: Chelsea Tymms
#   EDITOR: Jonathan Poisson | poissonj@ufl.edu

# Build-in packages
import argparse  # Command line use

# Add-on packages
import gffutils
import itertools


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--input", dest="databaseFile", action='store', required=True, help="Input database")
    parser.add_argument("--output", dest="outputFile", action='store', required=True, help="Output file name")
    parser.add_argument("--size", dest="junctionSize", action='store', type=int, required=True, help="Size of Junction")

    args = parser.parse_args()
    return args


def createExonArray(exonList, mergedExons):
    """ Group overlapping exons together as an list within a list """
    exonArray = []

    for mergedExon in mergedExons:
        exonSubset = []

        for exon in exonList:
            if exon.start >= mergedExon.start and exon.stop <= mergedExon.stop:
                exonSubset.append(exon)

        exonArray.append(exonSubset)

    return exonArray


def createJunctionArray(exonArray):
    """ Using itertools create all possible combination between groups of exons """
    iterProduct = []

    for i in range(0, len(exonArray)):
        for j in range(i + 1, len(exonArray)):
            for k in itertools.product(exonArray[i], exonArray[j]):
                iterProduct.append(k)

    return iterProduct


def defineRegion(index, junctionArray):
    """ Instead of keeping the entire exon, only keep up to specified junciton size on either side of the junction """
    if len(junctionArray[index][0]) >= args.junctionSize + 1:
        exon1Start = junctionArray[index][0].stop - args.junctionSize
        exon1Stop  = junctionArray[index][0].stop
    else:
        exon1Start = junctionArray[index][0].start
        exon1Stop  = junctionArray[index][0].stop

    if len(junctionArray[index][1]) >= args.junctionSize + 1:
        exon2Start = junctionArray[index][1].start
        exon2Stop  = junctionArray[index][1].start + args.junctionSize
    else:
        exon2Start = junctionArray[index][1].start
        exon2Stop  = junctionArray[index][1].stop

    return (exon1Start, exon1Stop, exon2Start, exon2Stop)


def testOverlap(e1, e2):
    # Sanity check to make sure exons are the same
    if e1.chrom == e2.chrom:
        if e1.start <= e2.start and e1.end >= e2.start:  # They overlap
            return True
        else:
            return False
    else:
        print "Exons {0} and {1} are not on the same chromosome, something is wrong.".format(e1.id, e2.id)
        raise Exception


def main():
    # Get the database
    db_fn = args.databaseFile
    db = gffutils.FeatureDB(db_fn)
    genes = db.features_of_type('gene', order_by='start')

    # Open the output BED file
    with open(args.outputFile, 'wb') as outputFile:
        for gene in genes:

            # Get all exons
            exonList = list(db.children(gene, featuretype='exon', order_by='start'))

            # Skip genes that only have 1 exon, or genes that don't have an exon e.g. miRNA
            if len(exonList) == 0 or len(exonList) == 1:
                continue

            # Convert to BED 0-based coordinate system. Only have to change start.
            for exon in exonList:
                exon.start = exon.start - 1

            # Make all pairwise combos
            combos = itertools.combinations(exonList, 2)

            junctionArray = []
            for combo in combos:
                if not testOverlap(*combo):
                    junctionArray.append(combo)

            # Create BED file
            for i in range(0, len(junctionArray)):
                # Only keep up to 100bp on each side of the junction
                (exon1Start, exon1Stop, exon2Start, exon2Stop) = defineRegion(i, junctionArray)

                # Construct various parts of the junction BED file
                totalStart = exon1Start  # Take exon1 start
                totalStop  = exon2Stop  # Take exon2 stop
                chrom      = junctionArray[i][0].chrom
                strand     = junctionArray[i][0].strand
                score      = junctionArray[i][0].score
                color      = "255,0,0"
                name       = junctionArray[i][0].id + '|' + junctionArray[i][1].id
                num        = "2"
                lengths    = str(len(range(exon1Start, exon1Stop))) + ',' + str(len(range(exon2Start, exon2Stop)))
                starts     = str(exon1Start - totalStart) + ',' + str(exon2Start - totalStart)
                bedArray   = [chrom, totalStart, totalStop, name, score, strand, totalStart, totalStop, color, num,
                              lengths, starts]

                # Output to BED file
                outputFile.write("\t".join(str(i) for i in bedArray) + "\n")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
