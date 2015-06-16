#!/usr/bin/env python

#	DESCRIPTION: This program uses the transcript FASTA file and extracts the junctions made by each pair of exons. 
#	The output is a bed file with all of these junctions. 
#
#	AUTHOR: Chelsea Tymms


import gffutils
import itertools

def start_pos(x):
    """ Little function to let me sort gffutils objects by start position """
    return x.start

def createExonArray(exonList,mergedExons):
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

    for i in range(0,len(exonArray)):
        for j in range (i+1,len(exonArray)):
            for k in itertools.product(exonArray[i],exonArray[j]):
                iterProduct.append(k)

    return iterProduct

def defineRegion(index,junctionArray):
    """ Instead of keeping the entire exon, only keep up to 100bp on either side of the junction """
    if len(junctionArray[index][0]) >= 101:
        exon1Start = junctionArray[index][0].stop - 100
        exon1Stop = junctionArray[index][0].stop
    else:
        exon1Start = junctionArray[index][0].start
        exon1Stop = junctionArray[index][0].stop

    if len(junctionArray[index][1]) >= 101:
        exon2Start = junctionArray[index][1].start 
        exon2Stop = junctionArray[index][1].start + 100
    else:
        exon2Start = junctionArray[index][1].start 
        exon2Stop = junctionArray[index][1].stop

    return (exon1Start,exon1Stop,exon2Start,exon2Stop)


def main():
    # Get the database
    db_fn='../flybase_files/dmel-all-r5.51.gff.db'
    db=gffutils.FeatureDB(db_fn)
    genes=db.features_of_type('gene')

    # Open the output BED file
    with open('../output/test_junctions.bed','wb') as outputFile:
        for gene in genes:

            # Get all exons
            exonList = list(db.children(gene,featuretype='exon'))
            
            # Skip genes that only have 1 exon, or genes that don't have an exon e.g. miRNA
            if len(exonList) == 0 or len(exonList) == 1:
                continue

            # Convert to BED 0-based coordinate system. Only have to change start.
            for exon in exonList:
                exon.start = exon.start - 1 

            # Create list of overlapping exons
            # Note: mod(mdg4) has exons on both strands, so ignore_strand needs
            # to be set for it
            mergedExons = list(db.merge_features(exonList,ignore_strand=True))

            # Sort merged exons by start postiion
            mergedExons.sort(key=start_pos)

            # Using the merged exons to identify overlapping exons create an array where
            # overlapping exons are grouped together
            exonArray = createExonArray(exonList,mergedExons)

            # Now create all possible junctions 
            junctionArray = createJunctionArray(exonArray)

            # Create BED file
            for i in range(0,len(junctionArray)):

                # Only keep up to 100bp on each side of the junction
                (exon1Start, exon1Stop, exon2Start, exon2Stop) = defineRegion(i,junctionArray)

                # Construct various parts of the junction BED file
                totalStart=exon1Start # Take exon1 start
                totalStop=exon2Stop  # Take exon2 stop
                chrom=junctionArray[i][0].chrom
                strand=junctionArray[i][0].strand
                score=junctionArray[i][0].score
                color="255,0,0"
                name=junctionArray[i][0].id + '|' + junctionArray[i][1].id
                num="2"
                lengths=str(len(range(exon1Start,exon1Stop)))+','+str(len(range(exon2Start,exon2Stop)))
                starts=str(exon1Start-totalStart)+','+str(exon2Start-totalStart)
                bedArray=[chrom,totalStart,totalStop,name,score,strand,totalStart,totalStop,color,num,lengths,starts]

                # Output to BED file
                outputFile.write("\t".join(str(i) for i in bedArray)+"\n")
                        
if __name__ == '__main__':
    main()

