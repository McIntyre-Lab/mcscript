#!/usr/bin/env python
#Standar Libraries
import os
import argparse
import numpy as np
import re
import logging

#AddOn Libraries
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getOptions():
	"""Function to pull in arguments"""
	parser = argparse.ArgumentParser(description="Takes a SAM file [either SE or PE] and splits its reads into different categories")
	parser.add_argument("-s", dest="sam", action='store', required=True, help="Name of the SAM file [required]")
	parser.add_argument("-fq1", dest="fq1", action='store', required=True, help="Name of the FastQ1 file [required]")
	parser.add_argument("-fq2", dest="fq2", action='store', required=False, default=False, help="Name of the FastQ2 file")
	parser.add_argument("-summname", dest="summname", action='store', required=False, help="Name in the summary file")
	parser.add_argument("--outdir", dest="odir" , action='store', required=True, help="Out directory path")
	parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file")
	#parser.add_argument("-o", dest="oname", action='store', required=True)
	args=parser.parse_args()
	return (args)

def readFastQ(fastq_path):
	if args.log:
		logger.info("Reading fastq'%s' " % (fastq_path))
	fastq_file = open(fastq_path,'r')
	fastq_generator=FastqGeneralIterator(fastq_file)
	readDict = {re.sub('/[1-2]','',header).split(' ')[0]:(seq,qual) for header,seq,qual in fastq_generator}
	if args.log:
		logger.info("Finished Reading fastq '%s' " % (fastq_path))   
	return (readDict)

def writeOutput (headList,readDict,OUTFILE):
	if args.log:
		logger.info("Printing fastq results '%s'" %(OUTFILE))   
	for head in headList:
		OUTFILE.write ('\n'.join(['@'+head,readDict[head][0],'+',readDict[head][1],'']))
	if args.log:
		logger.info("Finished fastq results")

def SplitSAMPE (fname,odir,summname):
	"""Function to split all the reads in PE SAM files"""
    
	#Setting flags
	flags_opositestrand = ["81","97","145","161"]
	flags_grayones = ["321","323","329","337","339","353","355","369","371","377","385","387","393","401","403","417","419","433","435","441"]
	flags_unmappedread1 = ["69","101","177"]
	flags_unmappedread2 = ["133","165","181"]
	flags_bothunmapped1 = ["77"]
	flags_bothunmapped2 = ["141"]
	flags_mapped1 = ["65","73","83","89","99","113","121"]
	flags_mapped2 = ["153","185","137","147","163","129","117"]

	#Setting counters
	total = 0
	counter_opositestrand = 0
	counter_grayones = 0
	counter_unmappedread1 = 0
	counter_unmappedread2 = 0
	counter_bothunmapped1 = 0
	counter_bothunmapped2 = 0
	counter_mapped1 = 0
	counter_mapped2 = 0
	counter_ambiguous1 = 0
	counter_ambiguous2 = 0

	#Lists for unmapped and ambiguous reads
	unmappedread1 = []
	unmappedread2 = []
	bothunmapped1 = []
	bothunmapped2 = []
	ambiguous1 = []
	ambiguous2 = []

	#Filename
	bname=os.path.basename(fname)
	name = os.path.splitext(bname)[0]

	#Open SAM file and output files in SAM format.
	SAM = open(fname,'r')
	GRAY = open(os.path.join(odir,name+'_gray.sam'),'w')
	OPOSITE = open(os.path.join(odir,name+'_opposite.sam'),'w')
	UNRECOGNIZED = open(os.path.join(odir,name+'_unrecognized.sam'),'w')
	MAPPED = open(os.path.join(odir,name+'_mapped.sam'),'w')
	AMBIGUOUS = open(os.path.join(odir,name+'_ambiguous.sam'),'w')

	#Open Sumary file
	SUMMARY = open(os.path.join(odir,name+'_summary.csv'),'w')

	#Reading line by line SAM file (except headers)
	for line in SAM:
		if line.startswith('@'):continue	
		line=line.replace('\n','')
		line=line.split('\t')

		#Getting unmapped reads
		if line[1] in flags_unmappedread1:
			unmappedread1.append(line[0])
			counter_unmappedread1 += 1
			total += 1
		elif line[1] in flags_unmappedread2:
			unmappedread2.append(line[0])
			counter_unmappedread2 += 1
			total += 1
		elif line[1] in flags_bothunmapped1:
			bothunmapped1.append(line[0])
			counter_bothunmapped1 += 1
			total += 1	
		elif line[1] in flags_bothunmapped2:
			bothunmapped2.append(line[0])
			counter_bothunmapped2 += 1	
			total += 1
		#Getting & printing "gray" reads
		elif line[1] in flags_grayones:
			print >> GRAY,'\t'.join(line)
			counter_grayones += 1
			total += 1
		#Getting & printing "OPOSITE" reads
		elif line[1] in flags_opositestrand:
			print >> OPOSITE,'\t'.join(line)
			counter_opositestrand += 1
			total += 1
		#Getting & printing AMBIGUOUS reads, those who are not ambiguous are store as mapped reads
		elif line[1] in flags_mapped1:
			if line[-1].startswith('SA:'):
				if int(line[-2].replace('XS:i:','')) - int(line[-3].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous1.append(line[0])
					counter_ambiguous1 += 1
					total += 1
				else:
					#mappedread1.append(line[0])
					print >> MAPPED,'\t'.join(line)
					counter_mapped1 += 1
					total += 1
			else:
				if int(line[-1].replace('XS:i:','')) - int(line[-2].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous1.append(line[0])
					counter_ambiguous1 += 1
					total += 1
				else:
					#mappedread1.append(line[0])
					print >> MAPPED,'\t'.join(line)
					counter_mapped1 += 1
					total += 1

		elif line[1] in flags_mapped2:
			if line[-1].startswith('SA:'):
				if int(line[-2].replace('XS:i:','')) - int(line[-3].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous2.append(line[0])
					counter_ambiguous2 += 1
					total += 1
				else:
					#mappedread2.append(line[0])
					print >> MAPPED,'\t'.join(line)
					counter_mapped2 += 1
					total += 1
			else:
				if int(line[-1].replace('XS:i:','')) - int(line[-2].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous2.append(line[0])
					counter_ambiguous2 += 1
					total += 1
				else:
					#mappedread2.append(line[0])
					print >> MAPPED,'\t'.join(line)
					counter_mapped2 += 1
					total += 1
		#If not in the previous categories then unknown
		else:
			print "Warning: "+line[0]+" key is not recognized"
			print >> UNRECOGNIZED,'\t'.join(line)
            
	if args.log:
		logger.info("Finished reading of '%s' " % (fname))
		logger.info("Printing summary")

	#Print summary
	count_names = ["name","total_reads","counter_oposite_strand_read","counter_grayones","counter_unmapped_read1","counter_unmapped_read2","counter_both_unmapped_read1","counter_both_unmapped_read2","counter_mapped_read1","counter_mapped_read2","counter_ambiguous_read1","counter_ambiguous_read2"] 
	count_values = [summname,total,counter_opositestrand,counter_grayones,counter_unmappedread1,counter_unmappedread2,counter_bothunmapped1,counter_bothunmapped2,counter_mapped1,counter_mapped2,counter_ambiguous1,counter_ambiguous2]
	count_values = map(str,count_values)
	print >> SUMMARY,','.join(count_names)
	print >> SUMMARY,','.join(count_values)


	#Clossing all files
	SAM.close()
	GRAY.close()
	OPOSITE.close()
	UNRECOGNIZED.close()
	MAPPED.close()
	SUMMARY.close()
	AMBIGUOUS.close()
    
	if args.log:
		logger.info("Finished printing of SAM files")
        
	#return(unmappedread1,unmappedread2)
	return(unmappedread1,unmappedread2,bothunmapped1,bothunmapped2,ambiguous1,ambiguous2)

def SplitSAMSE (fname,odir,summname):
	"""Function to split all the reads in PE SAM files"""

	#Setting flags
	flags_opositestrand = ["16"]
	flags_unmappedreads = ["4"]
	flags_mapped = ["0"]
	
	#Setting counters
	total = 0
	counter_opositestrand = 0
	counter_unmappedread = 0
	counter_mapped = 0
	counter_ambiguous = 0

	#Lists for mapped and ambiguous reads
	unmappedread = []
	ambiguous = []

	#Filename
	bname=os.path.basename(fname)
	name = os.path.splitext(bname)[0]

	#Open SAM file and output files in SAM format.
	SAM = open(fname,'r')
	OPOSITE = open(os.path.join(odir,name+'_oposite.sam'),'w')
	MAPPED = open(os.path.join(odir,name+'_mapped.sam'),'w')
	AMBIGUOUS = open(os.path.join(odir,name+'_ambiguous.sam'),'w')

	#Open Sumary file
	SUMARY = open(os.path.join(odir,name+'_summary.csv'),'w')

	#Reading line by line SAM file (except headers)
	for line in SAM:
		if line.startswith('@'):continue	
		line=line.replace('\n','')
		line=line.split('\t')
		#Getting unmapped reads
		if line[1] in flags_unmappedreads:
			unmappedread.append(line[0])
			counter_unmappedread += 1
			total += 1
		#Getting & printing "OPOSITE" reads
		elif line[1] in flags_opositestrand:
			print >> OPOSITE,'\t'.join(line)
			counter_opositestrand += 1
			total += 1
		#Getting & printing AMBIGUOUS reads, those who are not ambiguous are store as mapped reads
		elif line[1] in flags_mapped:
			if line[-1].startswith('SA:'):
				if int(line[-2].replace('XS:i:','')) - int(line[-3].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous.append(line[0])
					counter_ambiguous += 1
					total += 1
				else:
					print >> MAPPED,'\t'.join(line)
					counter_mapped += 1
					total += 1
			else:
				if int(line[-1].replace('XS:i:','')) - int(line[-2].replace('AS:i:','')) == 0:
					print >> AMBIGUOUS,'\t'.join(line)
					ambiguous.append(line[0])
					counter_ambiguous += 1
					total += 1
				else:
					#mappedread1.append(line[0])
					print >> MAPPED,'\t'.join(line)
					counter_mapped += 1
					total += 1
		#If not in the previous categories then unknown
		else:
			print "Warning: "+line[0]+" key "+line[1]+" is not recognized"

	if args.log:
		logger.info("Finished reading of '%s' " % (fname))
		logger.info("Printing summary")

	#Print summary
	count_names = ["name","total_reads","count_mapped_read_oposite_strand","count_unmapped_read","count_mapped_read","count_ambiguous_read"] 
	count_values = [summname,total,counter_opositestrand,counter_unmappedread,counter_mapped,counter_ambiguous]
	count_values = map(str,count_values)
	print >> SUMARY,','.join(count_names)
	print >> SUMARY,','.join(count_values)

	#Clossing all files
	SAM.close()
	OPOSITE.close()
	MAPPED.close()
	SUMARY.close()
	AMBIGUOUS.close()
    
	if args.log:
		logger.info("Finished printing of SAM files")
        
	#return(unmappedread1,unmappedread2)
	return(unmappedread,ambiguous)

def main(args):
	"""Here we call all other functions"""

	#Paths
	odir = os.path.abspath(args.odir)
	bname = os.path.basename(args.sam)
	name = os.path.splitext(bname)[0]
	fname = args.sam

	if not args.summname:
		summname = name
	else:
		summname = args.summname

	if args.log:
		logger.info("Reading '%s' " % (fname))

	if args.fq2:
		unmapped1,unmapped2,bothunmapped1,bothunmapped2,ambiguous1,ambiguous2 = SplitSAMPE(fname,odir,summname)
	        
		#Open PE output files
		UNMAPPED1 = open (os.path.join(odir, name + '_unmapped1.fq'),'w')
		AMBIGUOUS1 = open (os.path.join(odir, name + '_ambiguous1.fq'),'w')
		BOTHUNMAPPED1 = open (os.path.join(odir,name + '_both_unmapped1.fq'),'w')
		UNMAPPED2 = open (os.path.join(odir, name + '_unmapped2.fq'),'w')
		AMBIGUOUS2 = open (os.path.join(odir, name + '_ambiguous2.fq'),'w')
		BOTHUNMAPPED2 = open (os.path.join(odir,name + '_both_unmapped2.fq'),'w')

		#Print unMapped1, bothinmapped1 and ambiguous1
		fastQ1Dict = readFastQ(args.fq1)
		writeOutput (unmapped1,fastQ1Dict,UNMAPPED1)
		writeOutput (ambiguous1,fastQ1Dict,AMBIGUOUS1)
		writeOutput (bothunmapped1,fastQ1Dict,BOTHUNMAPPED1)
		del fastQ1Dict

		#Print unMapped1, bothinmapped2 and ambiguous2
		fastQ2Dict = readFastQ(args.fq2)
		writeOutput (unmapped2,fastQ2Dict,UNMAPPED2)
		writeOutput (ambiguous2,fastQ2Dict,AMBIGUOUS2)
		writeOutput (bothunmapped2,fastQ2Dict,BOTHUNMAPPED2)
		del fastQ2Dict

		#Close PE output files
		UNMAPPED1.close()
		AMBIGUOUS1.close()
		BOTHUNMAPPED1.close()
		UNMAPPED2.close()
		AMBIGUOUS2.close()
		BOTHUNMAPPED2.close()

	else:
		unmapped,ambiguous = SplitSAMSE(fname,odir,summname)

		#Open SE output files
		UNMAPPED = open (os.path.join(odir, name + '_unmapped.fq'),'w')
		AMBIGUOUS = open (os.path.join(odir, name + '_ambiguous.fq'),'w')

		#Print unMapped, and ambiguous
		fastQDict = readFastQ(args.fq1)
		writeOutput (unmapped,fastQDict,UNMAPPED)
		writeOutput (ambiguous,fastQDict,AMBIGUOUS)

		#Close SE output files
		UNMAPPED.close()
		AMBIGUOUS.close()

if __name__=='__main__':
	args = getOptions()
	if args.log:
		logging.basicConfig(filename=(os.path.abspath(args.log)),level=logging.DEBUG)
		logger = logging.getLogger()
	main(args)
	if args.log:
		logger.info("Script complete.")

