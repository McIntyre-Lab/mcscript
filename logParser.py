
#
# DESCRIPTION: This script parses the given input bowtie and/or LAST files and creates a csv row of their data in the given output csv.
#
# AUTHOR: Chelsea Tymms

 
import sys, os.path
import argparse

def getOptions():
	"""Function to pull in arguments from the command line"""
	description="""This script takes an input fasta file of fusions and identifies all of the identical fusions."""
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-bowtie", "--bowtie_log_names", dest="bowtie", action='store', required=False, nargs = '*', help="bowtie log file names [Optional]")
	parser.add_argument("-last", "--last_log_names", dest="last", action='store', required=False, help="LAST log file names [Optional]")
	parser.add_argument("-treatment","--treatment_name",dest="treatment",action='store',required=True,nargs= '*', help="Treatment variables [Required]")
	parser.add_argument("-o","--output_file",dest="output",action='store',required=True,help="Output file name [Required]")
	args = parser.parse_args()
	if not args.bowtie and not args.last: #The user should give at least one bowtie or last log argument; otherwise the program does nothing
	    parser.error('No input logs given; add -bowtie or -last')
	return(args)


def main():
    args=getOptions()
    
    treatmentArray=args.treatment
    
    firstBowtieTot=0
    finalBowtieUnaln=0
    uniqAln=0
    
    #If the output file already exists, we will append to it. If it does not, we will open it and write its header.
    if os.path.isfile(args.output): #we will append
        outputFile=open(args.output,'ab')
             
    else: #write the header
        outputFile=open(args.output,'w')
        for i in range(1,len(treatmentArray)+1):
            outputFile.write('t_var_'+str(i)+',')
        if args.bowtie:
            for i in range(1,len(args.bowtie)+1):
                bowtieNum='bowtie'+str(i)
                outputFile.write(','.join(bowtieNum+'_'+n for n in ['tot','aln','unaln','ambig','per_uniq','per_aln'])+',')
        if args.last:
            outputFile.write(','.join(['last_uniq','last_ambig','last_per_uniq','last_per_aln'])+',')
        outputFile.write('per_uniq_aln'+'\n')
        
     
     
    outputFile.write(','.join(str(i) for i in treatmentArray)+',') 
    if args.bowtie:
        #Get some important counts from the first and the final bowtie logs
        proc,aln,unaln,ambig=parseBowtieLog(args.bowtie[0])
        firstBowtieTot=proc
        proc,aln,unaln,ambig=parseBowtieLog(args.bowtie[-1])
        finalBowtieUnaln=ambig+unaln
       
        #Get and write the counts for each Bowtie log
        for bowtieLog in args.bowtie:
            proc,aln,unaln,ambig=(parseBowtieLog(bowtieLog))
            perUniq,perAln=0,0
            if proc!=0:
                perUniq=float(aln)/proc * 100
                perAln=(float(aln)+ambig)/proc * 100
            uniqAln=uniqAln+aln
            outputFile.write(','.join(str(i) for i in [proc,aln,unaln,ambig,perUniq,perAln])+',')
            
    #Get and write the counts for the LAST log     
    if args.last:
        lastLog=args.last     
        ambig,uniq=(parseLastLog(lastLog))
        lastPerUniq,lastPerAln = 0,0
        if finalBowtieUnaln!=0:
            lastPerUniq=float(uniq)/finalBowtieUnaln * 100           
            lastPerAln=float(ambig+uniq)/finalBowtieUnaln * 100
        uniqAln=uniqAln+uniq
        outputFile.write(','.join(str(i) for i in [uniq,ambig,lastPerUniq,lastPerAln])+',')
     
    perUniqAln= perUniqAln=float(uniqAln)/firstBowtieTot * 100 if firstBowtieTot!=0 else 0
    outputFile.write(str(perUniqAln)+'\n')
    outputFile.close()
    
    
def parseBowtieLog(fileName):
    """Function to parse a bowtie log file"""
    if not os.path.isfile(fileName):
        print "WARNING: " +fileName+" does not exist."
        return 0,0,0,0 

    processed,aligned,unaligned,ambig=0,0,0,0
    with open(fileName,'rb') as bowtieLogFile:
        for line in bowtieLogFile.readlines():
            if 'reads processed' in line:
                processed=line.split(':')[1].strip()
            elif 'reads with at least one reported alignment' in line:
                aligned=line.split(':')[1].split(' ')[1]
            elif 'reads that failed to align' in line:
                unaligned=line.split(':')[1].split(' ')[1]
            elif 'reads with alignments suppressed' in line:
                ambig=line.split(':')[1].split(' ')[1]

    return int(processed),int(aligned),int(unaligned),int(ambig)
            
        
        
def parseLastLog(fileName): 
    """Function to parse a LAST log file"""
    if not os.path.isfile(fileName):
        print "WARNING: " +fileName+" does not exist."
        return 0,0 
    lastAmbig=0
    lastUniq=0
    with open(fileName,'rb') as lastLogFile:
        for line in lastLogFile.readlines():
            if "Ambiguously Aligned Reads" in line:
                lastAmbig=line.split(':')[1].strip()
            elif "Uniquely Aligned Reads" in line:   
                lastUniq=line.split(':')[1].strip()
    return int(lastAmbig),int(lastUniq)
 

if __name__ == '__main__':
    main()
        
