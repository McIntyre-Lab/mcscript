#!/usr/bin/env python
import csv, sys, argparse, os,itertools,operator,collections
csv.field_size_limit(1000000000)


# Parse command line arguments
parser = argparse.ArgumentParser(description='Counting total and unique fqs')
parser.add_argument('-i','--input',dest='fq', action='store', required=True, help='A list of fq file [Required]')
parser.add_argument('-t','--table',dest='table',action='store',required=False,help='Output table of intermediate size information [Optional]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for counts in csv format [Required]')
args = parser.parse_args()


### Open fq
with open(args.fq,'rb') as input_fq:
# Get the list of all of the sequences
# Sequences are every 4 lines starting with line 1
   seqs=itertools.islice(input_fq,1,None,4)
# Get the count for each sequence
   counts=collections.Counter(seqs)
   total_num=sum(counts.values())
   uniq_num=len(counts)

percent_uniq = float(uniq_num) / (total_num)
with open(args.out,'wb') as dataout:
    dataout.write('total # fq is '+str(total_num)+'\n# unique sequences is '+str(uniq_num)+'\npercent unique is '+str(percent_uniq))

if args.table:
# Sort by count in descending order
        sorted_counts=sorted(counts.iteritems(),key=operator.itemgetter(1),reverse=True)

   	with open(args.table,'wb') as tableOut:
   	    for item in sorted_counts:
                tableOut.write(str(item[1])+' '+str(item[0]).strip()+'\n')
 


