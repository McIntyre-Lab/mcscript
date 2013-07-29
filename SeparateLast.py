import os;
import sys;
from Bio import SeqIO


if len(sys.argv) != 3:
    print "Error, wrong number of arguments. \nUsage: SeparateLast fastq_file num_files"
else:
    input_file=sys.argv[1];
    num_files=int(sys.argv[2]); 

    fileArray=[];
    for numFile in range(0, num_files):
        fileArray.append(open(str(numFile)+"-split.fq", 'w+'))

    writeFlag = 0

    with open(input_file) as read_input_file:
        for record in SeqIO.parse(read_input_file, 'fastq'):
            SeqIO.write(record, fileArray[writeFlag], 'fastq')
            if writeFlag < num_files-1:
                writeFlag += 1
            else:
                writeFlag = 0

    for outfile in fileArray:
        outfile.close()
