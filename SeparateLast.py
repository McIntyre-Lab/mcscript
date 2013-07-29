import os;
import sys;
if len(sys.argv) != 4:
	print "Error, wrong number of arguments. \nUsage: SeparateLast input_file num_files num_lines"
else:
	input_file=sys.argv[1];
	num_files=int(sys.argv[2]); 
	num_lines=int(sys.argv[3]);	
	
	fileArray=[];
	for numFile in range(0, num_files):
		fileArray.append(open(n+"-split.fq", 'w+'))
	readline=0
	writefile=0
 #read each line of the input file and write it to one of the output files. Write each four lines to the next output file, round robin style, until done. 
	with open(input_file) as read_input_file:
		for writeline in read_input_file.readlines():
			fileArray[writefile].write(writeline)
			readline=readline+1
			if readline%num_lines==0:
				writefile=(writefile+1)%num_lines
	for outfile in fileArray:
		outfile.close()
	
