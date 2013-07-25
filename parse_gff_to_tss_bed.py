#!/usr/local/bin/python2.7

import argparse, logging

def getOptions():
    """ Function to pull in arguments """

    description="""This script takes a GFF file from FlyBase and creates a BED file with a region that is up and downstream of the expected transcription start site (TSS). """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-f", "--gff", dest="gname", action='store', required=True, help="FlyBase GFF file [Required]",metavar="GFF_FILE")
    parser.add_argument("-b", "--bound", dest="bound", action='store', type=int, required=False, default=2000, help="The size of the region you want as a boundry [Default: 2000]", metavar="BOUNDRY")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    parser.add_argument("-o", "--out", dest="oname", action='store', required=True, help="Out File", metavar="OUT_FILE")
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def neg_chk(tss_start):
    """ There is a chance where the tss_start - 2kb would be less than 0. I need to correct for this. """

    if tss_start < 0:
        tss_start = 0
    return(tss_start)

def calc_tss_bound(record,bound):
    """ Calculates the up and downstream region from the estimated transcription start site (tss_pos) """

    chrom = record[0]
    start = int(record[3])
    end = int(record[4])
    strand = record[6]
    if strand == '+':
        tss_pos = start
        tss_start = tss_pos - bound
        tss_end = tss_pos + bound
        tss_start = neg_chk(tss_start)
    elif strand == '-':
        tss_pos = end
        tss_start = tss_pos - bound
        tss_end = tss_pos + bound
        tss_start = neg_chk(tss_start)
    else:
        logging.info('This record had no strand information:')
        logging.info(record)
        return()

    return(chrom, tss_start, tss_end)

def parse_gff_mrna(fname,bound):
    logging.info('Parsing GFF File')
    gdict = {}
    with open(fname, 'r') as GFF:
        for row in GFF:
            if not row.startswith("#"):
                record = row.strip().split('\t')
                if len(record) >= 7:
                    if record[1] == 'FlyBase' and record[2] == 'mRNA':
                        try:                                                            # There are ~6 instances where strand('+/-') is not known so this script skips these.
                            chrom, tss_start, tss_end = calc_tss_bound(record,bound)    # Calculate the up and down stream region for the TSS start site
                            loc = chrom + ':' + str(tss_start) + '-' + str(tss_end)     # I am creating a dictionary key so if there are any mRNA with the same starting location they should get filtered out here.
                            if not gdict.has_key(loc):
                                gdict[loc] = {}
                                gdict[loc]['chrom'] = chrom
                                gdict[loc]['start'] = tss_start
                                gdict[loc]['end'] = tss_end
                        except:
                            pass

    return(gdict)

def write_bed(gdict, oname):
    logging.info('Writing BED File')
    with open(oname, 'w') as OUT:
        for key in gdict:
            data = [ gdict[key]['chrom'], gdict[key]['start'], gdict[key]['end'] ]
            OUT.write('\t'.join(str(x) for x in data) + '\n')

def main():

    args = getOptions()

    if args.log:                                  # Turn on Logging if option -g was given
        setLogger(args.log,logging.INFO)

    gdict = parse_gff_mrna(args.gname,args.bound)
    write_bed(gdict ,args.oname)

if __name__ == '__main__':
    main()
