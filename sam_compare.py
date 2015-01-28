#!/usr/bin/env python

from __future__ import print_function

## @package sam_compare
#  This script compares two SAM files (sequence alignment maps) and
#  generates counts for similarly mapped reads.
#  Translated to python from an existing perl script.
#  @author Oleksandr Moskalenko <om@hpc.ufl.edu>
#  @version 0.1
#  Released: 2013-01-29
#
#  H. del Risco <hdelrisco@ufl.edu>
#  2014-03-21
#
#  Initial release modifications (HdR):
#  Made significant changes to the code and optimized where relevant
#  Used different approach for loading and processing data
#  Modified log handling and message display and added time performance tracking display

# Changed from sending data to stderr to file/stdout option

import sys, os, time, datetime, locale
import csv, cStringIO, logging, pprint, re, mmap, argparse, subprocess
import shlex, gc, getpass, random
from Bio import SeqIO

# This code is the translation, with modifications, of an existing perl script to python
#
__author__  = "Oleksandr Moskalenko <om@hpc.ufl.edu> and Hector del Risco <hdelrisco@ufl.edu>"
__version__ = '1.0.0'

locale.setlocale(locale.LC_ALL, '')
logger = logging.getLogger()

# parse command line options
def parse_options():
    ## Get Options and Arguments
    #  @brief Generate command-line switches and help, process options and arguments
    #  @param None
    #  @return options, arguments
    parser = argparse.ArgumentParser(description="Merge 2 sam files and count variant combinations.", epilog="""This tool keeps track of the
            number of alignments, whether the alignment is exact, and chromosome (RNAME)
            and position of the alignment. Only alignments in the top 'strata'
            are counted, meaning if there is 1 exact match and 3 inexact
            matches, the alignment count for that read will be 1. Strata is
            determined by subtracting the edit distance from the number of read
            nucleotides aligned. In the event of an exact match, this number is
            equal to the read length. Once both alignments processed, reads are
            compared one at a time and counted depending on the relationship
            between the two alignments.\n""")
    parser.add_argument("-l", "--length", dest="length", required=True, help="Read length")
    parser.add_argument("-q", "--fastq", dest="fastq", required=True, help="Source FastQ data file name")
    parser.add_argument("-A", "--sama", dest="sama", required=True, help="First SAM file (A)", metavar="SAM_A")
    parser.add_argument("-B", "--samb", dest="samb", required=True, help="Second SAM file (B)", metavar="SAM_B")
    parser.add_argument("-f", "--fusion", dest="fusion", required=True, help="Fusion file name, TSV or BED")
    parser.add_argument("-c", "--counts", dest="counts", help="Output file name for the counts - CSV file. Defaults to counts_username_date_time_rand.csv. To output to console, set to -o stdout")
    parser.add_argument("-t", "--totals", dest="totals", help="Output file name for the totals - text file. Defaults to totals_username_date_time_randint.txt. To output to console, set to -t stdout")
    parser.add_argument("-n", "--nofqids", dest="nofqids", default=False, action='store_true', help="Do not check SAM reads QNAME against the fastq sequence ids. Saves time if already known to be good. Must still set -q option.")
    parser.add_argument("-g", "--log", dest="log", help="Log file name. Defaults to samcompare_username_date_time_rand.log. To output to console, set to -g stdout. There is normally no stdout unless specifically requested via options.")
    parser.add_argument("-d", "--debug", dest="debug", default=False, action='store_true', help="Debugging mode (verbose). Includes elapsed time display for performance tracking.")
    opts = parser.parse_args()
    return opts

# set log output to either go to file or sys.stdout
# error level messages, and above, will be sent to sys.stderr
def setLogger(logfilename, level):
    global logger

    console = True if logfilename.lower() == "stdout" else False
    loghandler = logging.StreamHandler(stream=sys.stdout) if console else logging.FileHandler(logfilename)
    errhandler = logging.StreamHandler(stream=sys.stderr)
    logger.setLevel(level)
    loghandler.setLevel(level)
    errhandler.setLevel(logging.ERROR)

    logfmt = logging.Formatter('%(levelname)s - %(message)s') if console else logging.Formatter('%(levelname)s - %(message)s')
    errfmt = logging.Formatter('%(levelname)s - %(message)s')
    loghandler.setFormatter(logfmt)
    errhandler.setFormatter(errfmt)
    logger.addHandler(loghandler)
    logger.addHandler(errhandler)

# generate a temporary file name to use for all relevant files by forming 'filename_{base}.ext'
# will use for temporary files, output files, and log files if names not provides as arguments
def get_temp_basefilename():
    now = datetime.datetime.now()
    strdt = "{0:04}{1:02}{2:02}_{3:02}{4:02}{5:02}".format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    filename = "{0}_{1}_{2}".format(getpass.getuser(), strdt, random.randint(1,9999))
    return filename

# get ids from fastq file, use sed/cut for performance - avoids loading huge data file into memory
def get_fastq_ids(infilename, basefilename):
    fqids = {}

    outfilename = "fastqids_{0}".format(basefilename)
    cmd = "sed -n '1~4p' {0} | cut -f1 | cut -c2- > {1}".format(infilename, outfilename)
    logger.debug("FQIDs cmd: {0}".format(cmd))
    args = shlex.split(cmd)
    try:
        tstart = time.time()
        logger.debug("Loading FASTQ data from file...")
        # the check_call func has an issue, I think is due to the parsing of the args (using pipe, etc.) - figure it out later!!!
        #result = subprocess.check_call(args, shell=True)
        #if result == 0:
        output = subprocess.check_output(cmd, shell=True)
        if output == "":
            with open(outfilename) as f:
                ids = f.read().splitlines()
            logger.debug("Time to load FASTQ file into memory: {0:n} seconds, {1:n} lines".format(time.time() - tstart, len(ids)))
            if ids:
                logger.debug("Creating FASTQ dictionary...")
                tstart = time.time()
                fqids = dict((x,None) for x in ids)
                logger.debug("Time to create FASTQ ids dictionary: {0:n} seconds, {1:n} entries".format(time.time() - tstart, len(fqids)))
        else:
            logger.error("Unable to get FASTQ ids from file '{0}':\n{1}".format(infilename, output))
    except Exception as e:
        logger.error("Unable to get FASTQ ids from file '{0}':\n{1}".format(infilename, e))
        fqids = {}

    # remove temporary file, may not exist
    try:
        os.unlink(outfilename)
    except Exception as e:
        pass
    return fqids

# load sam file data into dictionary
def process_sam_file(samid, filename, fqids, read_length):
    samreads = {}
    samdata = []

    try:
        # load data from file into memory
        logger.debug("Loading SAM{0} data from file '{1}' ...".format(samid, os.path.basename(filename)))
        tstart = time.time()
        with open(filename) as f:
            samdata = f.read().splitlines()
        logger.debug("Time to load SAM{0} data from file: {1:n} seconds, {2:n} lines".format(samid, time.time() - tstart, len(samdata)))

        # process data records
        try:
            logger.debug("Processing SAM{0} data records...".format(samid))
            tstart = time.time()
            samreads = get_sam_reads(samdata, fqids, read_length, samid)
            logger.debug("Time to process SAM{0} data records: {1:n} seconds, {2:n} reads".format(samid, time.time() - tstart, len(samreads)))
        except Exception as e:
            logger.error("Unable to process SAM{0} data: {1}".format(samid, e))
            samreads = {}
    except Exception as e:
        logger.error("Unable to load SAM{0} data from file '{1}':\n{2}".format(samid, filename, e))

    return samreads

"""
SAM file format for reference
Col Field Type   Regexp/Range             Brief description
--- ----- ------ -----------------------  -----------------------------
  1 QNAME String [!-?A-~]{1,255}          Query template NAME
  2 FLAG  Int    [0,2^16-1]               bitwise FLAG
  3 RNAME String \*|[!-()+-<>-~][!-~]*    Reference sequence NAME
  4 POS   Int    [0,2^31-1]               1-based leftmost mapping POSition
  5 MAPQ  Int    [0,2^8-1]                MAPping Quality
  6 CIGAR String \*|([0-9]+[MIDNSHPX=])+  CIGAR string
  7 RNEXT String \*|=|[!-()+-<>-~][!-~]*  Ref. name of the mate/next read
  8 PNEXT Int    [0,2^31-1]               Position of the mate/next read
  9 TLEN  Int    [-2^31+1,2^31-1]         observed Template LENgth
 10 SEQ   String \*|[A-Za-z=.]+           segment SEQuence
 11 QUAL  String [!-~]+                   ASCII of Phred-scaled base QUALity+33

Optional fields:
 NM i Edit distance to the reference, including ambiguous bases but excluding clipping (NM:i:1)
"""

# process sam file data records
REQUIRED_SAMFIELDS = 11
def get_sam_reads(samdata, fqids, read_length, samid):
    total = 0
    nofastq = 0
    processed = 0
    dupqname = 0
    reads = {}
    re_mmc = re.compile('(\d+)M')

    # process all sam file input lines
    #print("    Count: 0", end="")
    #sys.stdout.flush()
    for line in samdata:
        #if total and (total % 1000000) == 0:
        #    print(",{0}M".format(int(total/1000000)), end="")
        #    sys.stdout.flush()
        total += 1

        # all file lines were loaded, only process data record lines
        if not line.startswith('@'):
          record = line.strip().split("\t")
          if record and len(record) >= REQUIRED_SAMFIELDS:
            name = record[0]
            if not fqids or name in fqids:
                # get reference name
                rname = record[2]

                # process this SAM record
                processed += 1

                # populate read record with default values
                if name in reads:
                    # "Reads/segments having identical QNAME are regarded to come from the same template"
                    dupqname += 1
                else:
                    # [0] RNAME, [1] count, [2] exact_match, [3] edit_distance (opts.length)
                    reads[name] = [rname, 0, 0, read_length]

                # get alignment matches 
                match_mismatch_count = 0
                for mmc in re_mmc.findall(record[5]):
                    match_mismatch_count += int(mmc)

                # get edit distance - mismatches, insertions, deletions, plus the difference between aligned length and total length
                edit_distance = 0
                try:
                    # don't waste time looking for NM in required fields and only deal with first find
                    findnm = (field for field in record[REQUIRED_SAMFIELDS:] if field.startswith("NM:i:"))
                    for nm in findnm:
                        edit_distance = int(nm.split(':')[2])
                except Exception as e:
                    pass # No 'NM:i:x' record found
                edit_distance += (len(record[9]) - match_mismatch_count)

                # check if exact match
                if edit_distance == 0:
                    # check if exact match flag already set and update count if so
                    if reads[name][2] == 1:
                        reads[name][1] += 1
                    else:
                        # set exact match flag and reset count
                        reads[name][1:] = [1, 1, 0]
                else:
                    # not an exact match, check if we got closer than previous (read_length initially)
                    if edit_distance < reads[name][3]:
                        # we got closer to a match, set new edit distance and reset count
                        reads[name][3] = edit_distance
                        reads[name][1] = 1
                    elif edit_distance == reads[name][3]:
                        # same edit distance, inc count
                        reads[name][1] += 1
            else:
                nofastq += 1

    # display debug stats
    logger.debug("Processed a total of {0:n} records out of {1:n} input lines".format(processed, total))
    if nofastq:
        logger.warning("Found {0:n} SAM{1} record(s) with no matching sequence in FASTQ file".format(nofastq, samid))
    if dupqname:
        logger.warning("Found {0:n} SAM{1} records with same QNAME (read/segment from same template)".format(dupqname, samid))
    return (reads)

# load fusion data into dictionary
def get_fusions(filename):
    fusions = {}

    try:
        with open(filename, "r") as f:
            lines = f.read().splitlines()
        for line in lines:
            fields = line.strip().split("\t")
            fusions[fields[0]] = None
        f.close()
    except Exception as e:
        fusions = {}
        logger.error("Unable to process fusion file '{0}':\n{1}".format(filename, e))
        if f:
            f.close()

    return fusions

# process data to get counts and totals
def process_read_counts(fusions, sama_reads, samb_reads):
    counts = {}
    totals = []

    # Totals
    a_single_exact = 0
    a_single_inexact = 0
    a_multi_exact = 0
    a_multi_inexact = 0
    b_single_exact = 0
    b_single_inexact = 0
    b_multi_exact = 0
    b_multi_inexact = 0
    both_unaligned = 0
    both_single_exact_same = 0
    both_single_exact_diff = 0
    both_single_inexact_same = 0
    both_single_inexact_diff = 0
    both_inexact_diff_equal = 0
    both_inexact_diff_a_better = 0
    both_inexact_diff_b_better = 0
    both_multi_exact = 0
    both_multi_inexact = 0
    a_se_b_si = 0
    a_si_b_se = 0
    a_se_b_me = 0
    a_me_b_se = 0
    a_se_b_mi = 0
    a_mi_b_se = 0
    a_si_b_me = 0
    a_me_b_si = 0
    a_si_b_mi = 0
    a_mi_b_si = 0
    a_me_b_mi = 0
    a_mi_b_me = 0
    total_count = 0

    # create template for count rows
    count_template = {'fusion':"", 'b_single_exact':0, 'b_single_inexact':0,
            'a_single_exact':0, 'a_single_inexact':0,
            'both_single_exact_same':0, 'both_single_exact_diff':0,
            'both_single_inexact_same':0, 'both_single_inexact_diff':0,
            'both_inexact_diff_equal':0, 'both_inexact_diff_a_better':0,
            'both_inexact_diff_b_better':0, 'a_exact_b_inexact':0, 'b_exact_a_inexact':0}

    # create count rows based on fusion name, add fusion name to each row for later display
    tstart = time.time()
    for fusion in fusions:
        if fusion not in counts:
            counts[fusion] = dict(count_template)
            counts[fusion]['fusion'] = fusion
    logger.debug("Time to create count rows: {0:n} seconds, {1:n} rows".format(time.time() - tstart, len(counts)))

    # get all the unique, read ids, QNAMEs for sama and samb
    tstart = time.time()
    readkeys_sama = set(sama_reads.keys())
    readkeys_samb = set(samb_reads.keys())
    readkeys = readkeys_sama.union(readkeys_samb)
    logger.debug("Readkey sets - A: {0:n}, B: {1:n}, Union: {2:n}".format(len(readkeys_sama), len(readkeys_samb), len(readkeys)))
    logger.debug("Time to create readkey sets: {0:n} seconds".format(time.time() - tstart))

    # process all reads
    cnt = 0
    addreadid = 0
    addedids = ""
    tstart = time.time()
    #print("    Count: 0", end="")
    #sys.stdout.flush()
    for read in readkeys:
        #if cnt and (cnt % 1000000) == 0:
        #    print(",{0}M".format(int(cnt/1000000)), end="")
        #    sys.stdout.flush()
        cnt += 1

        # setup a and b values accordingly
        a_pos = 0
        if read in readkeys_sama:
            a_rname = sama_reads[read][0]
            a_count = sama_reads[read][1]
            a_exact = sama_reads[read][2]
            a_edits = sama_reads[read][3]
        else:            
            a_rname = ""
            a_count = 0
            a_exact = 0
            a_edits = 0

        b_pos = 1
        if read in readkeys_samb:
            b_rname = samb_reads[read][0]
            b_count = samb_reads[read][1]
            b_exact = samb_reads[read][2]
            b_edits = samb_reads[read][3]
        else:
            b_rname = ""
            b_count = 0
            b_exact = 0
            b_edits = 0

        # TODO: check for reads that are aligned in both sam files, but to different fusions.
        # TODO: Output length-adjusted counts?

        # make sure read id already included, add if not
        # use rname from A if available, else use B - if A is "*" and B is not set could end up with '' - HdR
        read_id = b_rname if (a_rname == '' or a_rname == '*' ) else a_rname
        if read_id not in counts:
            addreadid += 1
            if addedids:
                addedids += ", {0}".format(read_id)
            else:
                addedids = read_id
            counts[read_id] = dict(count_template)

        # check for zero reads
        if a_count == 0 and b_count == 0:
            both_unaligned += 1
        elif a_count == 0:
            if b_count == 1:
                if b_exact == 1:
                    b_single_exact += 1
                    counts[read_id]['b_single_exact'] += 1
                elif b_exact == 0:
                    b_single_inexact += 1
                    counts[read_id]['b_single_inexact'] += 1
            elif b_count > 1:
                if b_exact == 1:
                    b_multi_exact += 1
                elif b_exact == 0:
                    b_multi_inexact += 1
            else:
                counts = {}
                logger.error("Invalid b_count, {0}, with a_count = 0 for read id '{1}'".format(b_count, read_id))
                break
        elif b_count == 0:
            if a_count == 1:
                if a_exact == 1:
                    a_single_exact += 1
                    counts[read_id]['a_single_exact'] += 1
                elif a_exact == 0:
                    a_single_inexact += 1
                    counts[read_id]['a_single_inexact'] += 1
            elif a_count > 1:
                if a_exact == 1:
                    a_multi_exact += 1
                elif a_exact == 0:
                    a_multi_inexact += 1
            else:
                counts = {}
                logger.error("Invalid a_count, {0}, with b_count = 0 for read id '{1}'".format(a_count, read_id))
                break
        else: #both a_count and b_count are not zero - either 1 or > 1
            if a_count == 1:
                if a_exact == 1:
                    if b_count == 1:
                        if b_exact == 1:
                            # Impossible as we set a_pos = 0 and b_pos = 1
                            # a=1=1/b=1=1
                            if a_pos != '' and b_pos != '' and a_pos == b_pos and a_rname == b_rname:
                                both_single_exact_same += 1
                                counts[read_id]['both_single_exact_same'] += 1
                            else:
                                both_single_exact_diff += 1
                                counts[read_id]['both_single_exact_diff'] += 1
                        else: # b_exact == 0
                            a_se_b_si += 1
                            counts[read_id]['a_exact_b_inexact'] += 1
                    else: # b_count > 1
                        if b_exact == 0:
                            a_se_b_mi += 1
                        else:
                            a_se_b_me += 1
                else: # a_exact == 0
                    if b_count == 1:
                        if b_exact == 0:
                            if a_pos != '' and b_pos != '' and a_pos == b_pos and a_rname == b_rname:
                                both_single_inexact_same += 1
                                counts[read_id]['both_single_inexact_same'] += 1
                            else:
                                both_single_inexact_diff += 1
                                counts[read_id]['both_single_inexact_diff'] += 1
                                if a_edits == b_edits:
                                    both_inexact_diff_equal += 1
                                    counts[read_id]['both_inexact_diff_equal'] += 1
                                elif a_edits > b_edits:
                                    both_inexact_diff_b_better += 1
                                    counts[read_id]['both_inexact_diff_b_better'] += 1
                                else: # a_edits < b_edits:
                                    both_inexact_diff_a_better += 1
                                    counts[read_id]['both_inexact_diff_a_better'] += 1
                        else: #b_exact == 1
                            a_si_b_se += 1
                            counts[read_id]['b_exact_a_inexact'] += 1
                    else: #b_count > 1
                        if b_exact == 0:
                            a_si_b_mi += 1
                        else: #b_exact == 1
                            a_si_b_me += 1
            else: #a_count > 1
                if a_exact == 0:
                    if b_count == 1:
                        if b_exact == 0:
                            a_mi_b_si += 1
                        else: #b_exact == 1
                            a_mi_b_se += 1
                    else: #b_count > 1
                        if b_exact == 0:
                            both_multi_inexact += 1
                        else: #b_exact == 1
                            a_mi_b_me += 1
                else: #a_exact == 1
                    if b_count == 1:
                        if b_exact == 0:
                            a_me_b_si += 1
                        else: #b_exact == 1
                            a_me_b_se += 1
                    else: #b_count > 1:
                        if b_exact == 0:
                            a_me_b_mi += 1
                        else: #b_exact == 1
                            both_multi_exact += 1

        # check this, not sure why chg to output to console - HdR!!!
        #if opts.output or opts.totals:
        #    if a_count > 1:
        #        a_rname = '*'
        #        a_pos = '*'
        #    if b_count > 1:
        #        b_rname = '*'
        #        b_pos = '*'

    # total things up if no errors
    if counts:
        logger.debug("Time to process all reads: {0:n} seconds, {1:n} records".format(time.time() - tstart, cnt))

        # calculate total count and do sanity check
        total_count = a_single_exact + a_single_inexact + a_multi_exact + a_multi_inexact + b_single_exact + b_single_inexact + b_multi_exact + b_multi_inexact + both_unaligned + both_single_exact_same + both_single_exact_diff + both_single_inexact_same + both_single_inexact_diff + both_multi_exact + both_multi_inexact + a_se_b_si + a_si_b_se + a_se_b_me + a_me_b_se + a_se_b_mi + a_mi_b_se + a_si_b_me + a_me_b_si + a_si_b_mi + a_mi_b_si + a_me_b_mi + a_mi_b_me
        if len(readkeys) != total_count:
            logger.error("Total count, {0:n}, did not match read count, {1:n}".format(total_count, len(readkeys)))
        if addreadid:
            logger.info("{0:n} records were added that did not have a matching fusion read key".format(addreadid))

        # create totals list
        totals = [a_single_exact, a_single_inexact, a_multi_exact, a_multi_inexact,
                  b_single_exact, b_single_inexact, b_multi_exact, b_multi_inexact,
                  both_single_exact_same, both_single_exact_diff,
                  both_single_inexact_same, both_single_inexact_diff,
                  both_inexact_diff_equal, both_inexact_diff_a_better,
                  both_inexact_diff_b_better, both_multi_exact, both_multi_inexact,
                  a_se_b_si, a_si_b_se, a_se_b_me, a_me_b_se, a_se_b_mi, a_mi_b_se,
                  a_si_b_me, a_me_b_si, a_si_b_mi, a_mi_b_si, a_me_b_mi, a_mi_b_me,
                  total_count]

    return (totals, counts)

# write read counts to file or console
def write_counts(filename, counts):
    #header = ['fusion', 'both_single_inexact_same', 'both_single_exact_diff', 'both_single_exact_same', 
    #          'both_single_inexact_diff', 'both_inexact_diff_equal', 'a_single_exact', 'b_single_exact', 
    #          'a_exact_b_inexact', 'b_exact_a_inexact','a_single_inexact', 'b_single_inexact', 
    #          'both_inexact_diff_a_better', 'both_inexact_diff_b_better']
    header = ['FUSION_ID','BOTH_EXACT','BOTH_INEXACT_EQUAL','SAM_A_ONLY_EXACT','SAM_B_ONLY_EXACT',
              'SAM_A_EXACT_SAM_B_INEXACT','SAM_B_EXACT_SAM_A_INEXACT',
              'SAM_A_ONLY_SINGLE_INEXACT','SAM_B_ONLY_SINGLE_INEXACT',
              'SAM_A_INEXACT_BETTER','SAM_B_INEXACT_BETTER']
    fields = ['fusion', 'both_single_exact_diff', 'both_inexact_diff_equal', 
              'a_single_exact', 'b_single_exact',
              'a_exact_b_inexact', 'b_exact_a_inexact',
              'a_single_inexact', 'b_single_inexact', 
              'both_inexact_diff_a_better', 'both_inexact_diff_b_better']

    # Note: these fields are not currently written out in pearl script:
    #       'both_single_exact_same', 'both_single_inexact_same', 'both_single_inexact_diff']
    # You may use the code below instead of doing output manually to match old file header and format exactly
    # csvwriter = csv.DictWriter(f, header, extrasaction='raise')
    # csvwriter.writeheader()
    # csvwriter.writerows(counts.values())

    try:
        f = open(filename, 'w') if filename.lower() != "stdout" else sys.stdout
        f.write(",".join(header))
        f.write("\n")
        for row in counts.itervalues():
            rowstr = ""
            for field in fields:
                if rowstr:
                    rowstr += ","
                rowstr += str(row[field])
            f.write("{0}\n".format(rowstr))
        if f != sys.stdout:
            f.close()
    except Exception as e:
        logger.error("Unable to write counts output to file '{0}':\n{1}".format(filename, e))
        if f and f != sys.stdout:
            f.close()

# write totals to file or console
def write_totals(filename, totals):
    header = ['a_single_exact', 'a_single_inexact', 'a_multi_exact', 'a_multi_inexact',
              'b_single_exact', 'b_single_inexact', 'b_multi_exact', 'b_multi_inexact', 'both_single_exact_same',
              'both_single_exact_diff', 'both_single_inexact_same', 'both_single_inexact_diff',
              'both_inexact_diff_equal', 'both_inexact_diff_a_better', 'both_inexact_diff_b_better', 'both_multi_exact',
              'both_multi_inexact', 'a_single_exact_b_single_inexact', 'a_single_inexact_b_single_exact', 'a_single_exact_b_multi_exact',
              'a_multi_exact_b_single_exact', 'a_single_exact_b_multi_inexact', 'a_multi_inexact_b_single_exact', 'a_single_inexact_b_multi_exact',
              'a_multi_exact_b_single_inexact', 'a_single_inexact_b_multi_inexact', 'a_multi_inexact_b_single_inexact', 'a_multi_exact_b_multi_inexact',
              'a_multi_inexact_b_multi_exact', 'total_count']

    try:
        f = open(filename, 'w') if filename.lower() != "stdout" else sys.stdout
        f.write("Count totals:\n")
        for idx in range(0, len(header)):
            f.write("{0}:\t{1}\t{2}\n".format(idx+1, header[idx], totals[idx]))
        if f != sys.stdout:
            f.close()
    except Exception as e:
        logger.error("Unable to write totals to file '{0}':\n{1}".format(filename, e))
        if f and f != sys.stdout:
           f.close()


def main():
    retcode = -1
    random.seed()
    tinitial = time.time()

    # get base temporary file name
    basefilename = get_temp_basefilename()

    # parse command line options and setup logger to handle all output messages
    opts = parse_options()
    filename = opts.log if opts.log else "samcompare_{0}.log".format(basefilename)
    level = logging.DEBUG if opts.debug else logging.INFO
    setLogger(filename, level)

    logger.info("SAM Compare version {0}".format(__version__))
    logger.debug("Opts: {0}".format(opts))

    # real files - temp!!!
    #opts.fastq = "/bio/ufhpc/om/projects/mcintyre/sam-compare/ase/s_1_combined.fq"
    #opts.sama = "/bio/ufhpc/om/projects/mcintyre/sam-compare/ase/rna_1_to_berlin-fusions39.sam"
    #opts.samb = "/bio/ufhpc/om/projects/mcintyre/sam-compare/ase/rna_1_to_c1674-fusions34.sam"
    #opts.fusion = "/bio/ufhpc/om/projects/mcintyre/sam-compare/ase/fb526-si-fusions.tsv"
    # use ./sam_compare.py -l 54 -f fn -q fn -A fn -B fn for testing
    #
#    opts.fastq = "/scratch/lfs/hdr/dev/Tms_3.fastq"
#    opts.sama = "/scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/aln_to_orthos_final_uniq_pipe/Tms_3_to_Tdu_ortho.sam"
#    opts.samb = "/scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/aln_to_orthos_final_uniq_pipe/Tms_3_to_Tpr_ortho.sam"
#    opts.fusion = "/scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/references/Tdu_common_ortho_clean.tsv"
#    opts.counts = "countstms3.csv"
#    opts.totals = "totalstms3.txt"
    # cat Tm?_*.fastq > /scratch/lfs/hdr/dev/Tm_?.fastq
    # use ./sam_compare.py -l 100 -f fn -q fn -A fn -B fn -n -g stdout -d  (for testing)

    """
    perl sam-compare-custom_readID_w_space.pl
    -l 100 -f /scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/references/Tdu_common_ortho_clean.tsv
    -r /scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/references/Tdu_common_ortho_clean.fasta
    -q /scratch/lfs/hdr/dev/Tms_3.fastq
    -A /scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/aln_to_orthos_final_uniq_pipe/Tms_3_to_Tdu_ortho.sam
    -B /scratch/lfs/malex/projects/mcintyre/sam-compare/test/2014-03-18/perl/trago_example/aln_to_orthos_final_uniq_pipe/Tms_3_to_Tpr_ortho.sam
    >Tms_3_2_Tdu_Tpr.csv 2> Tms_3_2_Tdu_Tpr.log
    """

    # FASTQ - the option to skip saves a significant amount of time, very useful for development/testing
    tstart = time.time()
    fqids = {}
    if opts.nofqids:
        logger.info("Skipping extracting IDs from fastq file, -n option specified")
    else:
        logger.info("Extracting IDs from fastq file...")
        fqids = get_fastq_ids(opts.fastq, basefilename)
        logger.debug("Time to get FASTQ file ids: {0:n} seconds".format(time.time() - tstart))

    # make sure we got a request to skip checking ids or we got them
    if opts.nofqids or fqids:
        gc.collect()
        logger.info("Processing SAMA file...")
        tstart = time.time()
        sama_reads = process_sam_file("A", opts.sama, fqids, opts.length)
        logger.debug("Time to process SAMA file: {0:n} seconds".format(time.time() - tstart))
        if sama_reads:
            logger.info("Processing SAMB file...")
            tstart = time.time()
            samb_reads = process_sam_file("B", opts.samb, fqids, opts.length)
            logger.debug("Time to process SAMB file: {0:n} seconds".format(time.time() - tstart))
            if samb_reads:
                fqids = {}
                gc.collect()
                # additional processing
                logger.info("Loading fusion data...")
                tstart = time.time()
                fusions = get_fusions(opts.fusion)
                logger.debug("Time to process fusion data: {0:n} seconds, {1:n} records".format(time.time() - tstart, len(fusions)))
                if fusions:
                    logger.info("Calculating read counts...")
                    tstart = time.time()
                    totals, counts = process_read_counts(fusions, sama_reads, samb_reads)
                    logger.debug("Time to calculate read counts: {0:n} seconds".format(time.time() - tstart))

                    if totals and counts:
                        # write counts file
                        logger.info("Writing counts file...")
                        tstart = time.time()
                        filename = opts.counts if opts.counts else "counts_{0}.csv".format(basefilename)
                        write_counts(filename, counts)
                        logger.debug("Time to write counts file: {0:n} seconds".format(time.time() - tstart))

                        # write totals file
                        logger.info("Writing totals file...")
                        tstart = time.time()
                        filename = opts.totals if opts.totals else "totals_{0}.csv".format(basefilename)
                        write_totals(filename, totals)
                        logger.debug("Time to write totals file: {0:n} seconds".format(time.time() - tstart))

                        # all is good
                        retcode = 0
                    else:
                        logger.error("Unable to get totals and counts")
                else:
                    logger.error("No fusion data available for processing")

    # it took this long...
    logger.info("All done. Total running time: {0:n} seconds".format(time.time() - tinitial))
    sys.exit(retcode)

if __name__=='__main__':
    main()
