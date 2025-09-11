#!/usr/bin/env python3

import os
import argparse
import pysam
pysam.set_verbosity(0) # make pysam not whine about not having indexes that aren't needed

'''
Get command line arguments for input/output file paths, the minimal read count and length for a sequence
to be output, and the number of nts to trim
'''
parser = argparse.ArgumentParser(description="Dereplicates and trims sam/bam/cram file reads with counts")

parser.add_argument("--in_file", type=str, help="merged cram file to trim both ends")
parser.add_argument("--in_file2", type=str, help="unmerged cram file to trim 5' only")
parser.add_argument("--out_file", required=True, type=str, help="output sam name")
parser.add_argument(
    "--min_count", type=int, default=1, help="min count for a sequence to be included.  default 1"
)
parser.add_argument(
    "--trim", type=int, default=0, help="Number of nucleotides to trim from ends of the reads.  default 0.  Negative numbers are treated as 0"
)
parser.add_argument(
    "--min_len", type=int, default=20, help="Minimum read length (after trimming) to keep.  default 20"
)

args = parser.parse_args()

# check if at least one input file is provided
if not (args.in_file or args.in_file2):
    print('No input files specified')
    exit(1)

# make sure output type is specified
wt = ''
if args.out_file.endswith(".sam"):
    wt = "w"
elif args.out_file.endswith(".bam"):
    wt = "wb"
elif args.out_file.endswith(".cram"):
    wt = "wc"
else:
    print('unknown output type.  Please specifiy sam/bam/cram')
    exit()

'''
# dictionary for pysam cigar int matched to cigar letter code
# not used in code, included for reference
cig_dict = {
    0 : "M",
    1 : "I",
    2 : "D",
    4 : "S",
    5 : "H",
    }
'''

seq_dict = {} # dictionary for each sequence passing the length cutoff
header_dict = None # stores header of first parsed file to keep for the output
total_reads = 0 # keeps the count of the original number of reads in the input
below_min = 0 # keeps the count of reads that don't pass the length cutoff

def parse_file(file, left, right):
    '''
    Called to parse input mapped filed to dereplicate and trim reads
    Parameters:
    file - str, name of the file to be parsed
    left - int, number of nts to trim from 5' of read
    right - int, number of nts to trim from 3' of read
    Functionality:
    Uses global vairables to store info from parsing passed file.  The file is first opened via pysam.
    If an error is thrown by the initial pysam reading, or an invalid extention is used, further
    parsing is skipped.  Each line is parsed by trimming and recalculating the cigar entry and
    reference start position.

    Returns nothing
    '''
    global seq_dict
    global header_dict
    global total_reads
    global below_min
    # check if file can be opened by pysam
    if file.endswith(".sam"):
        try:
            in_fh = pysam.AlignmentFile(file, "r")
        except:
            print("Unable to open {file} as a sam file")
            return
    elif file.endswith(".bam"):
        try:
            in_fh = pysam.AlignmentFile(file, "rb")
        except:
            print("Unable to open {file} as a bam file")
            return
    elif file.endswith(".cram"):
        try:
            in_fh = pysam.AlignmentFile(file, "rc")
        except:
            print("Unable to open {file} as a cram file")
            return
    else:
        print(f"{file} doesn't have known extention (sam/bam/cram), skipping")
        return

    try:
        if not header_dict: # get header if one doesn't exits for output
            header_dict = in_fh.header.to_dict()
        for read in in_fh:
            if not read.flag in (0, 16): # only consider primary mapping of read
                continue
            total_reads += 1
            length = read.query_length
            if length == 0: # Should never be the case with primay mapping, but included as fail safe
                length = read.infer_query_length()
            if left + right + args.min_len > length: # don't keep sequences that are too short
                below_min += 1
                continue
            seq = read.query_sequence
            
            ltrim = left
            rtrim = right
            if read.flag == 16:
                ltrim = right
                rtrim = left
            if ltrim > 0:
                seq = seq[ltrim:] # trim actually nt seq
            if rtrim > 0:
                seq = seq[:-rtrim] # trim actually nt seq
            try: # add to count of unique seq if exist and continue to skip unnecssary reprocessing
                seq_dict[seq]['count'] += 1
                continue
            except KeyError:
                cig = read.cigartuples
                if ltrim > 0: # 5' recalc
                    start_pos = read.reference_start + ltrim # initial recalc of start pos
                    # iterate through cigar tuple list to recalculate it for trimmed sequence
                    to_trim = ltrim
                    for i in range(0, len(cig)):
                        if to_trim == 0: # don't continue loop if all the trimming has been accounted for
                            break
                        if cig[i][0] == 5: # hard clipping entries removed
                            cig[i] = '' # functionally remove entry
                        elif cig[i][0] == 2: # deletion required start pos recalc and entry removal
                            start_pos += cig[i][1]
                            cig[i] = '' # functionally remove entry
                        elif cig[i][1] <= to_trim: # if the nt seq described by the cigar entry has been clipped
                            to_trim -= cig[i][1] # recalce remaining trim length to account for
                            if cig[i][0] in (1, 4): # insertions and soft clipping require start pos recalc
                                start_pos -= cig[i][1]
                            cig[i] = '' # functionally remove entry
                        else: # if nt clip doesn't cover all of the cigar entry
                            if cig[i][0] in (1, 4): # insertions and soft clipping require start pos recalc
                                start_pos-= to_trim
                            cig[i] = (cig[i][0], (cig[i][1] - to_trim)) # recalc cigar run length
                            break

                if rtrim > 0: # 3' recalc
                    # iterate through cigar tuple list in reverse to recalculate it for trimmed sequence
                    to_trim = rtrim
                    for i in range(len(cig)-1, -1, -1):
                        if cig[i][0] in (2, 5): # hard clipping and del entries removed
                            cig[i] = '' # functionally remove entry
                        elif cig[i][1] <= to_trim: # if the nt seq described by the cigar entry has been clipped
                            to_trim -= cig[i][1]  # recalce remaining trim length to account for
                            cig[i] = '' # functionally remove entry
                        else: # if nt clip doesn't cover all of the cigar entry
                            cig[i] = (cig[i][0], (cig[i][1] - to_trim)) # recalc cigar run length
                            break

                re_cig = []
                for entry in cig: # remove empty entries so pysam can use
                    if entry:
                        re_cig.append(entry)

                # storing values needed for output
                seq_dict[seq] = {'ref' : read.reference_id, 'start' : start_pos, 'score' : read.mapping_quality, 'cig' : re_cig, 'count' : 1}

    except Exception as err: # catch any error in line parsing and report with exit to prevent innacurate output
        print(f"Failed parsing of {file}")
        print(err)
        exit(1)
    in_fh.close()

if args.in_file: # if passed file for both end trimming
    if not os.path.isfile(args.in_file): # check file path
        print(f"Can't locate {args.in_file}")
    else:
        parse_file(args.in_file, args.trim, args.trim)

if args.in_file2: # if passed file for only 5' end trimming
    if not os.path.isfile(args.in_file2): # check file path
        print(f"Can't locate {args.in_file2}")
    else:
        parse_file(args.in_file2, args.trim, 0)

l_c = len(seq_dict)

print(f"{total_reads} sequences collected.  Dereplicated to {l_c} unique sequences, {below_min} discarded for length.  Writing output file")

with open(f'{args.out_file}.count', 'w') as count_out: # output file with read count in map output file to be used in greater workflow
    count_out.write(f"{l_c}")
# sort sequences by reference start pos
sorted_seqs = sorted(seq_dict.keys(), key=lambda x: int(seq_dict[x]['start']))

# write output file
with pysam.AlignmentFile(args.out_file, wt, header=header_dict) as out_fh:
    counter = 1
    for seq in sorted_seqs:
        if seq_dict[seq]['count'] >= args.min_count: # check unique count passes min
            read = pysam.AlignedSegment(out_fh.header) # create pysam read and assign required values
            read.query_name = f"{counter}-{seq_dict[seq]['count']}"
            read.flag = 0
            read.reference_id = seq_dict[seq]['ref']
            read.reference_start = seq_dict[seq]['start']
            read.mapping_quality = seq_dict[seq]['score']
            read.cigartuples = seq_dict[seq]['cig']
            read.template_length = 0
            read.query_sequence = seq
            read.query_qualities = "*"

            out_fh.write(read)
            counter += 1
