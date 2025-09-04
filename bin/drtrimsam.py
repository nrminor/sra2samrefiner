#!/usr/bin/env python3

import os
import argparse
import pysam
pysam.set_verbosity(0)

parser = argparse.ArgumentParser(description="Dereplicates and trims cram file reads with counts")

parser.add_argument("--in_file", type=str, help="merged cram file to trim both ends")
parser.add_argument("--in_file2", type=str, help="unmerged cram file to trim 5' only")
parser.add_argument("--out_file", type=str, help="output sam name")
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

if not (args.in_file or args.in_file2):
    print('No input files specified')
    exit(1)
if not args.out_file:
    print('No output files specified')
    exit(1)

cig_dict = {
    0 : "M",
    1 : "I",
    2 : "D",
    4 : "S",
    5 : "H",
    }

seq_dict = {}
total_reads = 0
header_dict = None
below_min = 0

def parse_file(file, left, right):
    global header_dict
    global total_reads
    global seq_dict
    global below_min

    if file.endswith(".sam"):
        try:
            in_fh = pysam.AlignmentFile(file, "rs")
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
        print("{file} doesn't have known extention (sam/bam/cram)")
        return
        
    try:
        if not header_dict:
            header_dict = in_fh.header.to_dict()
        for read in in_fh:
            if not read.flag in (0, 16):
                continue
            total_reads += 1
            length = read.query_length
            if length == 0:
                length = read.infer_query_length()
            if left + right + args.min_len > length:
                below_min += 1
                continue
            if left + right > 0:
                cig = read.cigartuples
                if left > 0:
                    read.query_sequence = read.query_sequence[left:]
                    to_trim = left
                    start_pos = read.reference_start + to_trim
                    for i in range(0, len(cig)):
                        if to_trim == 0:
                            break
                        if cig[i][0] == 5:
                            cig[i] = ''
                        elif cig[i][0] == 2:
                            start_pos += cig[i][1]
                            cig[i] = ''
                        elif cig[i][1] <= to_trim:
                            to_trim -= cig[i][1]
                            if cig[i][0] in (1, 4):
                                start_pos -= cig[i][1]
                            cig[i] = ''
                        else:
                            if cig[i][0] in (1, 4):
                                start_pos-= to_trim
                            cig[i] = (cig[i][0], (cig[i][1] - to_trim))
                            break

                if right > 0:
                    to_trim = right
                    read.query_sequence = read.query_sequence[:-right]
                    for i in range(len(cig)-1, -1, -1):
                        if cig[i][0] in (2, 5):
                            cig[i] = ''
                        elif cig[i][1] <= to_trim:
                            to_trim -= cig[i][1]
                            cig[i] = ''
                        else:
                            cig[i] = (cig[i][0], (cig[i][1] - to_trim))
                            break

                cigstr = ''
                n_s_c = 0
                for entry in cig:
                    if entry:
                        n_s_c += 1
                        cigstr += str(entry[1]) + cig_dict[entry[0]]

                try:
                    seq_dict[read.query_sequence]['count'] += 1
                except KeyError:
                    seq_dict[read.query_sequence] = {'ref' : read.reference_name, 'start' : start_pos, 'score' : read.mapping_quality, 'cig' : cigstr, 'count' : 1}
    except:
        print(f"Failed parsing of {file}")
        exit(1)
    in_fh.close()

if args.in_file:
    if not os.path.isfile(args.in_file):
        print(f"Can't locate {args.in_file}")
    else:
        parse_file(args.in_file, args.trim, args.trim)

if args.in_file2:
    if not os.path.isfile(args.in_file2):
        print(f"Can't locate {args.in_file2}")
    else:
        parse_file(args.in_file2, args.trim, 0)

print(
    f"{total_reads} sequences collected.  Dereplicated to {len(seq_dict)} unique sequences, {below_min} discarded for length.  Writing output file"
)
sorted_seqs = sorted(seq_dict.keys(), key=lambda x: int(seq_dict[x]['start']))

wt = ''
if args.out_file.endswith(".sam"):
    wt = "ws"
if args.out_file.endswith(".bam"):
    wt = "wb"
if args.out_file.endswith(".cram"):
    wt = "wc"
with pysam.AlignmentFile(args.out_file, wt, header=header_dict) as out_fh:
    counter = 1
    for seq in sorted_seqs:
        if seq_dict[seq]['count'] >= args.min_count:
            read = pysam.AlignedSegment(out_fh.header)
            read.query_name = f"{counter}-{seq_dict[seq]['count']}"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = seq_dict[seq]['start']
            read.mapping_quality = seq_dict[seq]['score']
            read.cigarstring = seq_dict[seq]['cig']
            read.template_length = 0
            read.query_sequence = seq
            read.query_qualities = "*"
            out_fh.write(read)
            counter += 1
