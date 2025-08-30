#!/usr/bin/env python3

import argparse


parser = argparse.ArgumentParser(description="Dereplicates and trims sam file reads with counts")

parser.add_argument("--in_file", type=argparse.FileType("r"), help="sam file to trim")
parser.add_argument("--out_file", type=argparse.FileType("w"), help="output sam name")
parser.add_argument(
    "--min_count", type=int, default=1, help="min count for a sequence to be included.  default 1"
)
parser.add_argument(
    "--ltrim", type=int, default=0, help="Number of nucleotides to trim from the 5' (left) end of the read.  default 0"
)
parser.add_argument(
    "--rtrim", type=int, default=0, help="Number of nucleotides to trim from the 3' (right) end of the read.  default 0"
)
parser.add_argument(
    "--min_len", type=int, default=20, help="Minimum read length (after trimming) to keep.  default 20"
)

args = parser.parse_args()

seq_dict = {}
total_reads = 0
for line in args.in_file:
    if not line.startswith("@"):
        line = line.split("\t")
        if not line[1] in ("0", "16"):
            continue

        if args.ltrim + args.rtrim + args.min_len > len(line[9]):
            continue
        if args.ltrim + args.rtrim > 0:
            cig = []
            ntcount = ''
            for c in line[5]:
                if c.isnumeric():
                    ntcount += c
                else:
                    cig.append([c,int(ntcount)])
                    ntcount = ''

            if args.ltrim > 0:
                line[9] = line[9][args.ltrim:]
                to_trim = args.ltrim
                line[3] = int(line[3]) + to_trim
                for i in range(0, len(cig)):
                    if to_trim == 0:
                        break
                    if cig[i][0] == "D":
                        line[3] += cig[i][1]
                        cig[i] = ''
                    elif cig[i][1] <= to_trim:
                        to_trim -= cig[i][1]
                        if cig[i][0] == "I":
                            line[3] -= cig[i][1]
                        cig[i] = ''
                    else:
                        if cig[i][0] == "I":
                            line[3] -= to_trim
                        cig[i][1] -= to_trim
                        break
            line[3] = str(line[3])
            if args.rtrim > 0:
                to_trim = args.rtrim
                line[9] = line[9][:-args.rtrim]
                for i in range(len(cig)-1, -1, -1):
                    if cig[i][0] == "D":
                        cig[i] = ''
                    elif cig[i][1] <= to_trim:
                        to_trim -= cig[i][1]
                        cig[i] = ''
                    else:
                        cig[i][1] -= to_trim
                        break

            line[5] = ''
            for entry in cig:
                if entry:
                    line[5] += str(entry[1]) + entry[0]


        total_reads += 1
        try:
            seq_dict[line[9]]['count'] += 1
        except KeyError:
            seq_dict[line[9]] = {'ref' : line[2], 'start' : line[3], 'score' : line[4], 'cig' : line[5], 'count' : 1}

args.in_file.close()

print(
    f"{total_reads} sequences collected.  Dereplicated to {len(seq_dict)} unique sequences.  Writing output file"
)
sorted_seqs = sorted(seq_dict.keys(), key=lambda x: int(seq_dict[x]['start']))

counter = 1
for seq in sorted_seqs:
    if seq_dict[seq]['count'] >= args.min_count:
        args.out_file.write(f"{counter}-{seq_dict[seq]['count']}\t0\t{seq_dict[seq]['ref']}\t{seq_dict[seq]['start']}\t{seq_dict[seq]['score']}\t{seq_dict[seq]['cig']}\t*\t0\t0\t{seq}\t*\n")
        counter += 1
    else:
        break
args.out_file.close()
