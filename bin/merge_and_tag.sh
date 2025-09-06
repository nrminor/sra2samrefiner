#!/usr/bin/env bash
# merge_and_tag.sh — Robust BBTools merging + header tagging + concat
# Produces: <ACCESSION>.merged.fastq.gz and <ACCESSION>_ihist_merge.txt
#
# Usage:
#   merge_and_tag.sh -1 R1.fastq.gz -2 R2.fastq.gz -o ACCESSION [-t THREADS]
#
# Example:
#   merge_and_tag.sh -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o SRR123456 -t 8
#
# Notes:
#   - Requires: bbmerge.sh, bbrename.sh, reformat.sh, gzip
#   - Safe if any stream is empty (creates valid empty .gz placeholders)
#   - Writes final outputs to the current working directory

set -euo pipefail

# --------------------------- arg parsing & help ---------------------------- #
usage() {
	cat <<EOF
Usage: $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -o ACCESSION [-t THREADS]

Required:
  -1  Path to R1 FASTQ(.gz)
  -2  Path to R2 FASTQ(.gz)
  -o  Output accession/prefix (final file: <ACCESSION>.merged.fastq.gz)

Optional:
  -t  Threads for BBTools (default: 4)

Outputs (in CWD):
  <ACCESSION>.merged.fastq.gz
  <ACCESSION>_ihist_merge.txt
EOF
}

R1=""
R2=""
OUT=""
THREADS=4

while getopts ":1:2:o:t:h" opt; do
	case "$opt" in
	1) R1="$OPTARG" ;;
	2) R2="$OPTARG" ;;
	o) OUT="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	h)
		usage
		exit 0
		;;
	\?)
		echo "Unknown option: -$OPTARG" >&2
		usage
		exit 2
		;;
	:)
		echo "Missing arg for -$OPTARG" >&2
		usage
		exit 2
		;;
	esac
done

[[ -z "${R1}" || -z "${R2}" || -z "${OUT}" ]] && {
	usage
	exit 2
}
[[ -r "$R1" ]] || {
	echo "Cannot read R1: $R1" >&2
	exit 1
}
[[ -r "$R2" ]] || {
	echo "Cannot read R2: $R2" >&2
	exit 1
}

# ------------------------------- checks ----------------------------------- #
need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing dependency: $1" >&2
	exit 127
}; }
need bbmerge.sh
need bbrename.sh
need reformat.sh
need gzip

# ------------------------------- paths ------------------------------------ #
FINAL="${OUT}.merged.fastq.gz"
IHIST="${OUT}_ihist_merge.txt"

tmpdir="$(mktemp -d -t merge_and_tag.XXXXXX)"
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

MERGED_READY="$tmpdir/merged.tagged.ready.fastq.gz"
U1_TMP="$tmpdir/unmerged_R1.tagged.tmp.fastq.gz"
U2_TMP="$tmpdir/unmerged_R2.tagged.tmp.fastq.gz"
U1_READY="$tmpdir/unmerged_R1.tagged.ready.fastq.gz"
U2_READY="$tmpdir/unmerged_R2.tagged.ready.fastq.gz"

# ------------------------------ pipeline ---------------------------------- #
# 1) Merge pairs; label streams in-flight.
bbmerge.sh \
	in1="$R1" \
	in2="$R2" \
	out=>(bbrename.sh in=stdin out="$MERGED_READY" prefix=MERGED_) \
	outu1=>(bbrename.sh in=stdin out="$U1_TMP" prefix=UNMERGED_) \
	outu2=>(bbrename.sh in=stdin out="$U2_TMP" prefix=UNMERGED_) \
	qtrim=t \
	ihist="$IHIST" \
	threads="$THREADS" \
	-eoom

# 2) Ensure the three outputs exist even if empty (valid empty gzip files).
for f in "$MERGED_READY" "$U1_TMP" "$U2_TMP"; do
	[[ -e "$f" ]] || gzip -c </dev/null >"$f"
done

# 3) Normalize unmerged names to /1 and /2 (handles empty inputs fine).
reformat.sh \
	in1="$U1_TMP" \
	in2="$U2_TMP" \
	out1="$U1_READY" \
	out2="$U2_READY" \
	addslash=t spaceslash=f \
	threads="$THREADS"

# 4) Build final file in a deterministic order; single redirection for atomicity.
cat \
	"$MERGED_READY" \
	"$U1_READY" \
	"$U2_READY" \
	>"$FINAL"

# Done.
echo "Wrote: $FINAL"
echo "Histogram: $IHIST"
