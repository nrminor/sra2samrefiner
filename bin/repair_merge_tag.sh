#!/usr/bin/env bash
# repair_merge_tag.sh — Repair, merge, and tag paired-end reads with BBTools
# Uses repair.sh to separate mated from orphan reads before merging
# Auto-detects R1/R2 files from accession prefix
#
# Produces: <ACCESSION>.merged.fastq.gz and <ACCESSION>_ihist_merge.txt
#
# Usage:
#   repair_merge_tag.sh -o ACCESSION [-d DIRECTORY] [-t THREADS]
#
# Example:
#   repair_merge_tag.sh -o SRR123456 -d /path/to/fastqs -t 8
#   repair_merge_tag.sh -o SRR123456  # looks in current directory
#
# Notes:
#   - Requires: repair.sh, bbmerge.sh, bbrename.sh, reformat.sh, gzip
#   - Auto-detects R1/R2 files matching: <ACCESSION>*[_.]R?1*.f*q* patterns
#   - Runs repair.sh first to separate mated from orphan reads
#   - Safe if any stream is empty (creates valid empty .gz placeholders)
#   - Writes final outputs to the current working directory

set -euo pipefail

# --------------------------- arg parsing & help ---------------------------- #
usage() {
	cat <<EOF
Usage: $(basename "$0") -o ACCESSION [-d DIRECTORY] [-t THREADS]

Required:
  -o  Output accession/prefix (will auto-find matching R1/R2 files)

Optional:
  -d  Directory to search for R1/R2 files (default: current directory)
  -t  Threads for BBTools (default: 4)

Outputs (in CWD):
  <ACCESSION>.merged.fastq.gz
  <ACCESSION>_ihist_merge.txt

The script will automatically find R1/R2 files matching patterns like:
  <ACCESSION>_R1.fastq.gz, <ACCESSION>_R2.fastq.gz
  <ACCESSION>_1.fastq.gz, <ACCESSION>_2.fastq.gz
  <ACCESSION>.R1.fq.gz, <ACCESSION>.R2.fq.gz
  etc.
EOF
}

OUT=""
SEARCH_DIR="."
THREADS=4

while getopts ":o:d:t:h" opt; do
	case "$opt" in
	o) OUT="$OPTARG" ;;
	d) SEARCH_DIR="$OPTARG" ;;
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

[[ -z "${OUT}" ]] && {
	echo "Error: accession/prefix (-o) is required" >&2
	usage
	exit 2
}

[[ -d "$SEARCH_DIR" ]] || {
	echo "Directory does not exist: $SEARCH_DIR" >&2
	exit 1
}

# ---------------------------- auto-detect R1/R2 --------------------------- #
# Find R1 and R2 files matching the accession prefix
# Priority order of patterns to check
find_reads() {
	local prefix="$1"
	local dir="$2"
	local read_num="$3" # 1 or 2

	# Array of patterns to try, in order of preference
	local patterns=(
		"${prefix}_R${read_num}.fastq.gz"
		"${prefix}_R${read_num}.fq.gz"
		"${prefix}_R${read_num}.fastq"
		"${prefix}_R${read_num}.fq"
		"${prefix}_${read_num}.fastq.gz"
		"${prefix}_${read_num}.fq.gz"
		"${prefix}_${read_num}.fastq"
		"${prefix}_${read_num}.fq"
		"${prefix}.R${read_num}.fastq.gz"
		"${prefix}.R${read_num}.fq.gz"
		"${prefix}.R${read_num}.fastq"
		"${prefix}.R${read_num}.fq"
		"${prefix}.${read_num}.fastq.gz"
		"${prefix}.${read_num}.fq.gz"
		"${prefix}.${read_num}.fastq"
		"${prefix}.${read_num}.fq"
	)

	for pattern in "${patterns[@]}"; do
		local file="${dir}/${pattern}"
		if [[ -f "$file" && -r "$file" ]]; then
			echo "$file"
			return 0
		fi
	done

	# If no exact match, try wildcards (less preferred)
	local wildcard_patterns=(
		"${prefix}*[_.]R${read_num}*.f*q*"
		"${prefix}*[_.]${read_num}*.f*q*"
	)

	for pattern in "${wildcard_patterns[@]}"; do
		local matches=("${dir}"/"${pattern}")
		if [[ -f "${matches[0]}" && -r "${matches[0]}" ]]; then
			echo "${matches[0]}"
			return 0
		fi
	done

	return 1
}

# Count reads in a FASTQ file (handles both .fastq and .fastq.gz)
count_reads() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "0"
		return
	fi
	if [[ "$file" == *.gz ]]; then
		echo $(( $(gunzip -c "$file" | wc -l) / 4 ))
	else
		echo $(( $(wc -l < "$file") / 4 ))
	fi
}

echo "Searching for R1/R2 files for accession: $OUT in directory: $SEARCH_DIR"

R1=$(find_reads "$OUT" "$SEARCH_DIR" 1) || {
	echo "Error: Could not find R1 file for accession '$OUT' in $SEARCH_DIR" >&2
	echo "Looked for patterns like: ${OUT}_R1.fastq.gz, ${OUT}_1.fastq.gz, etc." >&2
	exit 1
}

R2=$(find_reads "$OUT" "$SEARCH_DIR" 2) || {
	echo "Error: Could not find R2 file for accession '$OUT' in $SEARCH_DIR" >&2
	echo "Looked for patterns like: ${OUT}_R2.fastq.gz, ${OUT}_2.fastq.gz, etc." >&2
	exit 1
}

echo "Found R1: $R1"
echo "Found R2: $R2"

# Count input reads
echo ""
echo "=== INPUT READ COUNTS ==="
R1_COUNT=$(count_reads "$R1")
R2_COUNT=$(count_reads "$R2")
TOTAL_INPUT=$((R1_COUNT + R2_COUNT))
echo "R1 reads: $R1_COUNT"
echo "R2 reads: $R2_COUNT"
echo "Total input reads: $TOTAL_INPUT"
echo ""

# ------------------------------- checks ----------------------------------- #
need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing dependency: $1" >&2
	exit 127
}; }
need repair.sh
need bbmerge.sh
need bbrename.sh
need reformat.sh
need gzip

# ------------------------------- paths ------------------------------------ #
FINAL="${OUT}.merged.fastq.gz"
IHIST="${OUT}_ihist_merge.txt"

tmpdir="$(mktemp -d -t repair_merge_tag.XXXXXX)"
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

# Repair.sh intermediate files
REPAIRED_R1="$tmpdir/repaired_R1.fastq.gz"
REPAIRED_R2="$tmpdir/repaired_R2.fastq.gz"
ORPHAN_READS="$tmpdir/orphans.fastq.gz"

# BBMerge intermediate files
MERGED_READY="$tmpdir/merged.tagged.ready.fastq.gz"
U1_TMP="$tmpdir/unmerged_R1.tagged.tmp.fastq.gz"
U2_TMP="$tmpdir/unmerged_R2.tagged.tmp.fastq.gz"
U1_READY="$tmpdir/unmerged_R1.tagged.ready.fastq.gz"
U2_READY="$tmpdir/unmerged_R2.tagged.ready.fastq.gz"
ORPHAN_READY="$tmpdir/orphans.tagged.ready.fastq.gz"

# ------------------------------ pipeline ---------------------------------- #
# 1) Run repair.sh to separate mated pairs from orphan reads
echo "Running repair.sh to separate mated pairs from orphans..."
repair.sh \
	in1="$R1" \
	in2="$R2" \
	out1="$REPAIRED_R1" \
	out2="$REPAIRED_R2" \
	outs="$ORPHAN_READS" \
	repair=t \
	threads="$THREADS" \
	-eoom

# Ensure repair outputs exist even if empty
for f in "$REPAIRED_R1" "$REPAIRED_R2" "$ORPHAN_READS"; do
	[[ -e "$f" ]] || gzip -c </dev/null >"$f"
done

# Count reads after repair
echo ""
echo "=== POST-REPAIR READ COUNTS ==="
REPAIRED_R1_COUNT=$(count_reads "$REPAIRED_R1")
REPAIRED_R2_COUNT=$(count_reads "$REPAIRED_R2")
ORPHAN_COUNT=$(count_reads "$ORPHAN_READS")
TOTAL_AFTER_REPAIR=$((REPAIRED_R1_COUNT + REPAIRED_R2_COUNT + ORPHAN_COUNT))
echo "Repaired R1 reads: $REPAIRED_R1_COUNT"
echo "Repaired R2 reads: $REPAIRED_R2_COUNT"
echo "Orphan reads: $ORPHAN_COUNT"
echo "Total after repair: $TOTAL_AFTER_REPAIR"
if [[ $TOTAL_AFTER_REPAIR -ne $TOTAL_INPUT ]]; then
	echo "WARNING: Lost $((TOTAL_INPUT - TOTAL_AFTER_REPAIR)) reads during repair!"
fi
echo ""

# 2) Merge repaired pairs; label streams in-flight.
echo "Merging paired reads with bbmerge..."
bbmerge.sh \
	in1="$REPAIRED_R1" \
	in2="$REPAIRED_R2" \
	out=>(bbrename.sh in=stdin out="$MERGED_READY" prefix=MERGED_) \
	outu1=>(bbrename.sh in=stdin out="$U1_TMP" prefix=UNMERGED_) \
	outu2=>(bbrename.sh in=stdin out="$U2_TMP" prefix=UNMERGED_) \
	qtrim=t \
	ihist="$IHIST" \
	threads="$THREADS" \
	-eoom

# Count reads after merging (before tagging orphans)
echo ""
echo "=== POST-MERGE READ COUNTS ==="
MERGED_COUNT=$(count_reads "$MERGED_READY")
UNMERGED_R1_COUNT=$(count_reads "$U1_TMP")
UNMERGED_R2_COUNT=$(count_reads "$U2_TMP")
TOTAL_FROM_PAIRS=$((MERGED_COUNT + UNMERGED_R1_COUNT + UNMERGED_R2_COUNT))
echo "Merged pairs → single reads: $MERGED_COUNT (reduced count by design)"
echo "Unmerged R1 reads: $UNMERGED_R1_COUNT"
echo "Unmerged R2 reads: $UNMERGED_R2_COUNT"
echo "Total from paired reads: $TOTAL_FROM_PAIRS"

# Check accounting - merged reads should reduce total by the number merged
EXPECTED_AFTER_MERGE=$((REPAIRED_R1_COUNT + REPAIRED_R2_COUNT - MERGED_COUNT))
if [[ $TOTAL_FROM_PAIRS -eq $EXPECTED_AFTER_MERGE ]]; then
	echo "✓ Read accounting correct (${MERGED_COUNT} pairs merged into singles)"
else
	DISCREPANCY=$((EXPECTED_AFTER_MERGE - TOTAL_FROM_PAIRS))
	echo "⚠ WARNING: Discrepancy of $DISCREPANCY reads during merging!"
fi
echo ""

# 3) Tag orphan reads with UNMERGED_ prefix (no /1 or /2 suffix forced)
echo "Tagging orphan reads..."
if [[ -s "$ORPHAN_READS" ]]; then
	bbrename.sh \
		in="$ORPHAN_READS" \
		out="$ORPHAN_READY" \
		prefix=UNMERGED_ \
		threads="$THREADS"
else
	gzip -c </dev/null >"$ORPHAN_READY"
fi

# 4) Ensure all merge outputs exist even if empty (valid empty gzip files).
for f in "$MERGED_READY" "$U1_TMP" "$U2_TMP" "$ORPHAN_READY"; do
	[[ -e "$f" ]] || gzip -c </dev/null >"$f"
done

# 5) Normalize unmerged names to /1 and /2 (handles empty inputs fine).
echo "Normalizing read names..."
reformat.sh \
	in1="$U1_TMP" \
	in2="$U2_TMP" \
	out1="$U1_READY" \
	out2="$U2_READY" \
	addslash=t spaceslash=f \
	threads="$THREADS"

# Count reads after normalization (final components before concatenation)
echo ""
echo "=== PRE-CONCATENATION COUNTS ==="
NORMALIZED_U1_COUNT=$(count_reads "$U1_READY")
NORMALIZED_U2_COUNT=$(count_reads "$U2_READY")
TAGGED_ORPHAN_COUNT=$(count_reads "$ORPHAN_READY")
echo "Final merged reads: $MERGED_COUNT"
echo "Final unmerged R1: $NORMALIZED_U1_COUNT"
echo "Final unmerged R2: $NORMALIZED_U2_COUNT"
echo "Final orphans: $TAGGED_ORPHAN_COUNT"
EXPECTED_FINAL=$((MERGED_COUNT + NORMALIZED_U1_COUNT + NORMALIZED_U2_COUNT + TAGGED_ORPHAN_COUNT))
echo "Expected total in final file: $EXPECTED_FINAL"
echo ""

# 6) Build final file in a deterministic order; single redirection for atomicity.
echo "Concatenating all reads into final output..."
cat \
	"$MERGED_READY" \
	"$U1_READY" \
	"$U2_READY" \
	"$ORPHAN_READY" \
	>"$FINAL"

# Final read count
echo ""
echo "=== FINAL OUTPUT ==="
FINAL_COUNT=$(count_reads "$FINAL")
echo "Total reads in $FINAL: $FINAL_COUNT"
echo ""

# Summary
echo "=== SUMMARY ==="
echo "Input reads: $TOTAL_INPUT (R1: $R1_COUNT, R2: $R2_COUNT)"
echo "Output reads: $FINAL_COUNT"

# Calculate expected reduction due to merging
EXPECTED_REDUCTION=$MERGED_COUNT
EXPECTED_OUTPUT=$((TOTAL_INPUT - EXPECTED_REDUCTION))

if [[ $FINAL_COUNT -eq $EXPECTED_OUTPUT ]]; then
	echo "✓ All reads accounted for!"
	echo "  $MERGED_COUNT read pairs were merged into single reads (expected reduction)"
elif [[ $FINAL_COUNT -lt $EXPECTED_OUTPUT ]]; then
	UNEXPECTED_LOSS=$((EXPECTED_OUTPUT - FINAL_COUNT))
	echo "⚠ WARNING: Lost $UNEXPECTED_LOSS reads beyond expected merging!"
	echo "  Expected $EXPECTED_OUTPUT reads after merging, but got $FINAL_COUNT"
else
	UNEXPECTED_GAIN=$((FINAL_COUNT - EXPECTED_OUTPUT))
	echo "⚠ WARNING: Have $UNEXPECTED_GAIN more reads than expected!"
	echo "  Expected $EXPECTED_OUTPUT reads after merging, but got $FINAL_COUNT"
fi

# Show merging efficiency
if [[ $REPAIRED_R1_COUNT -gt 0 ]]; then
	MERGE_RATE=$((MERGED_COUNT * 100 / REPAIRED_R1_COUNT))
	echo ""
	echo "Merge rate: ${MERGE_RATE}% of paired reads were successfully merged"
fi

# Breakdown of final file
echo ""
echo "=== FINAL FILE COMPOSITION ==="
echo "  Merged pairs (now single reads): $MERGED_COUNT"
echo "  Unmerged R1 reads: $UNMERGED_R1_COUNT"
echo "  Unmerged R2 reads: $UNMERGED_R2_COUNT"
echo "  Orphan/singleton reads: $ORPHAN_COUNT"
echo "  ────────────────────────────"
echo "  Total reads in output: $FINAL_COUNT"
echo ""

# Done.
echo "Wrote: $FINAL"
echo "Histogram: $IHIST"
echo "Pipeline complete."
