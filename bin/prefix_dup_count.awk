#!/usr/bin/env -S awk -f

BEGIN {
	OFS = "\t"
}


# Pass through header lines unchanged
/^@/ {
	print
	next
}

{
	count = ""
	for (i = 12; i <= NF; i++) {
		if ($i ~ /^dc:i:/) {
			split($i, a, ":")
			count = a[3]
			break
		}
	}
	if (count != "") {
		$1 = $1 "-" count
	}
	print
}

# For each alignment, if a dc:i: tag exists, append it to QNAME for SAM Refiner to be able to get counts
