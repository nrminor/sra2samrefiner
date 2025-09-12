#!/usr/bin/env awk -f

BEGIN {
	OFS = "\t"
}


# Pass through header lines unchanged
/^@/ {
	print
	next
}

{
	dc = ""
	for (i = 12; i <= NF; i++) {
		if ($i ~ /^dc:i:/) {
			split($i, a, ":")
			dc = a[3]
			break
		}
	}
	if (dc != "") {
		$1 = dc "-" $1
	}
	print
}

# For each alignment, if a dc:i: tag exists, prepend it to QNAME
