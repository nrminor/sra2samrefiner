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
		$1 = count "-" $1
	}
	print
}

# For each alignment, if a dc:i: tag exists, prepend it to QNAME
