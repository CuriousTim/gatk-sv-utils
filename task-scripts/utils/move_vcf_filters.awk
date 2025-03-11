# Move BOTHSIDES_SUPPORT, HIGH_SR_BACKGROUND, and PESR_GT_OVERDISPERSION values
# from the FILTER field to the INFO field.
# Requires GNU awk
# usage: gawk -f move_vcf_filters.awk <vcf>
# If <vcf> is compressed, then the VCF must be piped and <vcf> should be '-'.

BEGIN {
	FS = "\t"
	OFS = "\t"
}

/^##fileformat/ {
	vcf_version = $0
	next
}

/^##contig/ {
	contigs[++ncontigs] = $0
	next
}

/^##FILTER=<ID=BOTHSIDES_SUPPORT,/ {
	next
}

/^##FILTER=<ID=HIGH_SR_BACKGROUND,/ {
	next
}

/^##FILTER=<ID=PESR_GT_OVERDISPERSION,/ {
	next
}

/^##INFO=<ID=BOTHSIDES_SUPPORT,/ {
	print "BOTHSIDES_SUPPORT is already in the INFO field" > "/dev/stderr"
	exit 1
}

/^##INFO=<ID=HIGH_SR_BACKGROUND,/ {
	print "HIGH_SR_BACKGROUND is already in the INFO field" > "/dev/stderr"
	exit 1
}

/^##INFO=<ID=PESR_GT_OVERDISPERSION,/ {
	print "PESR_GT_OVERDISPERSION is already in the INFO field" > "/dev/stderr"
	exit 1
}

/^##/ {
	other_headers[++nheaders] = $0
	next
}

/^#CHROM/ {
	other_headers[++nheaders] = "##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description=\"Variant has read-level support for both sides of breakpoint\">"
	other_headers[++nheaders] = "##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description=\"High number of SR splits in background samples indicating messy region\">"
	other_headers[++nheaders] = "##INFO=<ID=PESR_GT_OVERDISPERSION,Number=0,Type=Flag,Description=\"High PESR dispersion count\">"
	nheaders = asort(other_headers)
	print vcf_version
	for (i = 1; i <= nheaders; ++i) {
		print other_headers[i]
	}
	# we want to maintain the order of contig headers
	for (i = 1; i <= ncontigs; ++i) {
		print contigs[i]
	}
	print $0
	next
}

$7 ~ /BOTHSIDES_SUPPORT/ || $7 ~ /HIGH_SR_BACKGROUND/ || $7 ~ /PESR_GT_OVERDISPERSION/ {
	split($7, filters, /;/)
	for (i in filters) {
		if (filters[i] == "BOTHSIDES_SUPPORT") {
			$8 = $8 ";BOTHSIDES_SUPPORT"
			delete filters[i]
		} else if (filters[i] == "HIGH_SR_BACKGROUND") {
			$8 = $8 ";HIGH_SR_BACKGROUND"
			delete filters[i]
		} else if (filters[i] == "PESR_GT_OVERDISPERSION") {
			$8 = $8 ";PESR_GT_OVERDISPERSION"
			delete filters[i]
		}
	}

	new_filter = ""
	for (i in filters) {
		new_filter = new_filter ";" filters[i]
	}
	sub(/^;/, "", new_filter)

	if (length(new_filter) == 0) {
		new_filter = "PASS"
	}
	$7 = new_filter
}

1
