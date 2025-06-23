# usage: gawk -f move_vcf_filters.awk <vcf>
#
# Move BOTHSIDES_SUPPORT, HIGH_SR_BACKGROUND, and PESR_GT_OVERDISPERSION values
# from the FILTER field to the INFO field.
#
# This script can only process text so if the VCF is compressed, it must be
# piped as uncompressed e.g.
# bcftools view --output-type v | gawk -f move_vcf_filters.awk -
# The updated VCF will be written to stdout.

@include "logging"

BEGIN {
	FS = "\t"
	OFS = "\t"
	logging::log_info("moving VCF filters")
	Records = 0
}

/^##fileformat/ {
	Version = $0
	next
}

/^##contig/ {
	Contigs[++Ncontigs] = $0
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
	logging::log_err("BOTHSIDES_SUPPORT is already in the INFO field")
	exit 1
}

/^##INFO=<ID=HIGH_SR_BACKGROUND,/ {
	logging::log_err("HIGH_SR_BACKGROUND is already in the INFO field")
	exit 1
}

/^##INFO=<ID=PESR_GT_OVERDISPERSION,/ {
	logging::log_err("PESR_GT_OVERDISPERSION is already in the INFO field")
	exit 1
}

/^##/ {
	Headers[++Nheaders] = $0
	next
}

/^#CHROM/ {
	Headers[++Nheaders] = "##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description=\"Variant has read-level support for both sides of breakpoint\">"
	Headers[++Nheaders] = "##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description=\"High number of SR splits in background samples indicating messy region\">"
	Headers[++Nheaders] = "##INFO=<ID=PESR_GT_OVERDISPERSION,Number=0,Type=Flag,Description=\"High PESR dispersion count\">"
	asort(Headers)
	print Version
	for (i in Headers) {
		print Headers[i]
	}

	# don't sort contigs to maintain their order
	for (i in Contigs) {
		print Contigs[i]
	}
	print
	next
}

{
	++Records
}

Records % 1000 == 0 {
	logging::log_info(Records " records processed")
}

$7 ~ /BOTHSIDES_SUPPORT/ || $7 ~ /HIGH_SR_BACKGROUND/ || $7 ~ /PESR_GT_OVERDISPERSION/ {
	split($7, filters, /;/)
	for (i in filters) {
		if (filters[i] == "BOTHSIDES_SUPPORT") {
			$8 = $8 ";BOTHSIDES_SUPPORT"
			delete filters[i]
		} else if (filters[i] == "HIGH_SR_BACKGROUND") {
			$8 = $8  ";HIGH_SR_BACKGROUND"
			delete filters[i]
		} else if (filters[i] == "PESR_GT_OVERDISPERSION") {
			$8 = $8 ";PESR_GT_OVERDISPERSION"
			delete filters[i]
		}
	}
	sub(/^;/, "", $8)

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

END {
	logging::log_info("done: " Records " processed")
}
