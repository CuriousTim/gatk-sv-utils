# usage: gawk -f add_end2_to_vcf.awk <vcf>
#
# Add END2 to the INFO field of intrachromosomal BND and CTX records in a VCF.
# For these SVs, END2 is defined as POS + INFO/SVLEN.
#
# This script can only process text so if the VCF is compressed, it must be
# piped as uncompressed e.g.
# bcftools view --output-type v | gawk -f add_end2_to_vcf.awk -
# The updated VCF will be written to stdout.

@include "logging"

BEGIN {
	FS = "\t"
	OFS = "\t"
	logging::log_info("adding END2")
	Records = 0
}

/^##INFO=<ID=END2,/ {
	if ($0 ~ /Description="Position of breakpoint on CHR2"/) {
		logging::log_err("VCF contains an END2 header with a different description: ")
		print "    " $0 > "/dev/stderr"
		exit 1
	}

	next
}

/^##fileformat/ {
	Version = $0
	next
}

/^##contig/ {
	Contigs[++Ncontigs] = $0
	next
}

/^##/ {
	Headers[++Nheaders] = $0
	next
}

/^#CHROM/ {
	Headers[++Nheaders] = "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Position of breakpoint on CHR2\">"
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

# column 8 is INFO column
# don't overwrite an existing END2
$8 ~ /^END2=/ || $8 ~ /;END2=/ {
	print
	next
}

# column 5 is ALT column
$5 ~ /BND|CTX/ {
	chr2 = ""
	svlen = 0
	match($8, /CHR2=([^;]+)/, a)
	if (RSTART) {
		chr2 = a[1]
	}

	# column 1 is CONTIG column
	# POS + SVLEN is only valid for intrachromosomal events
	if ($1 != chr2) {
		print
		next
	}

	match($8, /SVLEN=([^;]+)/, a)
	if (RSTART) {
		svlen = a[1]
	} else {
		# missing SVLEN is probably a bug
		print
		next
	}

	# column 2 is POS column
	$8 = $8 ";END2=" ($2 + svlen)
	sub(/^;/, "", $8)
}

1

END {
	logging::log_info("done: " Records " processed")
}
