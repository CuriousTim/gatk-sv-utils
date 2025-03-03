# Add END2 to INFO field of BND and CTX variants in a VCF.
# The VCF must be streamed e.g.
# bcftools view --output-type v | awk -f add_end2_to_vcf.awk

BEGIN {
	FS = "\t"
	OFS = "\t"
}

/^##/ {
	print
	if ($0 ~ /##INFO=<ID=END2/) {
		header_seen = 1
	}
	next
}

# add a meta line for END2 if it didn't exist
/^#CHROM/ {
	if (!header_seen) {
		print "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Position of breakpoint on CHR2\">"
	}
	print
	next
}

# column 8 is INFO column
# don't overwrite an existing END2
$8 ~ /END2=/ {
	print
	next
}

# column 5 is ALT column
$5 ~ /BND|CTX/ {
	chr2 = ""
	svlen = 0
	match($8, /CHR2=[^;]+/)
	if (RSTART) {
		chr2 = substr($8, RSTART + 5, RLENGTH - 5)
	}

	# column 1 is CONTIG column
	# POS + SVLEN is only valid for intrachromosomal events
	if ($1 != chr2) {
		print
		next
	}

	match($8, /SVLEN=[^;]+/)
	if (RSTART) {
		svlen = substr($8, RSTART + 6, RLENGTH - 6)
	} else {
		# missing SVLEN is probably a bug
		print
		next
	}

	# column 2 is POS column
	$8 = $8 ";END2=" ($2 + svlen)
}

1
