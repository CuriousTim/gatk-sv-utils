# usage: gawk -f set_gt_null.awk <contigs> <samples> <vcf>
#
# Set genotypes on given contigs and samples to missing in a VCF.
#
# This script can only process text so if the VCF is compressed, it must be
# piped as uncompressed e.g.
# bcftools view --output-type v \
#   | gawk -f set_gt_to_null.awk <contigs> <samples> -
# The updated VCF will be written to stdout.

@include "logging"

BEGIN {
	FS = "\t"
	OFS = "\t"
}

# Read contigs
ARGIND == 1 && $0 {
	++Ncontigs
	Contigs[$1]
	next
}

# Read samples
ARGIND == 2 && $0 {
	++Nsamples
	Samples[$1]
	next
}

# Print headers
ARGIND == 3 && /^##/ {
	print
	next
}

ARGIND == 3 && /^#CHROM/ {
	# Sample IDs are from column 10 and on
	if (NF < 10) {
		logging::log_err("VCF does not have samples")
		exit 1
	}

	for (i = 10; i <= NF; ++i) {
		Vcf_samples_map[$i] = i
	}

	print
	next
}

# If no contigs are given, null the genotypes on all contigs.
ARGIND == 3 && Ncontigs > 0 && !($1 in Contigs) {
	print
	next
}

{
	format = $9
	split(format, a, /:/)
	# VCF spec says GT must always be first key, if present
	if (a[1] != "GT") {
		msg = sprintf("record on line %d does not have genotypes (FORMAT = %s)",
			FNR, format)
		logging::log_warn(msg)
		print
		next
	}

	if (NF < 10) {
		logging::log_warn("record on line " FNR " does not have samples")
		print
		next
	}

	for (s in Samples) {
		if (s in Vcf_samples_map) {
			i = Vcf_samples_map[s]
			split($i, a, /:/)
			gsub(/[0-9]+/, ".", a[1])
			$i = join_fmt(a)
		}
	}

	print
}

function join_fmt(a,    s, i) {
	s = ""
	for (i in a) {
		s = s a[i] ":"
	}
	sub(/:$/, "", s)

	return s
}
