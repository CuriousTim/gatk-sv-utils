# Set genotypes on given contigs and samples to missing in a VCF.
# Requires GNU awk.
# usage: gawk -f set_gt_null.awk <contigs> <samples> <vcf>
# If <vcf> is compressed, then the VCF must be piped and <vcf> should be '-'.
# The script does not check that the VCF is well-formed, so if there is an
# error in the VCF, the output will likely be garbage.

BEGIN {
	FS = "\t"
	OFS = "\t"
	ncontigs = 0
}

# Read contigs
ARGIND == 1 {
	++ncontigs
	contigs_to_null[$1]
	next
}

# Read samples
ARGIND == 2 {
	samples_to_null[$1]
	next
}

# Print headers
ARGIND == 3 && /^##/ {
	print
	next
}

ARGIND == 3 && /^#CHROM/ {
	if (NF < 10) {
		print "VCF does not have samples" > "/dev/stderr"
		exit 1
	}

	for (i = 10; i <= NF; ++i) {
		samples_map[$i] = i
	}

	print
	next
}

# If no contigs are given, null the genotypes on all contigs.
ARGIND == 3 && ncontigs > 0 && !($1 in contigs_to_null) {
	print
	next
}

{
	format = $9
	split(format, a, /:/)
	# VCF spec says GT must always be first key, if present
	if (a[1] != "GT") {
		printf "warning: record on line %d does not have genotypes (FORMAT = %s)",
		       FNR, format > "/dev/stderr"
		print
		next
	}

	if (NF < 10) {
		print "warning: record on line %d does not have samples", FNR > "/dev/stderr"
		print
		next
	}

	for (s in samples_to_null) {
		if (s in samples_map) {
			i = samples_map[s]
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
