# usage: gawk -f set_vcf_algorithms.awk <vcf>
#
# Set the ALGORITHMS INFO field of every site in a VCF to pesr so that
# SVConcordance will work.
#
# The file must be piped to the script as an uncompressed VCF.

@include "logging"

BEGIN {
	FS = "\t"
	OFS = "\t"
	New_algo = "ALGORITHMS=pesr"
}

/^##fileformat/ {
	Version = $0
	next
}

/^##contig/ {
	Contigs[++Ncontigs] = $0
	next
}

/^##INFO=<ID=ALGORITHMS,/ {
	Algo_seen = 1
	Headers[++Nheaders] = $0
	next
}

/^##/ {
	Headers[++Nheaders] = $0
	next
}

/^#CHROM/ {
	if (!Algo_seen)
		Headers[++Nheaders] = "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">"
	asort(Headers)
	print Version
	for (i in Headers)
		print Headers[i]

	for (i in Contigs)
		print Contigs[i]
	print
	next
}

{
	++Records
}

Records % 1000 == 0 {
	logging::log_info(Records " records processed")
}

$8 ~ /ALGORITHMS=[^;]*/ {
	sub(/ALGORITHMS=[^;]*/, New_algo, $8)
	print
	next
}

$8 == "." {
	$8 = New_algo
	print
	next
}

{
	$8 = $8 ";" New_algo
	print
}

END {
	logging::log_info("done: " Records " processed")
}
