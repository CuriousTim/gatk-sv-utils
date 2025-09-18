# usage: gawk -f set_vcf_svlens.awk <vcf>
#
# Set the SVLEN INFO field of every site in a VCF according to SVTYPE.
# * DEL, DUP, and INV: END - POS + 1
# * INS: 50
# * CPX, CTX, and CNV: ignore
#
# This is to satisfy the SVLEN requirement of SVConcordance, but with the SVLEN
# set this way, the minimum size similarity for INS must be set to 0 with a
# clustering config otherwise the matching will be inaccurate SV lengths for
# INS. Any existing SVLEN will be overwritten.
#
# The file must be piped to the script as an uncompressed VCF.

@include "logging"

BEGIN {
	FS = "\t"
	OFS = "\t"
}

/^##fileformat/ {
	Version = $0
	next
}

/^##contig/ {
	Contigs[++Ncontigs] = $0
	next
}

/^##INFO=<ID=SVLEN,/ {
	next
}

/^##/ {
	Headers[++Nheaders] = $0
	next
}

/^#CHROM/ {
	Headers[++Nheaders] = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
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

{
	match($8, /SVTYPE=([^;]*)/, a)
	match($8, /END=([^;]*)/, b)
	if (!RSTART) {
		logging::log_warn("record with ID '" $3 "' does not have INFO/END")
		next
	}

	if (a[1] == "DEL" || a[1] == "DUP" || a[1] == "INV")
		svlen = (b[1] + 0) - $2 + 1
	else if (a[1] == "INS")
		svlen = 50
	else
		next

	if ($8 ~ /SVLEN=[^;]*/)
		sub(/SVLEN=[^;]*/, "SVLEN=" svlen, $8)
	else
		$8 = $8 ";SVLEN=" svlen
	print
}

END {
	logging::log_info("done: " Records " processed")
}
