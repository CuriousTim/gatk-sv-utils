# usage: gawk -f <script> <thresholds> <vcf>
#
# Filter the genotypes in <vcf> by setting them to null for the samples with a
# GQ lower than the threshold. Any records that have all genotypes updated to
# missing are not printed.
#
# <thresholds> is a tab-delimited file with an SV type in the first column and
# a minimum GQ in the second column. Each SV type will be filtered using the
# corresponding GQ threshold.

@include "logging"
@include "vcf"

BEGIN {
	FS = "\t"
	OFS = "\t"

	Record_n = 0
	Gts_changed = 0
}

ENDFILE {
	if (ARGIND == 1 && length(Thresholds) == 0) {
		logging::log_err("no thresholds given")
		exit 1
	}
}

END {
	logging::log_info(Gts_changed " genotypes updated")
}

ARGIND == 1 {
	if (NF != 2) {
		logging::log_err("expected 2 columns in thresholds file, got " NF " [line: " FNR "]")
		exit 1
	}

	Thresholds[$1] = $2 + 0
}

ARGIND == 2 && $0 ~ /^#/ {
	print
	next
}

ARGIND == 2 {
	if (Record_n > 0 && Record_n % 1000 == 0) {
		logging::log_info("processed " Record_n " records")
	}

	++Record_n

	vcf::parse_info($8, info)
	if (!("SVTYPE" in info)) {
		logging::log_warn("record " Record_n " has no SVTYPE")
	}

	if (!(info["SVTYPE"] in Thresholds)) {
		print
		next
	}

	min_gq = Thresholds[info["SVTYPE"]]
	ncc = 0
	for (i = 10; i <= NF; ++i) {
		vcf::parse_format($9, $i, gt)
		if (!("GQ" in gt["gt"])) {
			logging::log_warn("record " Record_n " missing GQ")
			print
			next
		}

		if (!("GT" in gt["gt"])) {
			logging::log_warn("record " Record_n " missing GT")
			print
			next
		}

		if (gt["gt"]["GQ"] != "." && (gt["gt"]["GQ"] + 0) < min_gq) {
			gt["gt"]["GT"] = "./."
			$i = vcf::make_format(gt)
			++Gts_changed
		}

		if (gt["gt"]["GT"] == "./.") {
			++ncc
		}
	}

	# Only print the record if there is at least one non-missing GT
	if (ncc < NF - 9) {
		print
	} else {
		logging::log_warn("record " Record_n " has no genotypes after update")
	}
}
