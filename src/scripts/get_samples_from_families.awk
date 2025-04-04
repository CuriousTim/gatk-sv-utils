# usage: gawk -f <script> [--families <families> | --nfamilies <n>]
#                         <pedigree>
#
# Print family members from <pedigree>.
#
# options:
#   --families <families>  Get family IDs from file <families>. Takes
#                          precedence over --nfamilies.
#   --nfamilies <n>        Select up to <n> random families in <pedigree>.

@include "logging"
@include "random"
@include "sh"

BEGIN {
	init()

	pos = 0
	nargs = 0
	while (++pos < ARGC) {
		if (ARGV[pos] == "--families") {
			Options["families_file"] = ARGV[++pos]
			continue
		}

		if (ARGV[pos] == "--nfamilies") {
			Options["nfamilies"] = int(ARGV[++pos] + 0)
			continue
		}

		args[++nargs] = ARGV[pos]
	}

	Pedigree = args[1]
	logging::log_info("using pedigree at " Pedigree)
	if (Options["families_file"]) {
		logging::log_info("using families list at " Options["families_file"])
		ARGV[1] = Options["families_file"]
		ARGV[2] = Pedigree
	} else {
		logging::log_info(Options["nfamilies"] " random families requested")
		if (Options["nfamilies"] <= 0) {
			logging::log_err("number of families to select must be > 0")
			exit 1
		}
		ARGV[1] = Pedigree
		ARGV[2] = Pedigree
	}

	ARGC = 3
}

ENDFILE {
	if (FILENAME == Options["families_file"]) {
		if (length(Fams) == 0) {
			logging::log_err("no families found in family list")
			exit 1
		} else {
			logging::log_info("read " length(Fams) " families from family list")
		}
	}

	if (ARGIND == 1 && FILENAME == Pedigree) {
		if (length(PedFams) == 0) {
			logging::log_err("no families found in pedigree")
			exit 1
		} else {
			logging::log_info("read " length(PedFams) " families from pedigree")
			select_random_families()
		}
	}

	if (ARGIND == 1) {
		logging::log_info("printing samples")
	}
}

ARGIND == 1 && FILENAME == Options["families_file"] {
	Fams[$1]
}

ARGIND == 1 && FILENAME == Pedigree && !/^#/ {
	PedFams[$1]
}

ARGIND == 2 && !/^#/ && ($1 in Fams) {
	Did_print = 1
	print $2
}

END {
	if (Did_print) {
		logging::log_info("done")
	}
}

function select_random_families(    tmp, i) {
	random::sample_array(Options["nfamilies"], PedFams, tmp)
	delete Fams
	for (i in tmp) {
		Fams[tmp[i]]
	}
}

function init() {
	FS = "\t"
	OFS = "\t"
	Options["families_file"] = ""
	Options["nfamilies"] = 1000
	srand()
}
