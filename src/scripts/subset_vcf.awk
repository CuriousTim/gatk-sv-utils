# usage: gawk -f <script> [--samples <samples> | --nsamples <n>]
#                         [--keep-private-sites] [--update-af]
#                         <vcf> <outfile>
#
# Subset <vcf> with a list of samples and write to <outfile>.
#
# options:
#   --samples <samples>   Get sample IDs from file <samples>. Takes
#                         precedence over --nsamples.
#   --nsamples <n>        Select up to <n> random samples in <vcf>.
#   --keep-private-sites  Keep sites in VCF that do not have any non-ref
#                         genotypes after subsetting. Default is to remove.
#   --update-af           Update site INFO/AC and INFO/AN after subsetting.
#                         Default is to leave them.

@include "logging"
@include "random"
@include "sh"

BEGIN {
	init()

	pos = 0
	nargs = 0
	while (++pos < ARGC) {
		if (ARGV[pos] == "--samples") {
			Options["samples_file"] = ARGV[++pos]
			continue
		}

		if (ARGV[pos] == "--nsamples") {
			Options["nsamples"] = int(ARGV[++pos] + 0)
			continue
		}

		if (ARGV[pos] == "--keep-private-sites") {
			Options["rm_private_sites"] = 0
			continue
		}

		if (ARGV[pos] == "--update-af") {
			Options["keep_af"] = 0
			continue
		}

		args[++nargs] = ARGV[pos]
	}

	vcf = args[1]
	outfile = args[2]

	logging::log_info("subsetting VCF at " vcf)
	logging::log_info("writing subset VCF to " outfile)
	logging::log_info("sites without any non-ref genotypes after subsetting will be " (Options["rm_private_sites"] ? "removed" : "kept"))
	logging::log_info("INFO/AC and INFO/AN will " (Options["keep_af"] ? "not" : "") " be recomputed after subsetting")

	if (Options["samples_file"]) {
		logging::log_info("using samples list at " Options["samples_file"])
	} else {
		logging::log_info(Options["nsamples"] " random samples requested")
		if (Options["nsamples"] <= 0) {
			logging::log_err("number of samples to select must be > 0")
			exit 1
		}

		select_random_samples(vcf)
		if (length(Samples) > 0) {
			logging::log_info("selected " length(Samples) " samples from VCF")
		} else {
			logging::log_err("no samples found in VCF")
		}

		Tmpfile = mktemp()
		if (!Tmpfile) {
			logging::log_err("failed to create temporary file")
			exit 1
		}

		write_samples(Tmpfile)
		Options["samples_file"] = Tmpfile
	}

	rtn = subset_vcf(vcf, outfile, Options["samples_file"])
	if (rtn != 0) {
		logging::log_err("subsetting VCF returned non-zero exit code")
		exit 1
	} else {
		logging::log_info("subsetting VCF done")
	}

	logging::log_info("indexing subset VCF")
	rtn = system("bcftools index --tbi " sh::quote(outfile))
	if (rtn != 0) {
		logging::log_err("indexing VCF returned non-zero exit code")
		exit 1
	}

	logging::log_info("done")

	exit 0
}

END {
	if (Tmpfile) {
		system("rm -f " sh::quote(Tmpfile))
	}
}

function select_random_samples(vcf,    cmd, i, j, tmp, vcf_samples) {
	delete Samples
	# This relies on the list of samples from the VCF not having
	# duplicates, which should be the case.
	cmd = "bcftools query --list-samples " sh::quote(vcf)
	while ((cmd | getline line) > 0) {
		vcf_samples[++i] = line
	}
	close(cmd)

	random::sample(Options["nsamples"], i, tmp)
	i = 0
	for (j in tmp) {
		Samples[++i] = vcf_samples[tmp[j]]
	}
}

function subset_vcf(vcf, outfile, samples_file,    cmd) {
	cmd = "bcftools view --force-samples --samples-file " sh::quote(samples_file)
	if (Options["keep_af"]) {
		cmd = cmd " --no-update " sh::quote(vcf)
	}

	if (Options["rm_private_sites"]) {
		cmd = cmd " --output-type u"
		cmd = cmd " | bcftools view --output-type v --exclude 'SVTYPE != \"CNV\" && COUNT(GT = \"alt\") == 0'"
	} else {
		cmd = cmd " --output-type v"
	}

	cmd = cmd " | grep -v '^##bcftools' | bgzip -c > " sh::quote(outfile)
	print cmd

	return system(cmd)
}

function write_samples(file) {
	for (i in Samples) {
		print Samples[i] > file
	}
	close(file)
}

function mktemp(    cmd, tmp_path) {
	cmd = "mktemp -p " sh::quote(ENVIRON["PWD"])
	cmd | getline tmp_path
	close(cmd)

	return tmp_path
}

function init() {
	FS = "\t"
	OFS = "\t"
	Options["samples_files"] = ""
	Options["nsamples"] = 1000
	Options["rm_private_sites"] = 1
	Options["keep_af"] = 1
	srand()
}
