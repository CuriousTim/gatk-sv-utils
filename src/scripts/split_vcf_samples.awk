# usage: gawk -f <script> <vcf> <genotypes_per_split> <split_dir>
#
# Split the samples in a VCF so that each split will have approximately
# <genotypes_per_split> genotypes. The list of samples in each split will
# placed in <split_dir>.

@include "logging"
@include "sh"

BEGIN {
	FS = "\t"

	Vcf = ARGV[1]
	Gts_per_split = ARGV[2] + 0
	Splits_dir = ARGV[3]

	if (Gts_per_split <= 0) {
		logging::log_err("genotypes per split must be positive. got '" Gts_per_split "'")
		exit 1
	}

	records_n = count_vcf_records(Vcf)
	if (records_n == 0) {
		logging::log_err("no records in VCF")
		exit 1
	} else {
		logging::log_info("found " records_n " records in VCF")
	}

	samples_n = get_vcf_samples(Vcf, vcf_samples)
	if (samples_n == 0) {
		logging::log_err("no samples in VCF")
		exit 1
	} else {
		logging::log_info("found " samples_n " samples in VCF")
	}

	if (Gts_per_split <= records_n) {
		logging::log_warn("genotypes per split <= number of records")
		samples_per_split = 1
	} else {
		samples_per_split = int(Gts_per_split / records_n)
		remaining_samples = samples_n % samples_per_split
		excess_genotypes = (records_n * samples_n) - (records_n * remaining_samples)
		# We want to avoid having a split with a lot fewer genotypes
		# than the other ones.
		if ((4 * excess_genotypes) <= Gts_per_split) {
			full_splits = int(samples_n / samples_per_split)
			samples_per_split += int(remaining_samples / full_splits) + 1
		}
	}

	mkdir(Splits_dir)
	Made_splits_dir = 1
	splits_n = int(samples_n / samples_per_split) + (samples_n % samples_per_split > 0)
	logging::log_info("using " samples_per_split " samples per split")
	prefix_width = length(splits_n "")
	j = 1
	out = sprintf(Splits_dir "/part-%0" prefix_width "d.list", j)
	for (i = 1; i <= samples_n; i++) {
		print vcf_samples[i] > out
		if (i % samples_per_split == 0) {
			close(out)
			++j
			out = sprintf(Splits_dir "/part-%0" prefix_width "d.list", j)
		}
	}
	Wrote_splits = 1

	exit 0
}

END {
	if (!Wrote_splits && Made_splits_dir) {
		rmdir(Splits_dir)
	}
}

function count_vcf_records(vcf,    cmd, count, line, fields) {
	cmd = "bcftools index --stats " sh::quote(vcf)
	count = 0
	while ((cmd | getline line) > 0) {
		n = split(line, fields)
		if (n != 3) {
			logging::log_err("'bcftools index --stats' did not produce three fields")
			exit 1
		}
		count += fields[3] + 0
	}
	close(cmd)

	return count
}

function get_vcf_samples(vcf, samples,    cmd, i, line) {
	cmd = "bcftools query --list-samples " sh::quote(vcf)
	i = 0
	while ((cmd | getline line) > 0) {
		samples[++i] = line
	}
	close(cmd)

	return i
}

function mkdir(x,    rtn) {
	rtn = system("mkdir " sh::quote(x))
	if (rtn != 0) {
		logging::log_err("failed to create directory '" x "'")
		exit 1
	}
}

function rmdir(x) {
	system("rm -rf " sh::quote(x))
}
