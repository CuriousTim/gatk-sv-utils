# Benchmark the de novo pipeline.
# In order to run this script, SVConcordance needs to be run three times with
# the `--keep-all` option, which reports all matching sites in the truth VCF
# for each site in the evaluation VCF.
#
# SVConcordance Runs
# ==================
#
# eval vs truth
# -------------
# The evaluation VCF is the de novo VCF (results from the de novo pipeline) and
# the truth VCF is the true de novos VCF.
#
# truth vs eval
# -------------
# The evaluation VCF is the true de novos VCF and the truth VCF is the de novo
# VCF.
#
# truth vs start
# --------------
# The evaluation VCF is the true de novos VCF and the truth VCF is the start
# VCF (the input to the de novo pipeline).
#
# --keep-all
# ==========
# SVConcordance with the `--keep-all` option, will annotate each site in the
# evalution VCF with a comma-separated list of the IDs of all the sites in the
# truth VCF that match the site in the evaluation VCF. To check that a sample
# in the evaluation VCF has a matching SV in the truth VCF, the genotypes from
# all matching sites much be checked.
#
# usage: gawk -f benchmark_denovo.awk <eval_v_truth> <truth_v_start>
#   <truth_v_eval> <start>
#
# <eval_v_truth>: a TSV dump of the de novos VCF after (eval vs truth)
#   1. chromosome
#   2. start
#   3. end
#   4. SV type
#   5. site ID
#   6. comma-separated IDs of matching sites in the truth VCF, or "." for
#      missing
#   7+ all remaning columns are IDs of samples that carry the SV
#
# <truth_v_start>: a TSV dump of the true de novos VCF after (truth vs start)
#   1. site ID
#   2. comma-separated IDs of matching sites in the start VCF, or "." for
#      missing
#
# <truth_v_eval>: a TSV dump of the true de novos VCF after (truth vs eval) in
#   the same format as <eval_v_truth>
#
# <start>: a TSV dump of the start VCF
#   1. site ID
#   2+ all remaining columns are IDs of samples that carry the SV
# For memory efficiency, the site IDs should be restricted to those that appear
# in the second column of <truth>
#
# Five output files will be written.
#
# eval_vs_truth.tsv
# -----------------
# The first five columns of <eval_v_truth> plus:
# * count of true positives at that site
# * count of false positives at that site
#
# truth_vs_eval.tsv
# -----------------
# The first five columns of <truth_v_eval> plus:
# * count of genotypes in the truth VCF that matched a genotype in the start
#   VCF, but not the eval VCF (false negative type 1)
# * count of genotypes in the truth VCF that did not match a genotype in the
#   start VCF (false negative type 2)
#
# false_positives.tsv
# -------------------
# All false positive de novo calls.
#
# false_negatives-type1.tsv
# -------------------------
# All type 1 false negative de novo calls.
#
# false_negatives-type2.tsv
# -------------------------
# All type 2 false negative de novo calls.

BEGIN {
	FS = "\t"
	OFS = "\t"

	read_gts(ARGV[1], Eval_gts, 5, 7)
	read_gts(ARGV[3], Truth_gts, 5, 7)
	read_gts(ARGV[4], Start_gts, 1, 2)

	--ARGC

	print "chr", "start", "end", "svtype", "vid", "samples" > "false_positive.tsv"
	print "chr", "start", "end", "svtype", "vid", "samples" > "false_negatives-type1.tsv"
	print "chr", "start", "end", "svtype", "vid", "samples" > "false_negatives-type2.tsv"
}

# check eval vs truth
ARGIND == 1 {
	tp = 0
	fp = 0
	# this site in the eval VCF did not have any matching variants in the
	# truth VCF
	if ($5 == ".") {
		fp += NF - 6
		for (i = 7; i <= NF; ++i) {
			print $1, $2, $3, $4, $5, $i > "false_positives.tsv"
		}
	} else {
		split($6, vids, /,/)
		for (i = 7; i <= NF; ++i) {
			if (has_match($i, vids, Truth_gts)) {
				++tp
			} else {
				++fp
				print $1, $2, $3, $4, $5, $i > "false_positives.tsv"
			}
		}
	}

	print $1, $2, $3, $4, $5, tp, fp > "eval_vs_truth.tsv"
}

# read truth VCF sites with matching start VCF sites
ARGIND == 2 && $2 != "." {
	split($2, vids, /,/)
	for (i in vids)
		Truth_in_start[$1][i] = vids[i]
}

# check truth vs eval
ARGIND == 3 {
	fn1 = 0
	fn2 = 0
	if ($6 != ".")
		split($6, vids, /,/)

	for (i = 7; i <= NF; ++i) {
		# check if the truth genotype is in the start VCF first because
		# the eval VCF is a subset of the start VCF so it cannot be in
		# the eval VCF, but not in the start VCF
		if (!($5 in Truth_in_start) || !has_match($i, Truth_in_start[$5], Start_gts)) {
			++fn2
			print $1, $2, $3, $4, $5, $i > "false_negatives-type2.tsv"
			continue
		}

		if ($6 == "." || !has_match($i, vids, Eval_gts)) {
			print $1, $2, $3, $4, $5, $i > "false_negatives-type1.tsv"
			++fn1
		}
	}

	print $1, $2, $3, $4, $5, fn1, fn2 > "truth_vs_eval.tsv"
}

function read_gts(file, arr, vid_col, sid_col,    line, n, vid, fields, i) {
	while ((getline line < file) > 0) {
		n = split(line, fields, /\t/)
		vid = fields[vid_col]
		for (i = sid_col; i <= n; ++i)
			arr[vid][fields[i]]
	}
	close(file)
}

# Does a sample with sid have an alt genotype at any of the variants given in
# vids?
function has_match(sid, vids, gts,    i) {
	for (i in vids) {
		if (vids[i] in gts && sid in gts[vids[i]])
			return 1
	}

	return 0
}
