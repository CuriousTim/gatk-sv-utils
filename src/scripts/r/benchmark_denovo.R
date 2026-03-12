library(data.table)

denovo_calls <- fread("denovo_svs-outliers_flagged.tsv.gz")
eval_vs_truth_hits <- fread(
    "eval_in_truth.tsv.gz",
    header = FALSE,
    col.names = c("vid", "match_vid"),
    key = "vid"
)
truth_carriers <- fread(
    "truth_vcf_carriers.tsv.gz",
    header = FALSE,
    col.names = c(
        "chr", "start", "end", "svtype",
        "vid", "sample"
    ),
    key = "vid"
)
truth_vs_start_hits <- fread(
    "truth_in_start.tsv.gz",
    header = FALSE,
    col.names = c("vid", "match_vid"),
    key = "vid"
)
start_carriers <- fread(
    "start_vcf_carriers.tsv.gz",
    header = FALSE,
    col.names = c("start_vid", "samples"),
    key = "start_vid"
)
start_carriers[, start_samples := strsplit(samples, ",", fixed = TRUE)]
start_carriers[, samples := NULL]

truth_vcf_annotations <- fread(
    "truth_vcf_annotations.tsv",
    header = FALSE,
    col.names = c("vid", "genomic_context", "ovp_pc_gene"),
    key = "vid"
)

eval_benchmark <- eval_vs_truth_hits[denovo_calls, on = c(vid = "name")]
truth_carriers_by_vid <- truth_carriers[, list(truth_carriers = list(sample)), by = "vid"]
setnames(truth_carriers_by_vid, "vid", "truth_vid")
eval_benchmark <- truth_carriers_by_vid[eval_benchmark, on = c(truth_vid = "match_vid")]
eval_benchmark[, is_true_de_novo := sample %in% truth_carriers]
eval_benchmark[, c("truth_vid", "truth_carriers") := list(NULL)]
grp_by_cols <- colnames(eval_benchmark)[colnames(eval_benchmark) != "is_true_de_novo"]
eval_benchmark <- eval_benchmark[, list(is_true_de_novo = any(is_true_de_novo)), by = grp_by_cols]
fwrite(eval_benchmark, "denovo_svs-benchmark.tsv.gz", sep = "\t", quote = FALSE)

truth_benchmark <- truth_vs_start_hits[truth_carriers, on = "vid"]
truth_benchmark <- start_carriers[truth_benchmark, on = c(start_vid = "match_vid")]
eval_carriers_by_vid <- denovo_calls[is_de_novo == TRUE, list(eval_samples = list(sample)), by = "name"]
truth_benchmark <- eval_carriers_by_vid[truth_benchmark, on = c(`name` = "start_vid")]
truth_benchmark[, `:=`(start_match = sample %in% start_samples, eval_match = sample %in% eval_samples)]
truth_benchmark <- truth_benchmark[, list(in_eval = any(eval_match), in_start = any(start_match)),
                                   by = c("chr", "start", "end", "svtype", "name", "sample")]
setnames(truth_benchmark, "name", "vid")
truth_benchmark <- truth_vcf_annotations[truth_benchmark, by = "vid"]
fwrite(truth_benchmark, "truth_denovos-benchmark.tsv.gz", sep = "\t", quote = FALSE)
