# Usage: Rscript benchmark_denovo.R <denovo> <truth> <cleanvcf> <sample_table>

args <- commandArgs(trailingOnly = TRUE)
test_path <- args[[1]]
truth_path <- args[[2]]
cleanvcf_path <- args[[3]]
sample_table_path <- args[[4]]

# if (!require("R.utils", quietly = TRUE)) {
#     install.packages("R.utils")
# }
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

# Which SVs in query have a matching SV in subject
overlap <- function(query, subject, min_ovp = 0.5) {
    hits <- findOverlaps(query, subject)
    qr <- query[queryHits(hits)]
    sr <- subject[subjectHits(hits)]
    pi <- pintersect(qr, sr)
    qov <- width(pi) / width(qr)
    sov <- width(pi) / width(sr)

    ovp_pass <- qov >= min_ovp & sov >= min_ovp
    sample_pass <- mcols(qr)$sample_id == mcols(sr)$sample_id
    svtype_pass <- mcols(qr)$svtype == mcols(sr)$svtype

    sort(unique(queryHits(hits)[ovp_pass & sample_pass & svtype_pass]))
}

# Sample ID conversion --------------------------------------------------------

sample_table <- fread(sample_table_path,
                      sep = "\t",
                      select = c(
                          "entity:sample_id",
                          "other_id",
                          "cohort_short")
)
setnames(sample_table, "entity:sample_id", "sample_id")

truth <- fread(truth_path, sep = ",")
setnames(truth, c("#chrom", "pos", "end"), c("chr", "start", "end"))
truth[, svtype := sub("LINE1|SVA|ALU", "INS", svtype)]

test <- fread(test_path, sep = "\t")
setnames(test, "chrom", "chr")
test <- merge(test, sample_table, by.x = "sample", by.y = "sample_id", all.x = TRUE)

# bcftools query -i 'GT="alt" & INFO/SVTYPE != "CNV" & INFO/SVTYPE != "BND"' \
#   -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t[%SAMPLE,]\n' ssc_quads-annotated.vcf.gz \
#   | awk -F'\t' '{sub(/,$/, "", $5); split($5, a, /,/); for(i in a){$5=a[i]; print}}' OFS='\t' > clean_vcf_svs.tsv
cleanvcf <- fread(cleanvcf_path,
                  sep = "\t",
                  col.names = c("chr", "start", "end", "svtype", "sample"),
                  header = FALSE)
cleanvcf <- merge(cleanvcf, sample_table, by.x = "sample", by.y = "sample_id", all.x = TRUE)

# Remove non-SSC from CleanVcf
cleanvcf <- cleanvcf[cohort_short == "SSC", ]

cleanvcf_samples <- unique(cleanvcf$other_id)

# Subset the test set and truth set to SSC samples in CleanVcf
test <- test[other_id %in% cleanvcf_samples, ]
message(sprintf("keeping %d / %d calls from truth after matching samples in CleanVcf", sum(truth$sample %in% cleanvcf_samples), nrow(truth)))
truth <- truth[sample %in% cleanvcf_samples, ]

# Exclude CTX and CPX from comparisons and CNV from truth
test <- test[!SVTYPE %in% c("CPX", "CTX"), ]
truth <- truth[!svtype %in% c("CPX", "CTX", "CNV"), ]
dn <- test[is_de_novo == TRUE, ]
not_dn <- test[is_de_novo == FALSE, ]

chr <- paste0("chr", c(1:22, "X", "Y"))
test_gr <- GRanges(test$chr, IRanges(start = test$start, end = test$end),
                   sample_id = test$other_id, svtype = test$SVTYPE, is_de_novo = test$is_de_novo)
seqlevels(test_gr) <- chr
dn_gr <- GRanges(dn$chr, IRanges(start = dn$start, end = dn$end),
                 sample_id = dn$other_id, svtype = dn$SVTYPE, is_de_novo = dn$is_de_novo)
seqlevels(dn_gr) <- chr
not_dn_gr <- GRanges(not_dn$chr, IRanges(start = not_dn$start, end = not_dn$end),
                 sample_id = not_dn$other_id, svtype = not_dn$SVTYPE, is_de_novo = not_dn$is_de_novo)
seqlevels(not_dn_gr) <- chr
cleanvcf_gr <- GRanges(cleanvcf$chr, IRanges(start = cleanvcf$start, end = cleanvcf$end),
                       sample_id = cleanvcf$other_id, svtype = cleanvcf$svtype)
seqlevels(cleanvcf_gr) <- chr
truth_gr <- GRanges(truth$chr, IRanges(start = truth$start, end = truth$end),
                       sample_id = truth$sample, svtype = truth$svtype)
seqlevels(truth_gr) <- chr

# Which truth calls are in CleanVcf -------------------------------------------
truth_in_cleanvcf <- overlap(truth_gr, cleanvcf_gr)
message(sprintf("%d / %d truth calls matched in CleanVcf", length(truth_in_cleanvcf), length(truth_gr)))

# Which truth calls are not in CleanVcf ---------------------------------------
truth_not_in_cleanvcf <- setdiff(seq_along(truth_gr), truth_in_cleanvcf)

truth_in_cleanvcf_gr <- truth_gr[truth_in_cleanvcf]

# PPV (positive predictive value): which de novo calls are in truth -----------
# true positive / all positive
dn_in_truth <- overlap(dn_gr, truth_in_cleanvcf_gr)
message(sprintf("PPV: %d / %d = %0.3f", length(dn_in_truth), length(dn_gr), length(dn_in_truth) / length(dn_gr)))

# FDR (false discovery rate) which de novo calls are not in truth -------------
# false positive / all positive
dn_not_in_truth <- setdiff(seq_along(dn_gr), dn_in_truth)
message(sprintf("FDR: %d / %d = %0.3f", length(dn_not_in_truth), length(dn_gr), length(dn_not_in_truth) / length(dn_gr)))

# FOR (false omission rate): which not de novo calls are in truth -------------
# false negative / all negative
not_dn_in_truth <- overlap(not_dn_gr, truth_in_cleanvcf_gr)
message(sprintf("FOR: %d / %d = %0.3f", length(not_dn_in_truth), length(not_dn_gr), length(not_dn_in_truth) / length(not_dn_gr)))

# NPV (negative predictive value): which not de novo calls are not in truth ---
# true negative / all negative
not_dn_not_in_truth <- setdiff(seq_along(not_dn_gr), not_dn_in_truth)
message(sprintf("NPV: %d / %d = %0.3f", length(not_dn_not_in_truth), length(not_dn_gr), length(not_dn_not_in_truth) / length(not_dn_gr)))

# sensitivity: which truth calls are in de novo calls -------------------------
# captured positive / true positives
truth_in_dn <- overlap(truth_in_cleanvcf_gr, dn_gr)
message(sprintf("sensitivity: %d / %d = %0.3f", length(truth_in_dn), length(truth_in_cleanvcf_gr), length(truth_in_dn) / length(truth_in_cleanvcf_gr)))

# missed calls: which truth calls are not in de novo calls -----------------------------
truth_not_in_dn <- setdiff(seq_along(truth_in_cleanvcf_gr), truth_in_dn)

fwrite(truth[truth_in_cleanvcf, ], "truth_in_cleanvcf.tsv", sep = "\t")
fwrite(truth[truth_not_in_cleanvcf, ], "truth_not_in_cleanvcf.tsv", sep = "\t")
fwrite(dn[dn_in_truth, ], "denovo_in_truth.tsv", sep = "\t")
fwrite(dn[dn_not_in_truth, ], "denovo_not_in_truth.tsv", sep = "\t")
fwrite(not_dn[not_dn_in_truth, ], "not_denovo_in_truth.tsv", sep = "\t")
fwrite(not_dn[not_dn_not_in_truth, ], "not_denovo_not_in_truth.tsv", sep = "\t")
fwrite(truth[truth_in_dn, ], "truth_in_denovo.tsv", sep = "\t")
fwrite(truth[truth_not_in_dn, ], "truth_not_in_denovo.tsv", sep = "\t")
