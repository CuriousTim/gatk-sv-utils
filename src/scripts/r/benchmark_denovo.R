# Usage: Rscript benchmark_denovo.R <denovo> <truth> <cleanvcf> <sample_table>

# if (!require("R.utils", quietly = TRUE)) {
#     install.packages("R.utils")
# }
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
test_path <- args[[1]]
truth_path <- args[[2]]
cleanvcf_path <- args[[3]]
sample_table_path <- args[[4]]

svmatch <- function(query, subject, min_ovp, size_similarity, breakend_window) {
    browser()
    if (length(query) == 0 || length(subject) == 0) {
        return(NULL)
    }

    hits <- findOverlaps(query, subject, maxgap = breakend_window)
    qr <- query[queryHits(hits)]
    sr <- subject[subjectHits(hits)]
    pi <- pintersect(qr, sr)
    qov <- width(pi) / width(qr)
    sov <- width(pi) / width(sr)

    ovp_pass <- qov >= min_ovp & sov >= min_ovp
    sample_pass <- mcols(qr)$other_id == mcols(sr)$other_id
    size_pass <- pmin(width(qr), width(sr)) / pmax(width(qr), width(sr)) > size_similarity

    qr[ovp_pass & sample_pass & size_pass]
}

# Which SVs in query have a matching SV in subject
overlap <- function(query, subject) {
    # small DEL
    sm_del <- svmatch(
        query[mcols(query)$svtype == "DEL" & width(query) < 500],
        subject[mcols(subject)$svtype == "DEL"],
        0, 0, 300
    )
    # large DEL
    lg_del <- svmatch(
        query[mcols(query)$svtype == "DEL" & width(query) >= 500],
        subject[mcols(subject)$svtype == "DEL"],
        0.1, 0.5, 5000
    )
    # small INS
    sm_ins <- svmatch(
      query[mcols(query)$svtype == "INS" & width(query) < 500],
      subject[mcols(subject)$svtype == "INS"],
      0, 0, 300
    )
    # large INS
    lg_ins <- svmatch(
      query[mcols(query)$svtype == "INS" & width(query) >= 500],
      subject[mcols(subject)$svtype == "INS"],
      0.1, 0.5, 5000
    )
    # DUP
    dup <- svmatch(
      query[mcols(query)$svtype == "DUP"],
      subject[mcols(subject)$svtype == "DUP"],
      0.1, 0.5, 5000
    )
    # INV
    inv <- svmatch(
      query[mcols(query)$svtype == "INV"],
      subject[mcols(subject)$svtype == "INV"],
      0.1, 0.5, 5000
    )

    results <- list(sm_del, lg_del, sm_ins, lg_ins, dup, inv)

    sort(unique(unlist(sapply(Filter(Negate(is.null), results), \(x) mcols(x)$src_row))))
}

dt2gr <- function(dt) {
    gr <- GRanges(dt$chr,
                  IRanges(start = dt$start, end = dt$end),
                  other_id = dt$other_id, svtype = dt$svtype,
                  src_row = seq_len(nrow(dt)))
    seqlevels(gr) <- paste0("chr", c(1:22, "X", "Y"))
    gr
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
truth <- truth[(!svtype %in% c("CPX", "CTX", "CNV")) & (mosaic == FALSE), ]
setnames(truth, "sample", "other_id")

test <- fread(test_path, sep = "\t")
setnames(test, "svtype", "ALT")
setnames(test, c("chrom", "SVTYPE"), c("chr", "svtype"))
# Exclude CTX and CPX from comparisons
test <- test[!svtype %in% c("CPX", "CTX"), ]
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
message(sprintf("keeping %d / %d calls from truth after matching samples in CleanVcf", sum(truth$other_id %in% cleanvcf_samples), nrow(truth)))
truth <- truth[other_id %in% cleanvcf_samples, ]

dn <- test[is_de_novo == TRUE, ]
not_dn <- test[is_de_novo == FALSE, ]

dn_gr <- dt2gr(dn)
not_dn_gr <- dt2gr(not_dn)
cleanvcf_gr <- dt2gr(cleanvcf)
truth_gr <- dt2gr(truth)

# Which truth calls are in CleanVcf -------------------------------------------
truth_in_cleanvcf <- overlap(truth_gr, cleanvcf_gr)
message(sprintf("%d / %d truth calls matched in CleanVcf", length(truth_in_cleanvcf), length(truth_gr)))

# Which truth calls are not in CleanVcf ---------------------------------------
truth_not_in_cleanvcf <- setdiff(seq_along(truth_gr), truth_in_cleanvcf)

truth_in_cleanvcf_dt <- truth[truth_in_cleanvcf, ]
truth_in_cleanvcf_gr <- dt2gr(truth_in_cleanvcf_dt)

# PPV (positive predictive value): which de novo calls are in truth -----------
# true positive / all positive
dn_in_truth <- overlap(dn_gr, truth_in_cleanvcf_gr)
message(sprintf("PPV: %d / %d = %0.3f", length(dn_in_truth), length(dn_gr), length(dn_in_truth) / length(dn_gr)))

# FDR (false discovery rate): which de novo calls are not in truth ------------
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

# missed calls: which truth calls are not in de novo calls --------------------
truth_not_in_dn <- setdiff(seq_along(truth_in_cleanvcf_gr), truth_in_dn)

fwrite(truth[truth_in_cleanvcf, ], "truth_in_cleanvcf.tsv", sep = "\t")
fwrite(truth[truth_not_in_cleanvcf, ], "truth_not_in_cleanvcf.tsv", sep = "\t")
fwrite(dn[dn_in_truth, ], "denovo_in_truth.tsv", sep = "\t")
fwrite(dn[dn_not_in_truth, ], "denovo_not_in_truth.tsv", sep = "\t")
fwrite(not_dn[not_dn_in_truth, ], "not_denovo_in_truth.tsv", sep = "\t")
fwrite(not_dn[not_dn_not_in_truth, ], "not_denovo_not_in_truth.tsv", sep = "\t")
fwrite(truth_in_cleanvcf_dt[truth_in_dn, ], "truth_in_denovo.tsv", sep = "\t")
fwrite(truth_in_cleanvcf_dt[truth_not_in_dn, ], "truth_not_in_denovo.tsv", sep = "\t")

old_par <- par(no.readonly = TRUE)
jpeg("benchmark_venn.jpeg", width = 800, height = 450)
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
symbols(0.35, 0.5, circles = 0.45, inches = FALSE, add = TRUE, bg = "#1ECBE166")
symbols(0.35, 0.6, circles = 0.35, inches = FALSE, add = TRUE, bg = "#E11ECB66")
symbols(0.7, 0.4, circles = 0.35, inches = FALSE, add = TRUE, bg = "#CBE11E66")
text(0.9, 0.4, labels = length(truth_not_in_cleanvcf), cex = 2, adj = 0.5)
text(0.55, 0.2, labels = length(truth_not_in_dn), cex = 2, adj = 0.5)
text(0.5, 0.5, labels = length(dn_in_truth), cex = 2, adj = 0.5)
text(0.3, 0.7, labels = length(dn_not_in_truth), cex = 2, adj = 0.5)
text(0.2, 0.2, labels = length(not_dn_not_in_truth), cex = 2, adj = 0.5)
dev.off()
