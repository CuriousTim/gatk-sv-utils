library(rtracklayer)
library(GenomicRanges)

plot_venn <- function(eval_vs_truth, truth_vs_eval, gc = c("GN", "SR", "RM", "SD", "UN"), path = NULL) {
  gc <- match.arg(gc)
  tmp1 <- eval_vs_truth[mcols(eval_vs_truth)$gc == gc]
  tmp2 <- truth_vs_eval[mcols(truth_vs_eval)$gc == gc]
  if (!is.null(path)) {
    jpeg(path, width = 800, height = 450)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(0, 0, 0, 0))
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
  symbols(0.3, 0.5, circles = 0.45, inches = FALSE, add = TRUE, bg = "#94E6F166")
  symbols(0.7, 0.5, circles = 0.45, inches = FALSE, add = TRUE, bg = "#B8F19466")
  symbols(0.3, 0.65, circles = 0.30, inches = FALSE, add = TRUE, bg = "#CD94F166")
  text(-0.15, 0.9, labels = expand_gc(gc), cex = 1.5, adj = 0.5)
  text(0.01, 0.3, labels = "input VCF", cex = 1, adj = 0.5)
  text(0.9, 0.3, labels = "Belyeu", cex = 1, adj = 0.5)
  text(0.3, 0.9, labels = "predicted de novo", cex = 1, adj = 0.5)
  text(0.9, 0.5, labels = sum(mcols(tmp2)$fn2), cex = 2, adj = 0.5)
  text(0.5, 0.3, labels = sum(mcols(tmp2)$fn1), cex = 2, adj = 0.5)
  text(0.4, 0.6, labels = sum(mcols(tmp1)$tp), cex = 2, adj = 0.5)
  text(0.15, 0.7, labels = sum(mcols(tmp1)$fp), cex = 2, adj = 0.5)
}

expand_gc <- function(x) {
  if (x == "GN") {
    return("protein-coding gene")
  } else if (x == "SR") {
    return("simple repeats")
  } else if (x == "RM") {
    return("RepeatMasker")
  } else if (x == "SD") {
    return("segmental duplications")
  } else {
    return("unique regions")
  }
}

coverage <- function(query, subject) {
  ovps <- findOverlaps(query, subject, ignore.strand = TRUE)
  qh <- queryHits(ovps)
  sh <- subjectHits(ovps)
  cov <- double(length(query))
  int <- pintersect(query[qh], subject[sh])
  qhf <- factor(qh)
  tmp <- by(width(int) / width(query[qh]), qhf, sum)
  cov[as.integer(names(tmp))] <- tmp
  cov
}

evt <- read.table("eval_vs_truth.tsv", sep = "\t", col.names = c("chr", "start", "end", "vid", "tp", "fp"))
tve <- read.table("truth_vs_eval.tsv", sep = "\t", col.names = c("chr", "start", "end", "vid", "fn1", "fn2"))

evt_gr <- as(evt, "GRanges")
tve_gr <- as(tve, "GRanges")

sr_gr <- import("hg38_SR.bed.gz")
rm_gr <- import("hg38_RM.bed.gz")
sd_gr <- import("hg38_SD.bed.gz")
gencode <- import("gencode.v49.basic.annotation.gff3.gz")
pc_genes <- gencode[mcols(gencode)$type == "gene" & mcols(gencode)$gene_type == "protein_coding"]

mcols(evt_gr)$gc <- "UN"
mcols(evt_gr)$gc[overlapsAny(evt_gr, sr_gr, ignore.strand = TRUE)] <- "SR"
mcols(evt_gr)$gc[overlapsAny(evt_gr, rm_gr, ignore.strand = TRUE)] <- "RM"
mcols(evt_gr)$gc[overlapsAny(evt_gr, sd_gr, ignore.strand = TRUE)] <- "SD"
mcols(evt_gr)$gc[overlapsAny(evt_gr, pc_genes, ignore.strand = TRUE)] <- "GN"

mcols(tve_gr)$gc <- "UN"
mcols(tve_gr)$gc[overlapsAny(tve_gr, sr_gr, ignore.strand = TRUE)] <- "SR"
mcols(tve_gr)$gc[overlapsAny(tve_gr, rm_gr, ignore.strand = TRUE)] <- "RM"
mcols(tve_gr)$gc[overlapsAny(tve_gr, sd_gr, ignore.strand = TRUE)] <- "SD"
mcols(tve_gr)$gc[overlapsAny(tve_gr, pc_genes, ignore.strand = TRUE)] <- "GN"

plot_venn(evt_gr, tve_gr, "GN", if (interactive()) NULL else "denovo_benchmark_GN.jpg")
plot_venn(evt_gr, tve_gr, "SD", if (interactive()) NULL else "denovo_benchmark_SD.jpg")
plot_venn(evt_gr, tve_gr, "RM", if (interactive()) NULL else "denovo_benchmark_RM.jpg")
plot_venn(evt_gr, tve_gr, "SR", if (interactive()) NULL else "denovo_benchmark_SR.jpg")
plot_venn(evt_gr, tve_gr, "UN", if (interactive()) NULL else "denovo_benchmark_UN.jpg")

mcols(evt_gr)$gc <- "UN"
mcols(evt_gr)$gc[coverage(evt_gr, sr_gr) > 0.5] <- "SR"
mcols(evt_gr)$gc[coverage(evt_gr, rm_gr) > 0.5] <- "RM"
mcols(evt_gr)$gc[coverage(evt_gr, sd_gr) > 0.5] <- "SD"

mcols(tve_gr)$gc <- "UN"
mcols(tve_gr)$gc[coverage(tve_gr, sr_gr) > 0.5] <- "SR"
mcols(tve_gr)$gc[coverage(tve_gr, rm_gr) > 0.5] <- "RM"
mcols(tve_gr)$gc[coverage(tve_gr, sd_gr) > 0.5] <- "SD"

plot_venn(evt_gr, tve_gr, "SD", if (interactive()) NULL else "denovo_benchmark_SD_cov50.jpg")
plot_venn(evt_gr, tve_gr, "RM", if (interactive()) NULL else "denovo_benchmark_RM_cov50.jpg")
plot_venn(evt_gr, tve_gr, "SR", if (interactive()) NULL else "denovo_benchmark_SR_cov50.jpg")
plot_venn(evt_gr, tve_gr, "UN", if (interactive()) NULL else "denovo_benchmark_UN_cov50.jpg")
