library(data.table)

plot_venn <- function(eval_vs_truth, truth_vs_eval, gc = c("GN", "SR", "RM", "SD", "UN", "UN_or_GN"), path = NULL) {
    gc <- match.arg(gc)
    if (!is.null(path)) {
        jpeg(path, width = 800, height = 450)
        on.exit(dev.off(), add = TRUE)
    }
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    tp <- sum(eval_vs_truth$is_de_novo & eval_vs_truth$is_true_de_novo)
    fp <- sum(eval_vs_truth$is_de_novo & (!eval_vs_truth$is_true_de_novo))
    fn1 <- sum((!truth_vs_eval$in_eval) & truth_vs_eval$in_start)
    fn2 <- sum(!truth_vs_eval$in_start)
    sens1 <- if (tp + fn1 == 0) NA_real_ else tp / (tp + fn1)
    sens2 <- if (tp + fn1 + fn2 == 0) NA_real_ else tp / (tp + fn1 + fn2)
    prec <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)

    par(mar = c(0, 0, 0, 0))
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
    symbols(0.3, 0.5, circles = 0.45, inches = FALSE, add = TRUE, bg = "#94E6F166")
    symbols(0.7, 0.5, circles = 0.45, inches = FALSE, add = TRUE, bg = "#B8F19466")
    symbols(0.3, 0.65, circles = 0.30, inches = FALSE, add = TRUE, bg = "#CD94F166")
    text(-0.15, 0.9, labels = expand_gc(gc), cex = 1.5, adj = 0.5)
    text(0.01, 0.3, labels = "input VCF", cex = 1, adj = 0.5)
    text(0.9, 0.3, labels = "Belyeu", cex = 1, adj = 0.5)
    text(0.3, 0.9, labels = "predicted de novo", cex = 1, adj = 0.5)
    text(0.9, 0.5, labels = fn2, cex = 2, adj = 0.5)
    text(0.5, 0.3, labels = fn1, cex = 2, adj = 0.5)
    text(0.4, 0.6, labels = tp, cex = 2, adj = 0.5)
    text(0.15, 0.7, labels = fp, cex = 2, adj = 0.5)
    text(1.2, 0.8, labels = sprintf("Sensitivity 1:\n%0.2f", sens1), cex = 1.5, adj = 0.5)
    text(-0.2, 0.8, labels = sprintf("Sensitivity 2:\n%0.2f", sens2), cex = 1.5, adj = 0.5)
    text(1.2, 0.2, labels = sprintf("Precision:\n%0.2f", prec), cex = 1.5, adj = 0.5)
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
    } else if (x == "UN") {
        return("unique regions")
    } else if (x == "UN_or_GN") {
        return("unique or gene")
    }
}

evt <- readLines("eval_benchmarks.txt") |>
    lapply(\(x) fread(x, sep = "\t", header = TRUE)) |>
    rbindlist()

tve <- readLines("truth_benchmarks.txt") |>
    lapply(\(x) fread(x, sep = "\t", header = TRUE)) |>
    rbindlist()

plot_venn(
    evt[evt$genomic_context == "SD", ],
    tve[tve$genomic_context == "SD", ],
    "SD",
    if (interactive()) NULL else "denovo_benchmark_SD.jpg"
)
plot_venn(
    evt[evt$genomic_context == "RM", ],
    tve[tve$genomic_context == "RM", ],
    "RM",
    if (interactive()) NULL else "denovo_benchmark_RM.jpg"
)
plot_venn(
    evt[evt$genomic_context == "SR", ],
    tve[tve$genomic_context == "SR", ],
    "SR",
    if (interactive()) NULL else "denovo_benchmark_SR.jpg"
)
plot_venn(
    evt[evt$genomic_context == "UN", ],
    tve[tve$genomic_context == "UN", ],
    "UN",
    if (interactive()) NULL else "denovo_benchmark_UN.jpg"
)
plot_venn(
    evt[evt$ovp_pc_gene == TRUE, ],
    tve[tve$ovp_pc_gene == TRUE, ],
    "GN",
    if (interactive()) NULL else "denovo_benchmark_PCG.jpg"
)
plot_venn(
    evt[evt$ovp_pc_gene == TRUE | evt$genomic_context == "UN", ],
    tve[tve$ovp_pc_gene == TRUE | tve$genomic_context == "UN", ],
    "UN_or_GN",
    if (interactive()) NULL else "denovo_benchmark_UN_or_PCG.jpg"
)

fwrite(evt, "eval_benchmark.tsv.gz", sep = "\t", quote = FALSE)
fwrite(tve, "truth_benchmark.tsv.gz", sep = "\t", quote = FALSE)
