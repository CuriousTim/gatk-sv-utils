plot_venn <- function(eval_vs_truth, truth_vs_eval, gc = c("GN", "SR", "RM", "SD", "UN", "UN_or_GN"), path = NULL) {
    gc <- match.arg(gc)
    if (!is.null(path)) {
        jpeg(path, width = 800, height = 450)
        on.exit(dev.off(), add = TRUE)
    }
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    tp <- sum(eval_vs_truth$tp)
    fp <- sum(eval_vs_truth$fp)
    fn1 <- sum(truth_vs_eval$fn1)
    fn2 <- sum(truth_vs_eval$fn2)
    sens <- if (tp + fn1 == 0) NA_real_ else tp / (tp + fn1)
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
    text(1.2, 0.8, labels = sprintf("Sensitivity:\n%0.2f", sens), cex = 1.5, adj = 0.5)
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

evt <- read.table(
    "eval_bench.tsv.gz",
    sep = "\t",
    col.names = c("chr", "start", "end", "svtype", "vid", "context", "ovp_gene", "tp", "fp"),
    colClasses = c(
        "character",
        "integer",
        "integer",
        "character",
        "character",
        "character",
        "integer",
        "integer",
        "integer"
    )
)
# convert to 1-start
evt$start <- evt$start + 1L
tve <- read.table(
    "truth_bench.tsv.gz",
    sep = "\t",
    col.names = c("chr", "start", "end", "svtype", "vid", "context", "ovp_gene", "fn1", "fn2"),
    colClasses = c(
        "character",
        "integer",
        "integer",
        "character",
        "character",
        "character",
        "integer",
        "integer",
        "integer"
    )
)
tve$start <- tve$start + 1L

plot_venn(
    evt[evt$context == "SD", ],
    tve[tve$context == "SD", ],
    "SD",
    if (interactive()) NULL else "denovo_benchmark_SD.jpg"
)
plot_venn(
    evt[evt$context == "RM", ],
    tve[tve$context == "RM", ],
    "RM",
    if (interactive()) NULL else "denovo_benchmark_RM.jpg"
)
plot_venn(
    evt[evt$context == "SR", ],
    tve[tve$context == "SR", ],
    "SR",
    if (interactive()) NULL else "denovo_benchmark_SR.jpg"
)
plot_venn(
    evt[evt$context == "UN", ],
    tve[tve$context == "UN", ],
    "UN",
    if (interactive()) NULL else "denovo_benchmark_UN.jpg"
)
plot_venn(
    evt[evt$ovp_gene == TRUE, ],
    tve[tve$ovp_gene == TRUE, ],
    "GN",
    if (interactive()) NULL else "denovo_benchmark_PCG.jpg"
)
plot_venn(
    evt[evt$ovp_gene == TRUE | evt$context == "UN", ],
    tve[tve$ovp_gene == TRUE | tve$context == "UN", ],
    "UN_or_GN",
    if (interactive()) NULL else "denovo_benchmark_UN_or_PCG.jpg"
)

conn <- gzfile("eval_bench-with_header.tsv.gz", open = "w")
write.table(evt, file = conn, sep = "\t", row.names = FALSE, quote = FALSE)
close(conn)
conn <- gzfile("truth_bench-with_header.tsv.gz", open = "w")
write.table(tve, file = conn, sep = "\t", row.names = FALSE, quote = FALSE)
close(conn)
