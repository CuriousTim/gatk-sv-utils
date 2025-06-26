#!/usr/bin/env Rscript

library(optparse)

# Options parsing ------------------------------------------------------------
parser <- OptionParser()

# See `read_cnvs()`
parser <- add_option(parser, c("-b", "--cnvs"), type = "character",
                     metavar = "<path>",
                     help = "File with CNVs to plot")
# File should be two tab-separated columns without a header
# 1. sample ID
# 2. batch ID
parser <- add_option(parser, c("-s", "--sample-batches"), type = "character",
                     metavar = "<path>",
                     help = "Mapping between batch ID and sample ID for every sample in the cohort")
# File should be two tab-separated columns without a header
# 1. batch ID
# 2. path to coverage matrix
parser <- add_option(parser, c("-d", "--coverage-paths"), type = "character",
                     metavar = "<path>",
                     help = "Mapping between batch ID and coverage matrix")
# File should be two tab-separated columns without a header
# 1. batch ID
# 2. path to coverage medians
parser <- add_option(parser, c("-m", "--medians-paths"), type = "character",
                     metavar = "<path>",
                     help = "Mapping between batch ID and coverage medians")
parser <- add_option(parser, c("-o", "--output"), type = "character",
                     metavar = "<path>",
                     help = "Path to output directory")

opts <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
if (is.null(opts$cnvs)) {
    stop("-b,--cnvs is required")
}
if (is.null(opts$sample_batches)) {
    stop("-s,--sample-batches is required")
}
if (is.null(opts$coverage_paths)) {
    stop("-d,--coverage-paths is required")
}
if (is.null(opts$medians_paths)) {
    stop("-m,--medians-paths is required")
}
if (is.null(opts$output)) {
    stop("-o,--output is required")
}

# Constants -------------------------------------------------------------------

# The intervals in the bincov matrices are 100 bp so the matrices large CNVs
# take up a lot of memory. For these CNVs, we sample the median coverage at
# equally spaced windows across the SV and plot the samples.

# number of samples to take
LARGE_CNV_SAMPLE_COUNT <- 500L
# size of the sampling window
LARGE_CNV_SAMPLE_WINDOW_SIZE <- 2000L
# minimum CNV size to use sampling strategy
LARGE_CNV_SIZE <- LARGE_CNV_SAMPLE_COUNT * LARGE_CNV_SAMPLE_WINDOW_SIZE

# Seeing the coverage upstream and downstream of the CNV is useful for
# determining whether the breakpoints are accurate so we add some padding
# to the CNV when plotting.

# Fraction of CNV length to add to the CNV as padding
PAD_EXPANSION_FACTOR <- 0.5

# must be odd
SMOOTH_WINDOW <- 21

# number of intervals to plot
INTERVAL_PLOT_COUNT <- 26L

MAX_BACKGROUND_SAMPLES <- 500L

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

# Functions ------------------------------------------------------------------

#' Read the CNVs to plot.
#'
#' File should be six tab-separated columns without a header
#' 1. chromosome
#' 2. start (1-based, inclusive)
#' 3. end (1-based, inclusive)
#' 4. variant ID
#' 5. SV type
#' 6. sample ID
#'
#' @param path `character(1)` Path to the file.
#' @returns `data.table` The table in `path`.
read_cnvs <- function(path) {
    cnvs <- fread(path, sep = "\t", header = FALSE,
                  col.names = c("chr", "start", "end", "vid", "svtype", "samples"),
                  colClasses = c("character", "integer", "integer", "character",
                                 "character", "character"))
    if (nrow(cnvs) == 0) {
        stop("no variants to plot")
    }

    if (!all(cnvs$end >= cnvs$start)) {
        stop("all CNV ends must be greater than or equal to starts")
    }

    if (!all(grepl("DEL|DUP", cnvs$svtype))) {
        stop("only DELs and DUPs can be plotted")
    }

    if (!all(nzchar(cnvs$vid))) {
        stop("all variant IDs must be non-empty")
    }

    if (any(grepl(.Platform$file.sep, cnvs$vid, fixed = TRUE))) {
        stop("variant IDs must not contain path separators")
    }

    cnvs
}

#' Read a two-column, tab-separated file into a hash table.
#'
#' The first column contains the keys and the second column contains the
#' values.
#'
#' @param path `character(1)` Path to the file.
#' @returns `utils::hashtab` Hash table.
read_keyed_tsv <- function(path) {
    d <- fread(path, sep = "\t", header = FALSE,
                    col.names = c("V1", "V2"),
                    colClasses = c("character", "character"))
    if (nrow(d) == 0) {
        stop("at least one row must be given in mapping")
    }

    if (anyDuplicated(d$V1) != 0) {
        stop("TSV keys must be unique")
    }

    h <- hashtab(type = "identical", nrow(d))
    mapply(\(k, v) sethash(h, k, v), d$V1, d$V2)

    h
}

#' Retrieve the header of a bincov matrix.
#'
#' @param con `TabixFile` Connection to matrix.
#' @returns `character` Header of the matrix.
bincov_header <- function(con) {
    header <- headerTabix(con)[["header"]]
    if (is.null(header)) {
        stop(sprintf("coverage matrix at '%s' is missing a header", con$ath))
    }
    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) <= 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop(sprintf("coverage matrix at '%s' has an invalid header", con$path))
    }
    header[1:3] <- c("chr", "start", "end")

    header
}

#' Query a bincov matrix.
#'
#' @param con `TabixFile` Connection to bincov matrix.
#' @param header `character` The parsed header line from the bincov matrix.
#'   Assumes the first three elements are "chr", "start", "end" and the
#'   remaining elements are sample IDs.
#' @param samples `character` Samples to keep.
#' @param ranges `GRanges` Genomic ranges to check for overlap.
#' @returns `data.table` The bincov matrix. If `ranges` contains a single
#'   range, the overlapping bincov rows will be returned as is. Otherwise,
#'   the median coverage of each sample for the bins overlapping each
#'   range will be returned. Ranges without any overlapping bins will have
#'   `NA` coverage values unless all ranges do not have any overlapping bins
#'   in which case an error is signaled.
query_bincov <- function(con, header, samples, ranges) {
    tabix <- scanTabix(con, param = ranges)
    if (all(lengths(tabix) == 0)) {
        stop("no overlapping intervals in coverage matrix")
    }

    qcoords <- data.table(
        qchr = as.character(seqnames(ranges)),
        qstart = start(ranges),
        qend = end(ranges)
    )

    coltypes <- list("character" = 1L,
                     "integer" = c(2L, 3L),
                     "double" = seq.int(4, length(header)))
    cols_to_keep <- which(header %in% samples)
    samples_to_keep <- header[cols_to_keep]
    parse <- function(qchr, qstart, qend, lines) {
        if (length(lines) == 0) {
            return(NULL)
        }
        d <- fread(text = lines, sep = "\t", header = FALSE,
                   colClasses = coltypes, select = c(1:3, cols_to_keep))
        colnames(d) <- c("chr", "start", "end", samples_to_keep)
        # convert from 0-start, exclusive to 1-start, inclusive
        d[, start := start + 1L]
        d[, c("qchr", "qstart", "qend") := list(qchr, qstart, qend)]
        d
    }

    mat <- mapply(parse, qcoords$qchr, qcoords$qstart, qcoords$qend, tabix, SIMPLIFY = FALSE, USE.NAMES = FALSE) |>
        rbindlist(use.names = TRUE)

    # a single range indicates the query region is the span of the CNV while multiple regions
    # indicates interval sampling
    if (length(ranges) > 1) {
        mat <- mat[qcoords, on = c("qchr", "qstart", "qend")]
        # each query interval result should be an NA row, or not have any NA rows so computing
        # median without na.rm is ok
        mat <- mat[, lapply(.SD, median), .SDcols = samples_to_keep, by = c("qchr", "qstart", "qend")]
        setnames(mat, colnames(qcoords), c("chr", "start", "end"))
    } else {
        mat[, c("qchr", "qstart", "qend") := list(NULL)]
    }
    setkey(mat, chr, start, end)

    mat
}

#' Read coverage medians files into a named vector.
#'
#' @param paths `character` Paths to the coverage medians files.
#' @returns `double` The coverage medians. The names of the vector are the
#'   sample IDs.
read_medians <- function(paths) {
    input <- lapply(paths, \(x) readLines(x, n = 2L, ok = FALSE))
    ids <- unlist(strsplit(unlist(lapply(input, \(x) x[[1]])), split = "\t", fixed = TRUE))
    medians <- strsplit(unlist(lapply(input, \(x) x[[2]])), split = "\t", fixed = TRUE) |>
        unlist() |>
        as.double()

    names(medians) <- ids

    medians
}

#' Normalize coverage values.
#'
#' Normalization is done by dividing each sample's coverage values by the
#' genome-wide median coverage of the sample.
#'
#' @param x `data.table` Coverage matrix.
#' @param medians `double` Genome-wide median coverage for each sample in `x`.
#' @param `data.table` Normalized coverage matrix.
normalize_cov <- function(x, medians) {
    cov_cols <- colnames(x)
    cov_samples <- cov_cols[!cov_cols %in% c("chr", "start", "end")]
    medians_samples <- names(medians)
    common_samples <- intersect(cov_samples, medians_samples)
    if (length(cov_samples) != length(common_samples)) {
        stop("some samples in coverage matrix do not have a median coverage")
    }

    medians <- medians[common_samples]
    x[, names(.SD) := mapply(`/`, .SD, medians, SIMPLIFY = FALSE, USE.NAMES = FALSE), .SDcols = names(medians)]

    x
}

#' Tile the region of a CNV with equally spaced windows.
#'
#' The first window always begins at the start of the CNV and the last window
#' always ends at the end of the CNV. The intervening windows will be equally
#' spaced with excess gaps being distributed from back to front. Depends on
#' the LARGE_CNV_* constants being what they are.
#'
#' @param cnv `list` Information of the CNV to tile.
#' @returns `GRanges` Ranges of the windows.
tile_cnv <- function(cnv) {
    gaps <- LARGE_CNV_SAMPLE_COUNT - 1L
    total_gap_size <- cnv$end - cnv$start + 1L - LARGE_CNV_SIZE
    pad <- floor(total_gap_size / gaps)
    remaining_pad <- total_gap_size - pad * gaps
    pads <- rep(pad, gaps)
    i <- (gaps - remaining_pad + 1):gaps
    pads[i] <- pads[i] + 1L
    steps <- c(cnv$start, pads + LARGE_CNV_SAMPLE_WINDOW_SIZE)
    starts <- cumsum(steps)
    ends <- starts + LARGE_CNV_SAMPLE_WINDOW_SIZE - 1L

    GRanges(cnv$chr, IRanges(starts, ends))
}

# Like `tile_cnv()`, but different.
select_spaced_intervals <- function(n) {
    if (n <= INTERVAL_PLOT_COUNT) {
        return(1:n)
    }
    gaps <- INTERVAL_PLOT_COUNT - 1L
    total_gap_size <- n - INTERVAL_PLOT_COUNT
    pad <- floor(total_gap_size / gaps)
    remaining_pad <- total_gap_size - pad * gaps
    pads <- rep(pad, gaps)
    i <- (gaps - remaining_pad + 1):gaps
    pads[i] <- pads[i] + 1L
    steps <- c(1, pads + 1L)

    cumsum(steps)
}

rollmean <- function(x) {
    n <- length(x)
    if (n < SMOOTH_WINDOW) {
        return(frollmean(x, 1:n, adaptive = TRUE))
    }
    partial <- ceiling(SMOOTH_WINDOW / 2)
    a <- vapply(seq(partial, SMOOTH_WINDOW - 1),
                \(y) mean(x[1:y], na.rm = TRUE),
                double(1))
    b <- vapply(seq(SMOOTH_WINDOW - 1, partial),
                \(y) mean(x[(n - y + 1):n], na.rm = TRUE),
                double(1))
    # using algo = "fast" makes NA propagate even with na.rm = TRUE
    m <- frollmean(x, SMOOTH_WINDOW, algo = "exact", align = "center", na.rm = TRUE)
    m[1:(partial - 1)] <- a
    m[(n - partial + 2):n] <- b

    m
}

#' Load the coverage values over a CNV region.
#'
#' @param cnv `list` Information of the CNV.
#' @param paths `character` Paths to the coverage matrices.
#' @returns `data.table` Coverage matrix for the CNV region.
load_cnv_coverage <- function(cnv, paths) {
    if (cnv$end - cnv$start + 1L >= LARGE_CNV_SIZE) {
        ranges <- tile_cnv(cnv)
    } else {
        ranges <- GRanges(cnv$chr, IRanges(cnv$start, cnv$end))
    }

    cons <- TabixFileList(paths)
    headers <- lapply(cons, bincov_header)
    all_samples <- lapply(headers, \(x) x[-(1:3)]) |> unlist() |> unique()
    bg_samples <- all_samples[!all_samples %in% cnv$samples_split]
    bg_samples <- sample(bg_samples, min(MAX_BACKGROUND_SAMPLES, length(bg_samples)))
    keep_samples <- c(cnv$samples_split, bg_samples)

    mats <- mapply(query_bincov, cons, headers,
                   MoreArgs = list(samples = keep_samples, ranges = ranges),
                   SIMPLIFY = FALSE,
                   USE.NAMES = FALSE)

    # All of the coverage matrices should have the same intervals,
    # we merge just to be safe. This requires that `setkey()` be
    # called for every data.table.
    mats <- Reduce(merge, mats)

    if (nrow(mats) == 1) {
        stop("CNV must overlap at least two coverage intervals")
    }

    if (!all(cnv$samples_split %in% colnames(mats))) {
        stop("carrier samples are missing from the coverage matrix")
    }

    mats
}

#' Prettify the CNV size.
#'
#' @param cnv `list` Information of the CNV.
#' @returns `character(1)` A CNV size formatted as a pretty string.
pretty_cnv_size <- function(cnv) {
    cnv_size <- cnv$cnv_end - cnv$cnv_start + 1L
    if (cnv_size <= 1000) {
        size_pretty <- sprintf("%d b", cnv_size)
    } else if (cnv_size <= 100000) {
        size_pretty <- sprintf("%.2f Kb", cnv_size / 1000)
    } else {
        size_pretty <- sprintf("%.2f Mb", cnv_size / 1000000)
    }

    size_pretty
}

plot_singleton_intervals <- function(x, mids, col, cex) {
    idx <- which(!is.na(x))
    na_idx <- which(is.na(x))
    na_after <- idx[(idx + 1) %in% na_idx]
    na_before <- idx[(idx - 1) %in% na_idx]
    singleton_idx <- sort(intersect(na_after, na_before))
    points(mids[singleton_idx], x[singleton_idx], cex = cex, pch = 19, col = col)

    # Handle intervals at the start and end of vector
    if (length(x) < 2) {
        return()
    }
    if (!is.na(x[[1]]) && is.na(x[[2]])) {
        points(mids[[1]], x[[1]], cex = cex, pch = 19, col = col)
    }
    if (!is.na(x[[length(x)]]) && is.na(x[[length(x) - 1]])) {
        points(mids[[length(mids)]], x[[length(x)]], cex = cex, pch = 19, col = col)
    }
}

#' Plot coverage over a CNV region.
#'
#' @param cnv `list` Information of the CNV.
#' @param norm_cov `data.table` Normalized coverage matrix for the CNV region. Should be smoothed
#'   and sampled.
#' @param outfile The path to the plot.
make_plot <- function(cnv, norm_cov, outfile) {
    # make plot main
    size_pretty <- pretty_cnv_size(cnv)
    start_pretty <- formatC(cnv$cnv_start, big.mark = ",", format = "d")
    end_pretty <- formatC(cnv$cnv_end, big.mark = ",", format = "d")
    if (nchar(cnv$samples) > 30) {
        samples_pretty <- paste0(substr(cnv$samples, 1, 30), "...")
    } else {
        samples_pretty <- cnv$samples
    }
    main <- sprintf("%s:%s-%s (hg38)\n%s (%s)",
                    cnv$chr, start_pretty, end_pretty, samples_pretty, size_pretty)

    mids <- norm_cov$start + floor((norm_cov$end - norm_cov$start + 1) / 2)
    xlim <- c(min(norm_cov$start), max(norm_cov$end))
    jpeg(outfile, res = 300, width = 1800, height = 1800)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(6.1, 4.1, 4.1, 2.1))
    cnv_col <- if (cnv$svtype == "DEL") "red" else "blue"
    # plot background
    matplot(mids, as.matrix(norm_cov[, .SD, .SDcols = !patterns("chr|start|end")]),
            type = "l", lty = 1, col = "grey", lwd = 0.5,
            xlim = xlim,
            main = main,
            ylab = "Normalized Read Depth Ratio",
            xlab = "",
            xaxs = "i",
            xaxt = "n")
    # plot carriers
    matlines(mids, as.matrix(norm_cov[, .SD, .SDcols = cnv$samples_split]),
             lty = 1, lwd = 2,
             col = cnv_col)
    # lonely intervals - background
    bg_samples <- colnames(norm_cov)[!colnames(norm_cov) %in% c("chr", "start", "end", cnv$samples_split)]
    for (bg in norm_cov[, .SD, .SDcols = bg_samples]) {
        plot_singleton_intervals(bg, mids, "grey", 0.5)
    }
    # lonely intervals - carriers
    for (fg in norm_cov[, .SD, .SDcols = cnv$samples_split]) {
        plot_singleton_intervals(fg, mids, cnv_col, 2)
    }
    # add padding indicators
    ylim <- par("usr")[3:4]
    rect(xlim[[1]], ylim[[1]], cnv$cnv_start, ylim[[2]], col = "#FFAF0044", border = NA)
    rect(cnv$cnv_end, ylim[[1]], xlim[[2]], ylim[[2]], col = "#FFAF0044", border = NA)
    # add x-axis
    xticks <- axisTicks(xlim, log = FALSE, nint = 10)
    xlabs <- formatC(xticks, big.mark = ",", format = "d")
    axis(1, at = xticks, labels = xlabs, las = 2, cex.axis = 0.8)
    title(xlab = sprintf("%s Position (bp)", cnv$chr), line = 5)
    dev.off()
}

make_plot_path <- function(cnv, outdir) {
    plot_name <- sprintf("%s_%d-%d_%s_%s.jpg",
                         cnv$chr, cnv$cnv_start, cnv$cnv_end, cnv$vid, cnv$svtype)
    file.path(outdir, plot_name)
}

#' Add padding to CNV.
#'
#' The coordinates of the expanded region will be placed in `cnv$start` and
#' `cnv$end` and the original coordiantes will be placed in `cnv$cnv_start`
#' and `cnv$cnv_end`.
#'
#' @param cnv `list` Information of the CNV.
#' @returns `list` Updated CNV information.
expand_cnv <- function(cnv) {
    pad_size <- ceiling((cnv$end - cnv$start + 1L) * PAD_EXPANSION_FACTOR)
    cnv$cnv_start <- cnv$start
    cnv$cnv_end <- cnv$end
    cnv$start <- max(1L, cnv$cnv_start - pad_size)
    cnv$end <- cnv$cnv_end + pad_size

    cnv
}

# Main
visualize <- function(cnv, bincov_map, medians_map, sample_map, outdir) {
    cnv <- expand_cnv(cnv)
    batches <- unique(vapply(cnv$samples_split, \(x) gethash(sample_map, x), character(1)))
    if (any(sapply(batches, is.null))) {
        stop("some samples do not have a batch")
    }
    bincov_paths <- vapply(batches, \(x) gethash(bincov_map, x), character(1))
    if (any(sapply(bincov_paths, is.null))) {
        stop("some batches do not have a bincov path")
    }
    coverage <- load_cnv_coverage(cnv, bincov_paths)
    medians_paths <- vapply(batches, \(x) gethash(medians_map, x), character(1))
    if (any(sapply(medians_paths, is.null))) {
        stop("some batches do not have a medians path")
    }
    medians <- read_medians(medians_paths)
    norms <- normalize_cov(coverage, medians)

    # smooth
    norms[, names(.SD) := lapply(.SD, rollmean), .SDcols = !patterns("chr|start|end")]

    # sample intervals to reduce noise
    norms <- norms[select_spaced_intervals(.N), ]

    plot_path <- make_plot_path(cnv, outdir)
    make_plot(cnv, norms, plot_path)
}

message("reading CNVs")
records <- read_cnvs(opts$cnvs)
message("reading sample batches")
sample_map <- read_keyed_tsv(opts$sample_batches)
message("reading coverage paths")
bincov_map <- read_keyed_tsv(opts$coverage_paths)
message("reading median paths")
medians_map <- read_keyed_tsv(opts$medians_paths)
dir.create(opts$output)

for (i in seq_len(nrow(records))) {
    rec <- as.list(records[i, ])
    rec$samples_split <- strsplit(rec$samples, split = ",", fixed = TRUE)[[1]]
    message(sprintf("plotting '%s' (%s:%d-%d)", rec$vid, rec$chr, rec$start, rec$end))
    visualize(rec, bincov_map, medians_map, sample_map, opts$output)
}
message("done")
