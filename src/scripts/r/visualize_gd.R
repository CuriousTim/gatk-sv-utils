# Visualize Genomic Disorder Regions
#
# Usage:
# Rscript visualize_gd.R [options] <gd_regions> <sd_regions> <bincov> <medians> \
#   <samples> <sex_ploidy> <outdir>
# gd_regions      genomic disorder regions to visualize
# sd_regions      segmental duplications
# bincov          binned coverage matrix
# medians         coverage medians
# samples         list of samples to check for CNVs
# sex_ploidy      sex choromosome ploidy table, 0 for unknown, 1 for male, 2 for female
# outdir          output directory
#
# options
# --min-shift             minimum amount a sample's read depth ratio must shifted from 1 to
#                         considered a CNV carrier
# --pad                   fraction by which the genomic disorder region should be expanded
#                         for plotting
# --max-calls-per-sample  maximum number of calls per sample
# --violators             samples that had more than the max number of calls

# maximum length of a sequence that can be indexed with tabix
TABIX_MAX_SEQLEN <- 536870912L
# the window in number of bins to use for smoothing the binned coverage values before plotting
SMOOTHING_WINDOW <- 31
# maximum number of background samples to plot
MAX_BACKGROUND <- 200L
# the fraction of bins to plot, low mostly to reduce computation
BIN_PLOT_FRACTION <- 0.05
# the fraction of the region to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_PROP <- 0.01
# the minimum number of windows to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_MIN_BINS <- 101
# the assumed width of each bin in the bincov matrix
BIN_WIDTH <- 100
# minimum number of bins (assuming 100 bp bins) overlapping a NAHR GD region
# required to make a CNV call
NAHR_WINDOW_MIN_BINS <- 50

HG38_PRIMARY_CONTIGS <- paste0("chr", c(1:22, "X", "Y"))

assert_hg38_contigs <- function(x) {
    if (any(!x %in% HG38_PRIMARY_CONTIGS)) {
        stop("only chromosomes chr1-22,X,Y are permitted")
    }
}

assert_valid_contig_pos <- function(x) {
    if (any(x <= 0L | x > TABIX_MAX_SEQLEN)) {
        stop(sprintf("only chromosome positions between 1 and %d inclusive are permitted", TABIX_MAX_SEQLEN))
    }
}

assert_positive_range <- function(start, end) {
    if (any(end < start)) {
        stop("only positive genomic ranges are permitted")
    }
}

assert_valid_svtypes <- function(x) {
    if (any(!x %in% c("DUP", "DEL"))) {
        stop("only DUP and DEL SV types are permitted")
    }
}

assert_valid_nahr <- function(x) {
    if (any(!x %in% c("yes", "no"))) {
        stop("NAHR values must be 'yes' or 'no'")
    }
}

assert_valid_gd_terminal_type <- function(x) {
    if (any(!x %in% c("p", "q", "no"))) {
        stop("GD terminal type must be 'p', 'q', or 'no'")
    }
}

# Add an entry to the plots store to record that `s` has a plot at `path` and
# increment the count of plots for that sample. Plots from GDs at the same
# cluster count once per sample.
add_plot_to_store <- function(h, s, path, cluster) {
    val <- gethash(h, s)
    if (is.null(val)) {
        val <- list(clusters = hashtab(), paths = character(), count = 0L)
    }
    val[["paths"]] <- append(val[["paths"]], path)

    if (!is.na(cluster) && is.null(gethash(val[["clusters"]], cluster))) {
        sethash(val[["clusters"]], cluster, NULL)
        val[["count"]] <- val[["count"]] + 1L
    }

    if (is.na(cluster)) {
        val[["count"]] <- val[["count"]] + 1L
    }

    sethash(h, s, val)
}

# Is the sample a violator?
is_sample_violator <- function(h, s, max_plots) {
    val <- gethash(h, s)
    if (is.null(val)) {
        return(FALSE)
    }

    val[["count"]] > max_plots
}

# Remove all the plots from violators and record which ones were removed.
remove_violators_plots <- function(k, v, max_plots, con) {
    if (v[["count"]] > max_plots) {
        writeLines(k, con)
        file.remove(v[["paths"]])
    }
}

# Read the segdups in BED file.
read_sd <- function(path) {
    tmp <- fread(path, header = FALSE, select = 1:3, sep = "\t")

    tmp <- tmp[V1 %in% HG38_PRIMARY_CONTIGS, ]
    tmp[, c("V2", "V3") := list(as.integer(V2), as.integer(V3))]
    tmp[, V2 := V2 + 1L]

    assert_positive_range(tmp$V2, tmp$V3)

    tmp <- unique(tmp)
    reduce(GRanges(tmp[["V1"]], IRanges(tmp[["V2"]], tmp[["V3"]])))
}

# Read the bincov matrix header.
read_bincov_header <- function(con) {
    header <- headerTabix(TabixFile(path(con)))[["header"]]
    if (is.null(header) || length(header) == 0) {
        stop("bincov matrix is missing a header")
    }

    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) < 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop("bincov matrix has an invalid header")
    }
    header[1:3] <- c("chr", "start", "end")

    header
}

# Read the sample median coverages file.
read_medians_file <- function(path) {
    lines <- readLines(path, n = 2L, ok = FALSE)
    ids <- strsplit(lines[[1]], split = "\t", fixed = TRUE)[[1]]
    if (anyDuplicated(ids) != 0) {
        stop("duplicate samples in medians file")
    }
    medians <- as.double(strsplit(lines[[2]], split = "\t", fixed = TRUE)[[1]])

    stopifnot("median coverage values must be positive" = medians > 0)

    setNames(medians, ids)
}

# Read the file with samples to check for GDs.
read_samples_file <- function(path) {
    tmp <- readLines(path)
    if (length(tmp) == 0) {
        stop("at least one sample to check for CNVs must be given")
    }

    tmp
}

# Read the sex chromosome ploidy table.
read_sex_ploidy <- function(path) {
    tmp <- fread(path, header = FALSE, sep = "\t",
                 col.names = c("sample_id", "sex"),
                 colClasses = c("character", "integer"))
    if (any(!tmp$sex %in% c(0L, 1L, 2L))) {
        stop("permitted values for sex are 0, 1, and 2")
    }
    tmp[sex == 0, sex := NA_integer_]

    tmp
}

# Read the genomic disorders table
read_gd <- function(path) {
    tmp <- fread(path, header = TRUE, sep = "\t", strip.white = FALSE,
                 na.strings = "",
                 colClasses = c("character", "integer", "integer", "character", "character", "character", "character", "character"))
    req_header <- c("chr", "start_GRCh38", "end_GRCh38", "GD_ID", "svtype", "NAHR", "terminal", "cluster")
    disp_header <- paste0(req_header, collapse = ", ")
    if (!identical(colnames(tmp), req_header)) {
        stop(sprintf("genomic disorder regions file must have header: '%s'",
                     disp_header))
    }

    assert_hg38_contigs(tmp$chr)
    assert_valid_contig_pos(tmp$start_GRCh3)
    assert_valid_contig_pos(tmp$end_GRCh38)
    assert_positive_range(tmp$start_GRCh38, tmp$end_GRCh38)
    assert_valid_svtypes(tmp$svtype)
    assert_valid_nahr(tmp$NAHR)
    assert_valid_gd_terminal_type(tmp$terminal)

    tmp[, NAHR := NAHR == "yes"]

    tmp
}

# Return a function that can be used to query a bincov matrix
make_bincov_getter <- function(con, header, samples, medians) {
    ordered_medians <- medians[samples]

    f <- function(chr, start, end) {
        lines <- scanTabix(con, param = GRanges(chr, IRanges(start, end)))[[1]]
        if (length(lines) == 0) {
            return(NULL)
        }

        coltypes <- c(list(character(), integer(), integer()),
                      lapply(seq_len(length(header) - 3), \(x) double()))
        parsed <- scan(text = lines, what = coltypes, nmax = length(lines),
                       sep = "\t", quiet = TRUE)
        names(parsed) <- header
        tmp <- as.data.table(parsed)
        # bincov matrix coordinates are 0-start
        tmp[, start := start + 1L]

        regions <- GRanges(tmp$chr, IRanges(tmp$start, tmp$end))
        tmp[, c("chr", "start", "end") := list(NULL)]

        tmp[, names(.SD) := mapply(`/`, .SD, ..ordered_medians, SIMPLIFY = FALSE), .SDcols = samples]

        list(regions = regions, norm_bc = tmp)
    }
}

# Compute the expected read depth ratio for some samples for a chromosome given
# the sex of the samples (NA = unknown, 1 = male, 2 = female).
get_expected_rdr <- function(chr, sex_ploidy) {
    if (chr == "chrX" || chr == "chrY") {
        if (chr == "chrY") {
            sex_ploidy[sex_ploidy == 2] <- 0L
        }
        sex_ploidy / 2
    } else {
        rep_len(1, length(sex_ploidy))
    }
}

# Compute the sliding window size for non-NAHR GD regions.
non_nahr_window_size <- function(gdstart, gdend) {
    max(NON_NAHR_WINDOW_MIN_BINS, round((gdend - gdstart  + 1) * NON_NAHR_WINDOW_PROP / BIN_WIDTH))
}

# Return a function that compares the expected and actual read depth ratios and
# determines if the actual is sufficiently shifted from the expected based on
# SV type and a minimum shift. The expected and actual arguments should be two
# parallel vectors for multiple samples or a scalar for the expected and vector
# of multiple actuals for a single sample.
make_rdr_shift_filter <- function(svtype, min_shift) {
    op <- if (svtype == "DUP") {
        function(expected, actual) { actual >= expected + min_shift }
    } else {
        function(expected, actual) { actual <= expected - min_shift }
    }

    function(expected, actual) {
        !is.nan(actual) & !is.na(actual) & op(expected, actual)
    }
}

# Predict CNVs carriers at GD regions.
# For NAHR-mediated GDs, the test for CNVs is a shifted median read depth ratio
# over the entire region.
#
# For non-NAHR-mediated GDs, the test for CNVs depends on the GD terminal type:
# For p-arm GDs, the first window (smaller genomic position) in the region must
# have a shifted median read depth ratio.
# For q-arm GDs, the last window (larger genomic position) in the region must
# have a shifted median read depth ratio.
# For non-terminal GDs, at least one window in the region must have a shifted
# median read depth ratio.
predict_gd_carriers <- function(gd, bc, sex_ploidy, min_shift) {
    expected_rdr <- get_expected_rdr(gd$chr, sex_ploidy)
    filter_rdr <- make_rdr_shift_filter(gd$svtype, min_shift)

    if (gd$NAHR) {
        m <- vapply(bc, \(x) median(x, na.rm = TRUE), double(1))
        names(m[filter_rdr(expected_rdr, m)])
    } else {
        window_size <- non_nahr_window_size(gd$start_GRCh38, gd$end_GRCh38)
        if (window_size %% 2 == 0) {
            window_size <- window_size + 1
        }
        switch(gd$terminal,
               "p" = {
                    m <- bc[seq_len(window_size), vapply(.SD, \(x) median(x, na.rm = TRUE), double(1))]
                    names(m[filter_rdr(expected_rdr, m)])
               },
               "q" = {
                    m <- bc[seq(nrow(bc) - window_size, nrow(bc)), vapply(.SD, \(x) median(x, na.rm = TRUE), double(1))]
                    names(m[filter_rdr(expected_rdr, m)])
               },
               {
                    # ignore partial windows
                    m <- lapply(bc, \(x) runmed(x, window_size))
                    half_window <- window_size %/% 2
                    full_windows <- seq(1 + half_window, nrow(bc) - half_window)
                    filtered <- mapply(\(x, y) filter_rdr(x, y[full_windows]), expected_rdr, m, SIMPLIFY = FALSE)
                    pass <- vapply(filtered, any, logical(1))
                    names(m)[pass]
               })
    }
}

# Generate a vector of best-effort equally spaced integers from 1 to n.
spaced_intervals <- function(n) {
    nintervals <- ceiling(n * BIN_PLOT_FRACTION)
    if (nintervals < 2) {
        return(seq_len(n))
    }
    gaps <- nintervals - 1L
    total_gap_size <- n - nintervals
    remaining_pad <- total_gap_size %% gaps
    pads <- rep(floor(total_gap_size / gaps), gaps)
    i <- seq_len(remaining_pad)
    pads[i] <- pads[i] + 1L
    steps <- c(1, pads + 1L)

    cumsum(steps)
}

# Expand genomic disorder regions by a fraction of the region size. Modifies
# the input and returns it.
expand_gd_regions <- function(x, prop) {
    pad_size <- x[, ceiling((end_GRCh38 - start_GRCh38 + 1L) * prop)]
    x[, c("qstart ", "qend") := list(pmax(1L, start_GRCh38 - ..pad_size), pmin(TABIX_MAX_SEQLEN, end_GRCh38 + ..pad_size))]

    x
}

# Return a function that will make plots.
make_plotter <- function(gd, regions, bincov, bg_samples, sd_regions) {
    mids <- start(regions) + (width(regions) %/% 2)
    bins_to_plot <- spaced_intervals(length(mids))
    mids_to_plot <- mids[bins_to_plot]
    bincov_to_plot <- bincov[bins_to_plot, ]
    size_pretty <- format_size(gd$end_GRCh38 - gd$start_GRCh38 + 1)
    carrier_col <- if (gd$svtype == "DEL") "red" else "blue"
    bg_bincov <- bincov_to_plot[, ..bg_samples]
    sd_ovp <- sd_regions[queryHits(findOverlaps(sd_regions, GRanges(gd$chr, IRanges(gd$qstart, gd$qend))))]

    function(carrier, path) {
        carrier_pretty <- format_sample_id(carrier)
        main <- sprintf("%s (hg38)\n%s (%s)", gd$GD_ID, carrier_pretty, size_pretty)
        old_par <- par(no.readonly = TRUE)
        jpeg(path, res = 100, width = 960, height = 540)
        par(mar = c(3.1, 4.1, 4.1, 2.1))
        plot(NULL, main = main, xlim = range(mids), ylim = c(0, 3),
             ylab = "Normalized Read Depth Ratio", xlab = "", xaxs = "i", xaxt = "n")
        for (i in seq_along(bg_bincov)) {
            lines(mids_to_plot, bg_bincov[[i]], col = "grey", lwd = 0.5)
        }
        lines(mids_to_plot, bincov_to_plot[[carrier]], col = carrier_col, lwd = 0.5)
        if (length(sd_ovp) > 0) {
            rect(start(sd_ovp), 0.1, end(sd_ovp), 0.2, col = "brown4", border = NA)
        }

        axlims <- par("usr")
        rect(axlims[[1]], axlims[[3]], gd$start_GRCh38, axlims[[4]], col = "#FFAF0044", border = NA)
        rect(gd$end_GRCh38, axlims[[3]], axlims[[2]], axlims[[4]], col = "#FFAF0044", border = NA)
        box(lwd = 2)
        xticks <- axisTicks(axlims[1:2], log = FALSE, nint = 10)
        xlabs <- formatC(xticks, big.mark = ",", format = "d")
        axis(1, at = xticks, labels = FALSE, las = 2, cex.axis = 0.8, tcl = -1)
        mtext(xlabs, 1, at = xticks, adj = 1.05, cex = 0.6)
        title(xlab = sprintf("%s Position (bp)", gd$chr), line = 2)
        par(old_par)
        dev.off()
    }
}

# Format a length in bases to a pretty string.
format_size <- function(size) {
    if (size <= 1000) {
        size_pretty <- sprintf("%d b", size)
    } else if (size <= 100000) {
        size_pretty <- sprintf("%.2f Kb", size / 1000)
    } else {
        size_pretty <- sprintf("%.2f Mb", size / 1000000)
    }

    size_pretty
}

# Format a sample ID into a pretty string for plotting.
format_sample_id <- function(x) {
    if (nchar(x) > 30) {
        paste0(substr(x, 1, 30), "...")
    } else {
        x
    }
}

# Validate the script arguments.
validate_args <- function(x) {
    if (!is.finite(x$min_shift) || x$min_shift < 0) {
        stop("minimum shift must be a non-negative number")
    }

    if (!is.finite(x$pad) || x$pad < 0) {
        stop("padding must be a non-negative number")
    }

    if (x$max_calls_per_sample <= 0) {
        stop("max calls per sample must be a positive number")
    }

    x
}

# Parse the script arguments.
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    pos_args <- vector("list", 7)
    opts <- list(min_shift = 0.3, pad = 0.5,
                 max_calls_per_sample = 3, violators = NULL)
    i <- 1
    j <- 1
    repeat {
        if (i > length(args)) {
            break
        }

        if (args[[i]] == "--min-shift") {
            i <- i + 1
            opts$min_shift <- as.double(args[[i]])
        } else if (args[[i]] == "--pad") {
            i <- i + 1
            opts$pad <- as.double(args[[i]])
        } else if (args[[i]] == "--max-calls-per-sample") {
            i <- i + 1
            opts$max_calls_per_sample <- trunc(as.double(args[[i]]))
        } else if (args[[i]] == "--violators") {
            i <- i + 1
            opts$violators <- args[[i]]
        } else {
            pos_args[[j]] <- args[[i]]
            j <- j + 1
        }

        i <- i + 1
    }

    if (any(sapply(pos_args, is.null))) {
        stop("incorrect number of arguments")
    }
    names(pos_args) <- c("gd_regions", "sd_regions", "bincov", "medians",
                         "samples", "sex_ploidy", "outdir")

    validate_args(append(pos_args, opts))
}

# Main ------------------------------------------------------------------------

argv <- parse_args()

message("loading libraries")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))

message("reading genomic disorder regions")
gd_regions <- read_gd(argv$gd_regions)
message("reading segmental duplication regions")
sd_regions <- read_sd(argv$sd_regions)
message("reading median coverages file")
cov_medians <- read_medians_file(argv$medians)
message("reading list of target samples")
target_samples <- read_samples_file(argv$samples)
message("reading sex chromosome ploidy table")
sex_ploidy <- read_sex_ploidy(argv$sex_ploidy)

message("reading bincov matrix header")
bincov_con <- TabixFile(argv$bincov)
bincov_header <- read_bincov_header(bincov_con)

if (!all(target_samples %in% bincov_header)) {
    stop("all requested samples must be in the bincov matrix")
}

if (!all(target_samples %in% names(cov_medians))) {
    stop("all requested samples must be in the median coverage file")
}
dir.create(argv$outdir)

sex_ploidy <- sex_ploidy[target_samples, nomatch = NA, on = "sample_id"]
setkey(sex_ploidy, sample_id)

# set up hash map to store information on how many plots have been made per
# sample
plots_store <- hashtab()

gd_regions <- expand_gd_regions(gd_regions, argv$pad)
active_samples <- target_samples
query_bincov <- make_bincov_getter(bincov_con, bincov_header, target_samples, cov_medians)

message("getting to work")
for (i in seq_len(nrow(gd_regions))) {
    if (length(active_samples) == 0) {
        message("all samples have produced too many GD predictions")
        break
    }

    cur <- as.list(gd_regions[i, ])
    message(sprintf("checking for genomic disorder at %s:%d-%d (%s)",
                    cur$chr, cur$start_GRCh38, cur$end_GRCh38, cur$GD_ID))

    bc <- query_bincov(cur$chr, cur$qstart, cur$qend)
    if (is.null(bc)) {
        message("no bincov intervals overlap GD region")
        next
    }

    gd_idx <- queryHits(findOverlaps(bc$regions, GRanges(cur$chr, IRanges(cur$start_GRCh38, cur$end_GRCh38))))
    if ((cur$NAHR && length(gd_idx) < NAHR_WINDOW_MIN_BINS) || length(gd_idx) < NON_NAHR_WINDOW_MIN_BINS) {
        message("insufficient bincov intervals overlap GD region")
        next
    }

    sd_idx <- queryHits(findOverlaps(bc$regions, sd_regions))
    active_bc <- bc$norm_bc[, ..active_samples]
    # exclude bins overlapping segdups from median calculations
    active_bc[sd_idx, names(.SD) := list(NA_real_)]

    carriers <- predict_gd_carriers(cur, active_bc[gd_idx, ], sex_ploidy[active_samples, sex], argv$min_shift)

    if (length(carriers) == 0) {
        message("no potential CNV carriers")
        next
    } else {
        message(sprintf("found %d potential CNV carriers", length(carriers)))
    }

    bg_samples <- sample(target_samples, min(length(target_samples), MAX_BACKGROUND))

    if (length(bc$regions) > SMOOTHING_WINDOW) {
        bc$norm_bc[, names(.SD) := lapply(.SD, \(x) runmed(x, SMOOTHING_WINDOW, na.action = "fail"))]
    }

    message("making plots")
    make_plot <- make_plotter(cur, bc$regions, bc$norm_bc, bg_samples, sd_regions)
    for (carrier in carriers) {
        plot_name <- sprintf("%s_%d-%d_%s_%s_%s.jpg",
                             cur$chr, cur$start_GRCh38, cur$end_GRCh38, cur$GD_ID, cur$svtype, carrier)
        plot_path <- file.path(argv$outdir, plot_name)
        make_plot(carrier, plot_path)
        add_plot_to_store(plots_store, carrier, plot_path, cur$cluster)
    }
    violators <- carriers[vapply(carriers, \(x) is_sample_violator(plots_store, x, argv$max_calls_per_sample), logical(1))]
    active_samples <- setdiff(active_samples, violators)
}

violators_con <- file(if (is.null(argv$violators)) nullfile() else argv$violators)
maphash(plots_store, \(k, v) remove_violators_plots(k, v, argv$max_calls_per_sample, violators_con))
message("done")
