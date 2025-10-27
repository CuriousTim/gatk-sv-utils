TABIX_MAX_SEQLEN <- 536870912L
ROLLING_MEAN_WINDOW <- 21
MAX_BACKGROUND <- 200L
BIN_PLOT_FRACTION <- 0.1

read_sd <- function(path) {
  tmp <- fread(path, header = FALSE, select = 1:3, sep = "\t")
  tmp <- tmp[V1 %in% paste0("chr", c(1:22, "X", "Y")), ]
  tmp <- unique(tmp)
  reduce(GRanges(tmp[["V1"]], IRanges(tmp[["V2"]], tmp[["V3"]])))
}

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

read_medians_file <- function(path) {
    lines <- readLines(path, n = 2L, ok = FALSE)
    ids <- strsplit(lines[[1]], split = "\t", fixed = TRUE)[[1]]
    if (anyDuplicated(ids) != 0) {
        stop("duplicate samples in medians file")
    }
    medians <- as.double(strsplit(lines[[2]], split = "\t", fixed = TRUE)[[1]])

    names(medians) <- ids

    medians
}

query_bincov <- function(con, header, contig, start, end) {
    lines <- scanTabix(con, param = GRanges(contig, IRanges(start, end)))[[1]]
    if (length(lines) == 0) {
        return(NULL)
    }

    coltypes <- c(list(character(), integer(), integer()),
                  lapply(seq_len(length(header) - 3), \(x) double()))
    parsed <- scan(text = lines, what = coltypes, nmax = length(lines),
                   sep = "\t", quiet = TRUE)
    names(parsed) <- header
    tmp <- as.data.table(parsed)
    tmp[, start := start + 1L]
}

normalize_bincov <- function(x, medians) {
    samples <- colnames(x)
    if (!all(samples %in% names(medians))) {
        stop("all samples in the bincov matrix must have a median coverage")
    }

    medians <- medians[samples]
    assay(x, 1) <- sweep(assay(x, 1), 2, medians, "/", check.margin = FALSE)

    x
}

rollmean <- function(x) {
    binwidth <- x[, end - start + 1]
    if (length(unique(binwidth)) != 1) {
        stop("bin widths must be identical")
    }
    binwidth <- binwidth[[1]]
    window <- floor(ROLLING_MEAN_WINDOW / 2)
    bins <- x[["start"]] - x[1, start] / binwidth

    assay(x, 1) <- as.data.table(lapply(assay(x, 1),
                                        \(x) slider::slide_index_mean(x, bins, before = window, after = window)))

    x
}

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

get_sd_overlaps <- function(sds, chr, start, end) {
    query <- GRanges(chr, IRanges(start, end))
    ovp <- findOverlaps(query, sds)
    sds[subjectHits(ovp)]
}

# Main ------------------------------------------------------------------------

# Usage:
# Rscript visualize_gd.R <gd_regions> <sd_regions> <bincov> <medians> <samples> <min_shift> <pad> <outdir>
# gd_regions  genomic disorder regions
# sd_regions  segemental duplications
# bincov      binned coverage matrix
# medians     coverage medians
# samples     list of samples to check for CNVs
# min_shift   minimum amount a sample's read depth ratio must shifted from 1 to
#             considered a CNV carrier
# pad         fraction by which the genomic disorder region should be expanded
#             for plotting
# outdir      output directory

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(slider))

argv <- commandArgs(trailingOnly = TRUE)

gd_regions <- fread(argv[[1]], header = TRUE, sep = "\t")
sd_regions <- read_sd(argv[[2]])
bincov_con <- TabixFile(argv[[3]])
cov_medians <- read_medians_file(argv[[4]])
samples <- readLines(argv[[5]])
min_shift <- as.double(argv[[6]])
pad <- as.double(argv[[7]])
outdir <- argv[[8]]
dir.create(outdir)

if (pad < 0) {
    stop("padding must be a non-negative number")
}
header <- read_bincov_header(bincov_con)
if (!all(samples %in% header)) {
    stop("all requested samples must be in the bincov matrix")
}

if (!all(samples %in% names(cov_medians))) {
    stop("all requested samples must be in the median coverage file")
}

pad_size <- gd_regions[, ceiling((end_GRCh38 - start_GRCh38 + 1L) * ..pad)]
expanded_gd <- gd_regions[, list(chr = chr, start = pmax(1L, start_GRCh38 - pad_size), end = pmin(TABIX_MAX_SEQLEN, end_GRCh38 + pad_size))]
for (i in seq_len(nrow(gd_regions))) {
    chr <- gd_regions[i, chr]
    gdstart <- gd_regions[i, start_GRCh38]
    gdend <- gd_regions[i, end_GRCh38]
    gdid <- gd_regions[i, GD_ID]
    svtype <- gd_regions[i, svtype]
    qstart <- expanded_gd[i, start]
    qend <- expanded_gd[i, end]
    message(sprintf("checking for genomic disorder at %s:%d-%d (%s)", chr, gdstart, gdend, svtype))
    if (svtype != "DUP" && svtype != "DEL") {
        stop("SV type must be 'DUP' or 'DEL'")
    }
    bc <- query_bincov(bincov_con, header, chr, qstart, qend)
    if (is.null(bc) || nrow(bc) < 2) {
        message("giving up due to insufficient number of bincov intervals")
        next
    }

    bc_coords <- bc[, list(chr, start, end)]
    bc[, c("chr", "start", "end") := list(NULL)]
    bc <- bc[, ..samples]

    # normalize
    norms <- bc[, mapply(`/`, .SD, ..cov_medians[names(.SD)], SIMPLIFY = FALSE)]

    # potential carriers
    rdr_medians <- vapply(norms, median, double(1))
    if (svtype == "DUP") {
        carriers <- rdr_medians[rdr_medians >= 1 + min_shift]
    } else {
        carriers <- rdr_medians[rdr_medians <= 1 - min_shift]
    }
    if (length(carriers) == 0) {
        message("no potential CNV carriers")
        next
    }
    bg_samples <- setdiff(samples, names(carriers))
    bg_samples_to_plot <- sample(bg_samples, min(length(bg_samples), MAX_BACKGROUND))

    if (nrow(bc_coords) > ROLLING_MEAN_WINDOW) {
        binwidth <- bc_coords[, end - start + 1]
        if (length(unique(binwidth)) != 1) {
            stop("bincov interval widths must be identical")
        }
        binwidth <- binwidth[[1]]
        window <- floor(ROLLING_MEAN_WINDOW / 2)
        bins <- (bc_coords[["start"]] - bc_coords[1, start]) / binwidth
        norms[, names(.SD) := lapply(.SD, \(x) slide_index_mean(x, ..bins, before = ..window, after = ..window))]
    }

    mids <- bc_coords[["start"]] + floor(bc_coords[, end - start + 1] / 2)
    bins_to_plot <- spaced_intervals(length(mids))
    size_pretty <- format_size(gdend - gdstart + 1)
    sd_overlaps <- get_sd_overlaps(sd_regions, chr, qstart, qend)
    for (j in seq_along(carriers)) {
        carrier_id <- names(carriers)[[j]]
        cnv_col <- if (svtype == "DEL") "red" else "blue"
        if (nchar(carrier_id) > 30) {
            carrier_pretty <- paste0(substr(carrier_id, 1, 30), "...")
        } else {
            carrier_pretty <- carrier_id
        }
        main <- sprintf("%s (hg38)\n%s (%s)", gdid, carrier_pretty, size_pretty)
        plot_name <- sprintf("%s_%d-%d_%s_%s_%s.jpg",
                             chr, gdstart, gdend, gdid, svtype, carrier_id)
        plot_path <- file.path(outdir, plot_name)
        jpeg(plot_path, res = 100, width = 960, height = 540)
        par(mar = c(3.1, 4.1, 4.1, 2.1))
        plot(NULL, main = main, xlim = range(mids), ylim = c(0, 3), ylab = "Normalized Read Depth Ratio", xlab = "", xaxs = "i", xaxt = "n")
        mids_to_plot <- mids[bins_to_plot]
        rd_to_plot <- norms[bins_to_plot, ]
        rd_bg <- rd_to_plot[, ..bg_samples_to_plot]
        for (k in seq_along(rd_bg)) {
            lines(mids_to_plot, rd_bg[[k]], col = "grey", lwd = 0.5)
        }
        lines(mids_to_plot, rd_to_plot[[carrier_id]], lty = 1, lwd = 2, col = cnv_col)

        axlims <- par("usr")
        if (length(sd_overlaps) > 0) {
            rect(start(sd_overlaps), 0.1, end(sd_overlaps), 0.2, col = "brown4", border = NA)
        }

        rect(axlims[[1]], axlims[[3]], gdstart, axlims[[4]], col = "#FFAF0044", border = NA)
        rect(gdend, axlims[[3]], axlims[[2]], axlims[[4]], col = "#FFAF0044", border = NA)
        box(lwd = 2)
        xticks <- axisTicks(axlims[1:2], log = FALSE, nint = 10)
        xlabs <- formatC(xticks, big.mark = ",", format = "d")
        axis(1, at = xticks, labels = FALSE, las = 2, cex.axis = 0.8, tcl = -1)
        mtext(xlabs, 1, at = xticks, adj = 1.05, cex = 0.6)
        title(xlab = sprintf("%s Position (bp)", chr), line = 2)
        dev.off()
    }
}
