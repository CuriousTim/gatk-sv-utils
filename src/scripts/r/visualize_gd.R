# Visualize Genomic Disorder Regions
#
# Usage:
# Rscript visualize_gd.R [options] <gd_regions> <sd_regions> <bincov> <medians> \
#   <samples> <ploidy> <outdir>
# gd_regions      genomic disorder regions to visualize
# sd_regions      segmental duplications
# bincov          binned coverage matrix
# medians         coverage medians
# samples         list of samples to check for CNVs
# ploidy          ploidy table, 0 for unknown, 1 for male, 2 for female
# outdir          output directory
#
# options
# --min-shift             minimum amount a sample's read depth ratio must shifted from 1 to
#                         considered a CNV carrier
# --pad                   fraction by which the genomic disorder region should be expanded
#                         for plotting
# --min-shifted-bins      number of consecutive bins that must have a mean read depth
#                         shifted from 1, overriding the default of mean over the
#                         entire region
# --max-calls-per-sample  maximum number of calls per sample
# --violators             samples that had more than the max number of calls

TABIX_MAX_SEQLEN <- 536870912L
ROLLING_MEDIAN_WINDOW <- 31
MAX_BACKGROUND <- 200L
BIN_PLOT_FRACTION <- 0.05

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

read_ploidy <- function(path) {
    tmp <- fread(path, header = FALSE, sep = "\t",
                 col.names = c("sample_id", "sex"),
                 colClasses = c("character", "integer"))
    if (any(!tmp$sex %in% c(0L, 1L, 2L))) {
        stop("permitted values for sex are 0, 1, and 2")
    }
    tmp[sex == 0, sex := NA_integer_]
    tmp
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

predict_carriers <- function(bc, chr, ploidy, window, svtype, min_shift) {
    if (nrow(bc) < 2) {
        return(NA_character_)
    }

    if (chr == "chrY") {
        males <- ploidy[sex == 1L, sample_id]
        if (length(males) == 0) {
            return(character())
        }
        bc <- bc[, ..males]
        if (svtype == "DUP") {
            op <- function(x) { x >= 0.5 + min_shift}
        } else {
            op <- function(x) { x <= 0.5 - min_shift}
        }
    } else if (chr == "chrX") {
        with_ploidy <- ploidy[!is.na(sex), sample_id]
        if (length(with_ploidy) == 0) {
            return(character())
        }
        bc <- bc[, colnames(bc) %in% with_ploidy, with = FALSE]
        expected_rdr <- ploidy[colnames(bc), sex] / 2
        if (svtype == "DUP") {
            op <- function(x) { x >= expected_rdr + min_shift }
        } else {
            op <- function(x) { x <= expected_rdr - min_shift }
        }
    } else {
        if (svtype == "DUP") {
            op <- function(x) { x >= 1 + min_shift }
        } else {
            op <- function(x) { x <= 1 - min_shift}
        }
    }

    if (!is.null(window)) {
        if (nrow(bc) < window) {
            return(NA_character_)
        }
        keep <- seq.int(ceiling(window / 2), nrow(bc) - floor(window / 2))
        roll <- lapply(bc, \(x) runmed(x, window)[keep])
        pass <- vapply(roll, \(x) any(op(x), na.rm = TRUE), logical(1))
        pass <- pass[pass == TRUE]
    } else {
        m <- vapply(bc, median, double(1))
        pass <- m[op(m)]
    }

    names(pass)
}

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    pos_args <- vector("list", 7)
    opts <- list(min_shift = 0.3, pad = 0.5, min_shifted_bins = NULL, max_calls_per_sample = 3, violators = NULL)
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
        } else if (args[[i]] == "--min-shifted-bins") {
            i <- i + 1
            opts$min_shifted_bins <- as.integer(args[[i]])
        } else if (args[[i]] == "--max-calls-per-sample") {
            i <- i + 1
            opts$max_calls_per_sample <- as.integer(args[[i]])
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

    append(pos_args, opts)
}

add_plot_to_store <- function(h, s, path, cluster) {
    val <- gethash(h, s)
    if (is.null(val)) {
        val <- list(clusters = hashtab(), paths = character(), count = 0L)
    }
    val[["paths"]] <- append(val[["paths"]], path)

    if (!is.na(cluster) && is.null(gethash(val[["clusters"]], cluster))) {
        val[["clusters"]] <- sethash(val[["clusters"]], cluster, NULL)
        val[["count"]] <- val[["count"]] + 1L
    }

    if (is.na(cluster)) {
        val[["count"]] <- val[["count"]] + 1L
    }
}

remove_outliers <- function(k, v, max_count, con) {
    if (v[["count"]] > max_count) {
        writeLines(sprintf("%s\t%d", k, v[["count"]]), con)
        file.remove(v[["paths"]])
    }
}

# Main ------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))

argv <- parse_args()

gd_regions <- fread(argv[[1]], header = TRUE, sep = "\t", strip.white = FALSE, na.strings = "")
sd_regions <- read_sd(argv[[2]])
bincov_con <- TabixFile(argv[[3]])
cov_medians <- read_medians_file(argv[[4]])
samples <- readLines(argv[[5]])
ploidy <- read_ploidy(argv[[6]])
outdir <- argv[[7]]
min_shift <- argv[["min_shift"]]
pad <- argv[["pad"]]
min_shifted_bins <- argv[["min_shifted_bins"]]
max_calls_per_sample <- argv[["max_calls_per_sample"]]
violators <- argv[["violators"]]

dir.create(outdir)

if (!is.finite(pad) || pad < 0) {
    stop("padding must be a non-negative number")
}
if (!is.finite(min_shift) || min_shift < 0) {
    stop("min shift must be a non-negative number")
}
if (!is.null(min_shifted_bins) && (!is.finite(min_shifted_bins) || min_shifted_bins < 2 || min_shifted_bins %% 2 == 0)) {
    stop("min shifted bins must be an odd integer greater than 2")
}
if (length(samples) == 0) {
    stop("at least one sample must be given")
}

header <- read_bincov_header(bincov_con)
if (!all(samples %in% header)) {
    stop("all requested samples must be in the bincov matrix")
}

if (!all(samples %in% names(cov_medians))) {
    stop("all requested samples must be in the median coverage file")
}

ploidy <- ploidy[samples, nomatch = NA, on = "sample_id"]
setkey(ploidy, sample_id)

plots_store <- hashtab()

pad_size <- gd_regions[, ceiling((end_GRCh38 - start_GRCh38 + 1L) * ..pad)]
expanded_gd <- gd_regions[, list(chr = chr, start = pmax(1L, start_GRCh38 - pad_size), end = pmin(TABIX_MAX_SEQLEN, end_GRCh38 + pad_size))]
for (i in seq_len(nrow(gd_regions))) {
    chr <- gd_regions[i, chr]
    gdstart <- gd_regions[i, start_GRCh38]
    gdend <- gd_regions[i, end_GRCh38]
    gdid <- gd_regions[i, GD_ID]
    svtype <- gd_regions[i, svtype]
    cluster <- gd_regions[i, cluster]
    qstart <- expanded_gd[i, start]
    qend <- expanded_gd[i, end]
    message(sprintf("checking for genomic disorder at %s:%d-%d (%s)", chr, gdstart, gdend, svtype))
    if (svtype != "DUP" && svtype != "DEL") {
        stop("SV type must be 'DUP' or 'DEL'")
    }
    bc <- query_bincov(bincov_con, header, chr, qstart, qend)
    if (is.null(bc)) {
        message("no coverage intervals overlap query region")
        next
    }

    bc_coords <- GRanges(bc[["chr"]], IRanges(bc[["start"]], bc[["end"]]))
    bc[, c("chr", "start", "end") := list(NULL)]
    bc <- bc[, ..samples]
    nonsd_idx <- setdiff(seq_along(bc_coords), queryHits(findOverlaps(bc_coords, sd_regions)))
    gd_idx <- queryHits(findOverlaps(bc_coords, GRanges(chr, IRanges(gdstart, gdend))))

    # normalize
    norms <- bc[, mapply(`/`, .SD, ..cov_medians[names(.SD)], SIMPLIFY = FALSE)]

    carriers <- predict_carriers(norms[intersect(nonsd_idx, gd_idx), ], chr, ploidy, min_shifted_bins, svtype, min_shift)

    if (length(carriers) == 0) {
        message("no potential CNV carriers")
        next
    } else if (length(carriers) == 1 && is.na(carriers)) {
        message("insufficient number of bins overlapping region")
        next
    }

    if (length(samples) - length(carriers) < MAX_BACKGROUND) {
        bg_samples <- samples
    } else {
        bg_samples <- setdiff(samples, carriers)
    }
    bg_samples_to_plot <- sample(bg_samples, min(length(bg_samples), MAX_BACKGROUND))

    if (length(bc_coords) > ROLLING_MEDIAN_WINDOW) {
        norms[, names(.SD) := lapply(.SD, \(x) runmed(x, ROLLING_MEDIAN_WINDOW, na.action = "fail"))]
    }

    mids <- start(bc_coords) + floor(width(bc_coords) / 2)
    bins_to_plot <- spaced_intervals(length(mids))
    size_pretty <- format_size(gdend - gdstart + 1)
    sd_overlaps <- get_sd_overlaps(sd_regions, chr, qstart, qend)
    for (j in seq_along(carriers)) {
        carrier_id <- carriers[[j]]
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
        add_plot_to_store(plots_store, carrier_id, plot_path, cluster)
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

violators_con <- file(if (is.null(violators)) nullfile() else violators)
maphash(plots_store, \(k, v) remove_outliers(k, v, max_calls_per_sample, violators_con))
