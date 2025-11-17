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
# the window to use for smoothing the binned coverage values before plotting
ROLLING_MEDIAN_WINDOW <- 31
# maximum number of background samples to plot
MAX_BACKGROUND <- 200L
# the fraction of bins to plot, low mostly to reduce computation
BIN_PLOT_FRACTION <- 0.05
# the fraction of the region to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_PROP <- 0.01
# the assumed width of each bin in the bincov matrix
BIN_WIDTH <- 100
# minimum number of bins (assuming 100 bp bins) overlapping a GD region required
# to make a CNV call
MIN_GD_OVP_BINS <- 50

assert_coords_in_range <- function(x) {
    if (any(x <= 0L | x > TABIX_MAX_SEQLEN)) {
        stop(sprintf(
            "only chromosome coordinates between 0 and %d are permitted",
            TABIX_MAX_SEQLEN
        ))
    }
}

assert_hg38_chr <- function(x) {
    if (any(!grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", x))) {
        stop("only chromosomes chr1-22,X,Y are permitted")
    }
}

assert_positive_ranges <- function(start, end) {
    if (any(end - start + 1L <= 0)) {
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

    setNames(medians, ids)
}

read_samples_file <- function(path) {
    tmp <- readLines(path)
    if (length(tmp) == 0) {
        stop("at least one sample to check for CNVs must be given")
    }

    tmp
}

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

read_gd <- function(path) {
    tmp <- fread(path, header = TRUE, sep = "\t", strip.white = FALSE,
                 na.strings = "",
                 colClasses = c("character", "integer", "integer", "character", "character", "character", "character"))
    req_header <- c("chr", "start_GRCh38", "end_GRCh38", "GD_ID", "svtype", "NAHR", "cluster")
    disp_header <- paste0(req_header, collapse = ", ")
    if (!identical(colnames(tmp), req_header)) {
        stop(sprintf("genomic disorder regions file must have header: '%s'",
                     disp_header))
    }

    assert_hg38_chr(tmp$chr)
    assert_coords_in_range(tmp$start_GRCh3)
    assert_coords_in_range(tmp$end_GRCh38)
    assert_positive_ranges(tmp$start_GRCh38, tmp$end_GRCh38)
    assert_valid_svtypes(tmp$svtype)
    assert_valid_nahr(tmp$NAHR)

    tmp[, NAHR := NAHR == "yes"]

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
    ovp <- findOverlaps(GRanges(chr, IRanges(start, end)), sds)
    sds[subjectHits(ovp)]
}

make_rd_shift_checker <- function(bc, chr, sex_ploidy, svtype, min_shift) {
    if (chr == "chrX" || chr == "chrY") {
        ploidy <- sex_ploidy[colnames(bc), sex]
        if (chr == "chrY") {
            ploidy[ploidy == 2] <- 0L
        }
        expected_rdr <- ploidy / 2
    } else {
        expected_rdr <- 1
    }

    if (svtype == "DUP") {
        op <- function(x) { !is.nan(x) & !is.na(x) & !is.na(expected_rdr) & x >= expected_rdr + min_shift }
    } else {
        op <- function(x) { !is.nan(x) & !is.na(x) & !is.na(expected_rdr) & x <= expected_rdr - min_shift }
    }

    op
}

predict_carriers <- function(bc, sd_idx, gd_idx, window, filter) {
    bc <- copy(bc)
    bc[sd_idx, names(.SD) := list(NA_real_)]
    gd_bc <- bc[gd_idx, ]

    if (!is.null(window)) {
        keep <- seq.int(ceiling(window / 2), nrow(gd_bc) - floor(window / 2))
        roll <- lapply(gd_bc, \(x) runmed(x, window)[keep])
        pass <- vapply(roll, \(x) any(filter(x), na.rm = TRUE), logical(1))
        pass <- pass[pass == TRUE]
    } else {
        m <- vapply(gd_bc, \(x) median(x, na.rm = TRUE), double(1))
        pass <- m[filter(m)]
    }

    names(pass)
}

validate_args <- function(x) {
    if (!is.finite(x$min_shift) || x$min_shift < 0) {
        stop("minimum shift must be a non-negative number")
    }

    if (!is.finite(x$pad) || x$pad < 0) {
        stop("padding must be a non-negative number")
    }

    if (x$max_calls_per_sample < 0) {
        stop("max calls per sample must be a non-negative number")
    }

    x
}

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

gd_regions <- read_gd(argv$gd_regions)
sd_regions <- read_sd(argv$sd_regions)
bincov_con <- TabixFile(argv$bincov)
cov_medians <- read_medians_file(argv$medians)
target_samples <- read_samples_file(argv$samples)
sex_ploidy <- read_sex_ploidy(argv$sex_ploidy)

dir.create(argv$outdir)

header <- read_bincov_header(bincov_con)
if (!all(target_samples %in% header)) {
    stop("all requested samples must be in the bincov matrix")
}

if (!all(target_samples %in% names(cov_medians))) {
    stop("all requested samples must be in the median coverage file")
}

# sex chromosome ploidy is needed to call CNVs on chrX and chrY
sex_ploidy <- sex_ploidy[target_samples, nomatch = NA, on = "sample_id"]
setkey(sex_ploidy, sample_id)

plots_store <- hashtab()

pad_size <- gd_regions[, ceiling((end_GRCh38 - start_GRCh38 + 1L) * argv$pad)]
expanded_gd <- gd_regions[, list(chr = chr, start = pmax(1L, start_GRCh38 - pad_size), end = pmin(TABIX_MAX_SEQLEN, end_GRCh38 + pad_size))]
for (i in seq_len(nrow(gd_regions))) {
    chr <- gd_regions[i, chr]
    gdstart <- gd_regions[i, start_GRCh38]
    gdend <- gd_regions[i, end_GRCh38]
    gdid <- gd_regions[i, GD_ID]
    svtype <- gd_regions[i, svtype]
    nahr <- gd_regions[i, NAHR]
    gdcluster <- gd_regions[i, cluster]
    qstart <- expanded_gd[i, start]
    qend <- expanded_gd[i, end]
    message(sprintf("checking for genomic disorder at %s:%d-%d (%s)",
                    chr, gdstart, gdend, svtype))

    bc <- query_bincov(bincov_con, header, chr, qstart, qend)
    if (is.null(bc)) {
        message("no coverage intervals overlap query region")
        next
    }

    bc_coords <- GRanges(bc[["chr"]], IRanges(bc[["start"]], bc[["end"]]))
    sd_idx <- queryHits(findOverlaps(bc_coords, sd_regions))
    gd_idx <- queryHits(findOverlaps(bc_coords, GRanges(chr, IRanges(gdstart, gdend))))
    if (length(gd_idx) < MIN_GD_OVP_BINS) {
        message("insufficient coverage intervals overlap query region")
        next
    }
    bc[, c("chr", "start", "end") := list(NULL)]
    bc <- bc[, ..target_samples]

    norms <- bc[, mapply(`/`, .SD, ..cov_medians[names(.SD)], SIMPLIFY = FALSE)]
    if (nahr) {
        window <- NULL
    } else {
        window <- round((gdend - gdstart + 1) * NON_NAHR_WINDOW_PROP / BIN_WIDTH)
        if (window %% 2 == 0) {
            window <- window + 1
        }
    }
    filter_func <- make_rd_shift_checker(norms, chr, sex_ploidy, svtype, argv$min_shift)
    carriers <- predict_carriers(norms, sd_idx, gd_idx, window, filter_func)

    if (length(carriers) == 0) {
        message("no potential CNV carriers")
        next
    } else if (length(carriers) == 1 && is.na(carriers)) {
        message("insufficient number of bins overlapping region")
        next
    }

    if (length(target_samples) - length(carriers) < MAX_BACKGROUND) {
        bg_samples <- target_samples
    } else {
        bg_samples <- setdiff(target_samples, carriers)
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
        plot_path <- file.path(argv$outdir, plot_name)
        add_plot_to_store(plots_store, carrier_id, plot_path, gdcluster)
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

violators_con <- file(if (is.null(argv$violators)) nullfile() else argv$violators)
maphash(plots_store, \(k, v) remove_outliers(k, v, argv$max_calls_per_sample, violators_con))
