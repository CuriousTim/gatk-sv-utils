# Predict CNV carriers at GD regions.
#
# The general strategy for predicting carriers is to check the read depth
# ratio of samples over a genomic region and pick samples that have elevated
# or depressed values. Whether the read depth ratio is elevated or depressed
# depends on the expected value and the minimum shift. The expected value
# depends on the genomic region and sex of the sample. For autosomes and PAR
# regions, the expected ratio is 2. For males, it is 1 on chrX and chrY. For
# females, it is 2 on chrX and 0 on chrY. The minimum shift is a input
# parameter to the algorithm which determines the required deviation of the
# read depth ratio from the expectation.
#
# For NAHR-mediated GDs, the entire region is tested.
#
# For non-NAHR-mediated GDs, the test depends on the type.
# For p-arm GDs, the first window (smaller genomic position) in the region is
# tested.
# For q-arm GDs, the last window (larger genomic position) in the region is
# tested.
# For non-terminal GDs, a rolling window across the region is tested.

# the fraction of the region to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_PROP <- 0.01
# the minimum number of windows to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_MIN_BINS <- 101
# the fraction of terminal non-NAHR regions that should be used as the window
# size
TERMINALS_WINDOW_PROP <- 0.2
# the assumed width of each bin in the bincov matrix
BIN_WIDTH <- 100
# minimum number of bins (assuming 100 bp bins) overlapping a NAHR GD region
# required to make a CNV call
NAHR_WINDOW_MIN_BINS <- 50

# HG38 pseudoautosomal regions
CHRX_PAR1_START <- 10001
CHRX_PAR1_END <- 2781479
CHRX_PAR2_START <- 155701383
CHRX_PAR2_END <- 156030895
CHRY_PAR1_START <- 10001
CHRY_PAR1_END <- 2781479
CHRY_PAR2_START <- 56887903
CHRY_PAR2_END <- 57217415

new_gdcaller <- function(gd, bincov_mat, segdups, sex_ploidy) {
    structure(
        list(
            gd = gd,
            bincov_mat = bincov_mat,
            segdups = segdups,
            sex_ploidy = sex_ploidy
        ),
        class = "gdcaller"
    )
}

as_gdplotter <- function(x) {
    new_gdplotter(x$gd, x$bincov_mat, x$segdups, x$carriers)
}

call_gds <- function(x, targets, min_shift) {
    UseMethod("call_gds")
}

#' @export
call_gds.gdcaller <- function(x, targets, min_shift) {
    x$targets <- targets
    x$min_shift <- min_shift

    run_gd_prediction(x)
}

run_gd_prediction <- function(x) {
    carriers <- character()
    gd_idx <- S4Vectors::queryHits(GenomicRanges::findOverlaps(
        x$bincov_mat$ranges,
        GenomicRanges::GRanges(
            x$gd$chr,
            IRanges::IRanges(x$gd$start_GRCh38, x$gd$end_GRCh38)
        )
    ))
    if (x$gd$NAHR) {
        min_bins <- NAHR_WINDOW_MIN_BINS
        caller <- nahr_caller
    } else {
        min_bins <- NON_NAHR_WINDOW_MIN_BINS
        caller <- non_nahr_caller
    }

    if (length(gd_idx) >= min_bins) {
        gd_bincov <- x$bincov_mat[gd_idx, x$targets]
        sd_idx <- S4Vectors::queryHits(GenomicRanges::findOverlaps(
            gd_bincov$ranges,
            x$segdups
        ))
        gd_bincov$rd[sd_idx, ] <- NA_real_
        sex <- NULL
        sexes <- x$sex_ploidy[colnames(gd_bincov$rd), sex]
        expected_rdr <- get_expected_rdr(
            x$gd$chr,
            x$gd$start_GRCh38,
            x$gd$end_GRCh38,
            sexes
        )
        rdr_comparator <- make_rdr_comparator(x$gd$svtype, x$min_shift)

        carriers <- caller(x, gd_bincov, sd_idx, expected_rdr, rdr_comparator)
    }

    x$carriers <- carriers

    x
}

nahr_caller <- function(x, gd_bincov, sd_idx, expected_rdr, rdr_comparator) {
    if (length(gd_bincov$ranges) == length(sd_idx)) {
        character()
    } else {
        m <- matrixStats::colMedians(gd_bincov$rd, na.rm = TRUE)
        names(m[rdr_comparator(expected_rdr, m)])
    }
}

non_nahr_caller <- function(
    x,
    gd_bincov,
    sd_idx,
    expected_rdr,
    rdr_comparator
) {
    if (x$gd$terminal == "no") {
        non_terminal_non_nahr_caller(
            x,
            gd_bincov,
            sd_idx,
            expected_rdr,
            rdr_comparator
        )
    } else {
        terminal_non_nahr_caller(x, gd_bincov, expected_rdr, rdr_comparator)
    }
}

terminal_non_nahr_caller <- function(
    x,
    gd_bincov,
    expected_rdr,
    rdr_comparator
) {
    carriers <- character()
    window_size <- non_nahr_window_size(x)

    if (x$gd$terminal == "p") {
        bins <- seq_len(window_size)
    } else {
        bins <- seq.int(nrow(gd_bincov$rd) - window_size, nrow(gd_bincov$rd))
    }

    if (nrow(gd_bincov$rd) >= length(bins)) {
        m <- matrixStats::colMedians(gd_bincov$rd, rows = bins, na.rm = TRUE)
        if (all(!is.na(m))) {
            carriers <- names(m[rdr_comparator(expected_rdr, m)])
        }
    }

    carriers
}

non_terminal_non_nahr_caller <- function(
    x,
    gd_bincov,
    sd_idx,
    expected_rdr,
    rdr_comparator
) {
    carriers <- character()
    window_size <- non_nahr_window_size(x)
    tmp <- rep.int(0, length(gd_bincov$ranges))
    tmp[sd_idx] <- 1L
    sums <- data.table::frollsum(
        tmp,
        window_size,
        align = "center",
        hasNA = FALSE
    )
    callable_windows <- which(!is.na(sums) & sums < window_size)
    if (length(callable_windows) > 0) {
        m <- apply(gd_bincov$rd, 2, \(x) stats::runmed(x, window_size))
        m <- m[callable_windows, , drop = FALSE]
        pass <- vapply(
            seq_len(ncol(m)),
            \(i) any(rdr_comparator(expected_rdr[[i]], m[, i])),
            logical(1)
        )
        carriers <- colnames(m)[pass]
    }

    carriers
}

# Return a function that accepts an vector of expected read-depth ratios and
# actual read-depth ratios, and reports which actual ratios are higher or lower
# than expected, given the expected SV type and minimum shift.
make_rdr_comparator <- function(svtype, min_shift) {
    op <- if (svtype == "DUP") {
        function(expected, actual) {
            actual >= expected + min_shift
        }
    } else {
        function(expected, actual) {
            actual <= expected - min_shift
        }
    }

    function(expected, actual) {
        !is.nan(actual) &
            !is.na(actual) &
            !is.na(expected) &
            op(expected, actual)
    }
}

# Compute the expected read depth ratio for some samples of the given sexes,
# given a genomic region.
# `sex_ploidy` can be (NA = unknown, 1 = male, 2 = female).
get_expected_rdr <- function(chr, start, end, sex_ploidy) {
    tmp <- if (chr == "chrX") {
        if (ovp_chrx_par(start, end)) {
            sex_ploidy[sex_ploidy == 1L] <- 2L
        }
        sex_ploidy
    } else if (chr == "chrY") {
        sex_ploidy[sex_ploidy == 2L] <- 0L
        if (ovp_chry_par(start, end)) {
            sex_ploidy[sex_ploidy == 1L] <- 2L
        }
        sex_ploidy
    } else {
        rep_len(2, length(sex_ploidy))
    }

    tmp / 2
}

coordinates_ovp <- function(x1, y1, x2, y2) {
    x1 <= y2 & y1 >= x2
}

# Does a region on chrY overlap PAR?
ovp_chry_par <- function(start, end) {
    coordinates_ovp(start, end, CHRY_PAR1_START, CHRY_PAR1_END) |
        coordinates_ovp(start, end, CHRY_PAR2_START, CHRY_PAR2_END)
}

# Does a region on chrX overlap PAR?
ovp_chrx_par <- function(start, end) {
    coordinates_ovp(start, end, CHRX_PAR1_START, CHRX_PAR1_END) |
        coordinates_ovp(start, end, CHRX_PAR2_START, CHRX_PAR2_END)
}

# Compute the window size, in bins, for non-NAHR GD regions.
non_nahr_window_size <- function(x) {
    frac <- if (x$gd$terminal == "p" || x$gd$terminal == "q") {
        TERMINALS_WINDOW_PROP
    } else {
        NON_NAHR_WINDOW_PROP
    }
    w <- max(
        NON_NAHR_WINDOW_MIN_BINS,
        round((x$gd$end_GRCh38 - x$gd$start_GRCh38 + 1) * frac / BIN_WIDTH)
    )
    if (w %% 2 == 0) {
        w <- w + 1
    }

    w
}
