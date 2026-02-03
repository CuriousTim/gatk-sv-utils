# the fraction of the region to use as the window size for non-NAHR GDs
NON_NAHR_WINDOW_PROP <- 0.01
# the minimum number of bins to use as the window size for non-NAHR GDs
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

new_gdcaller <- function(gd, rd, sex_ploidy) {
    structure(
        list(
            gd = gd,
            rd = rd,
            segdups = segdups_gr,
            sex_ploidy = sex_ploidy
        ),
        class = "gdcaller"
    )
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
        x$rd$ranges,
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
        gd_rd_mat <- x$rd$mat[gd_idx, x$targets, drop = FALSE]
        gd_rd_ranges <- x$rd$ranges[gd_idx]

        sd_idx <- S4Vectors::queryHits(GenomicRanges::findOverlaps(
            gd_rd_ranges,
            x$segdups
        ))
        gd_rd_mat[sd_idx, ] <- NA_real_

        sex <- NULL
        sexes <- x$sex_ploidy[colnames(gd_rd_mat), sex]
        expected_rd <- get_expected_rd(
            x$gd$chr,
            x$gd$start_GRCh38,
            x$gd$end_GRCh38,
            sexes
        )
        rd_comparator <- make_rd_comparator(x$gd$svtype, x$min_shift)

        carriers <- caller(
            x,
            gd_rd_ranges,
            gd_rd_mat,
            sd_idx,
            expected_rd,
            rd_comparator
        )
    }

    x$carriers <- carriers

    x
}

nahr_caller <- function(
    x,
    gd_rd_ranges,
    gd_rd_mat,
    sd_idx,
    expected_rd,
    rd_comparator
) {
    if (length(gd_rd_ranges) == length(sd_idx)) {
        character()
    } else {
        m <- matrixStats::colMedians(gd_rd_mat, na.rm = TRUE)
        names(m[rd_comparator(expected_rd, m)])
    }
}

non_nahr_caller <- function(
    x,
    gd_rd_ranges,
    gd_rd_mat,
    sd_idx,
    expected_rd,
    rd_comparator
) {
    if (x$gd$terminal == "no") {
        non_terminal_non_nahr_caller(
            x,
            gd_rd_ranges,
            gd_rd_mat,
            sd_idx,
            expected_rd,
            rd_comparator
        )
    } else {
        terminal_non_nahr_caller(
            x,
            gd_rd_ranges,
            gd_rd_mat,
            expected_rd,
            rd_comparator
        )
    }
}

terminal_non_nahr_caller <- function(
    x,
    gd_rd_ranges,
    gd_rd_mat,
    sd_idx,
    expected_rd,
    rd_comparator
) {
    carriers <- character()
    window_size <- non_nahr_window_size(x)

    if (length(gd_rd_ranges) >= window_size) {
        if (x$gd$terminal == "p") {
            bins <- seq_len(window_size)
        } else {
            bins <- seq.int(
                length(gd_rd_ranges) - window_size + 1,
                length(gd_rd_ranges)
            )
        }

        medians <- matrixStats::colMedians(gd_rd_mat, rows = bins, na.rm = TRUE)
        if (all(!is.nan(medians))) {
            carriers <- names(medians[rd_comparator(expected_rd, medians)])
        }
    }

    carriers
}

non_terminal_non_nahr_caller <- function(
    x,
    gd_rd_ranges,
    gd_rd_mat,
    sd_idx,
    expected_rd,
    rd_comparator
) {
    carriers <- character()
    window_size <- non_nahr_window_size(x)

    tmp <- rep.int(0, length(gd_rd_ranges))
    tmp[sd_idx] <- 1L
    sums <- data.table::frollsum(
        tmp,
        window_size,
        align = "center",
        hasNA = FALSE
    )

    # find windows with at least one bin that is not overlapping a segdup
    callable_windows <- which(!is.na(sums) & sums < window_size)
    if (length(callable_windows) > 0) {
        m <- apply(gd_rd_mat, 2, \(x) stats::runmed(x, window_size))
        m <- m[callable_windows, , drop = FALSE]
        pass <- vapply(
            seq_len(ncol(m)),
            \(i) any(rd_comparator(expected_rd[[i]], m[, i])),
            logical(1)
        )
        carriers <- colnames(m)[pass]
    }

    carriers
}

# Return a function that accepts vectors of expected read-depth and
# actual read-depth, and reports which actual read-depth values are higher or
# lower than expected, given the expected SV type and minimum shift.
make_rd_comparator <- function(svtype, min_shift) {
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

# Compute the expected read-depth for some samples of the given sexes,
# given a genomic region. `sex_ploidy` can be
# (NA = unknown, 1 = male, 2 = female).
get_expected_rd <- function(chr, start, end, sex_ploidy) {
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

    # window size must be odd for rolling median
    if (w %% 2 == 0) {
        w <- w + 1
    }

    w
}
