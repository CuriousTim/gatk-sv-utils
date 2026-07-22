#' Forage for genomic disorders
#'
#' Predict genomic disorder CNVs using read depth evidence from GATK-SV.
#'
#' @details
#' # Algorithm Overview
#'
#' Given a set of genomic disorder regions, inspect the binned read-depth
#' over each region for each sample and determine if a sample potentially
#' carries a CNV by comparing its normalized read-depth (NRD) to its expected
#' value.  Samples with sufficiently different NRD are considered to be
#' carriers. The threshold for "sufficiently different" is given by
#' `min_shift`. For XX samples, the expected NRD is 1 for chr1-22 and chrX and
#' 0 for chrY. For XY, the expected NRD is 1 for chr1-22 and the
#' pseudoatuosomal regions of chrX and chrY and 0.5 for the rest of chrX and
#' chrY. For all other samples, the expected NRD is 1 for chr1-22, but chrX and
#' chrY are ignored. There are slightly different tests depending on the type
#' of the genomic disorder, but any base that overlaps a segmental duplication
#' is excluded from consideration.
#'
#' ## DEL vs DUP
#'
#' For duplications, the sample NRD must be greater than its expected NRD,
#' while for deletions, the sample NRD must be less than its expected NRD.
#'
#' ## NAHR-mediated
#'
#' The median of the NRD over the entire region is compared to
#' to the expected value for each sample. There must be at least 50 read-depth
#' bins overlapping the region to make a call.
#'
#' ## Terminal non-NAHR mediated
#'
#' A window of 20% of the region (minimum 101 bins) anchored to one end of
#' the region (upstream for p-arm, downstream for q-arm) is used to compute the
#' median NRD and compared to the expected NRD.
#'
#' ## Non-terminal non-NAHR mediated
#'
#' A window of 1% of the region (minimum 101 bins) is used to compute a rolling
#' median of NRD over the region and the value of each window is compared to the
#' expected NRD. If a sample is sufficiently different from expectation at any
#' window, it is considered a carrier.
#'
#' # Output
#'
#' The algorithm will create a visualization of the read-depth evidence for
#' each call it makes.
#'
#' @param gds_path Path to the genomic disorders table file. See
#'   [read_gdtable()]
#' @param rd_path Path to the binned coverage matrix.
#' @param medians_path Path to the median coverages file.
#' @param samples_path Path to a list of samples, one-per-line, to process.
#' @param sex_ploidy_path Path to a ploidy table. The format is a tab-delimited
#'   file with the following columns:
#'   1. sample ID
#'   2. chrX ploidy with 0 = unknown, 1 = male, 2 = female
#'   The file must not have column headers.
#' @param outdir_path Path the output directory. It will be created if it
#'   doesn't exist and any files that match the name of a plot created by the
#'   workflow will be silently overwritten.
#' @param min_shift Minimum amount by which the normalized read-depth of a
#'   sample must be shifted from expectation to be considered a carrier.
#' @param pad Proportion of the GD region that should be included on either
#'   side when plotting.
#' @param max_calls_per_sample Maximum number of GD calls a sample can have
#'   before it is marked as an outlier and all of its calls are treated as false.
#' @param outliers_log_path Path to a file where the IDs of any outlier samples
#'   will be written. If `NULL`, outlier samples are not recorded.
#' @export
forage <- function(
    gds_path,
    rd_path,
    medians_path,
    samples_path,
    sex_ploidy_path,
    outdir_path,
    min_shift = 0.3,
    pad = 0.5,
    max_calls_per_sample = 10,
    outliers_log_path = NULL
) {
    gds <- read_gdtable(gds_path)
    medians <- read_median_coverages(medians_path)
    rd_handle <- rd_file(rd_path, medians)
    target_samples <- readLines(samples_path)
    if (length(target_samples) == 0) {
        stop("at least one sample to check for GDs must be given")
    }

    if (!all(target_samples %in% rd_handle$header)) {
        stop("all requested samples must be in the bincov matrix")
    }

    sex_ploidy <- read_sex_ploidy(sex_ploidy_path)
    dir.create(outdir_path)

    active_samples <- target_samples
    expand_gdtable(gds, pad)

    calls_record <- utils::hashtab()

    for (i in seq_len(nrow(gds))) {
        if (length(active_samples) == 0) {
            break
        }

        gd <- gds[i, ]
        rd <- query(rd_handle, gd$chr, gd$qstart, gd$qend)
        caller <- new_gdcaller(gd, rd, sex_ploidy)
        caller <- call_gds(caller, active_samples, min_shift)

        if (length(caller$carriers) == 0) {
            next
        }

        make_gd_plots(caller, outdir_path, calls_record)
        outliers <- caller$carriers[is_outlier(
            calls_record,
            caller$carriers,
            max_calls_per_sample
        )]
        active_samples <- setdiff(active_samples, outliers)
    }

    cleanup_outliers(calls_record, max_calls_per_sample, outliers_log_path)
}

# Helper function to generate plots for all the GD carriers after calling.
make_gd_plots <- function(caller, outdir, calls_record) {
    plotter <- new_gd_plotter(caller)
    for (i in seq_along(caller$carriers)) {
        plot_name <- sprintf(
            "%s~~%s.jpg",
            caller$gd$GD_ID,
            caller$carriers[[i]]
        )
        plot_path <- file.path(outdir, plot_name)
        grDevices::jpeg(plot_path, res = 100, width = 960, height = 540)
        plot(plotter, carrier = i)
        grDevices::dev.off()
        add_call_to_record(
            calls_record,
            caller$carriers[[i]],
            plot_path,
            caller$gd$cluster
        )
    }
}

# Remove plots from samples with excessive number of calls.
cleanup_outliers <- function(calls_record, threshold, outliers_log_path) {
    con <- file(
        if (is.null(outliers_log_path)) nullfile() else outliers_log_path
    )
    f <- function(k, v) {
        if (v[["count"]] > threshold) {
            writeLines(k, con)
            file.remove(v[["paths"]])
        }
    }
    utils::maphash(calls_record, f)
}

# Add a call/plot to the calls record.
# The calls record records the number of calls each sample has and the paths to
# the plots for those calls. For GDs that are part of a cluster, all calls in
# the cluster only count once.
add_call_to_record <- function(calls_record, sid, plot_path, cluster) {
    val <- utils::gethash(calls_record, sid)
    if (is.null(val)) {
        val <- list(
            clusters = utils::hashtab(),
            paths = character(),
            count = 0L
        )
    }
    val[["paths"]] <- append(val[["paths"]], plot_path)

    if (
        !is.na(cluster) && is.null(utils::gethash(val[["clusters"]], cluster))
    ) {
        utils::sethash(val[["clusters"]], cluster, NULL)
        val[["count"]] <- val[["count"]] + 1L
    }

    if (is.na(cluster)) {
        val[["count"]] <- val[["count"]] + 1L
    }

    utils::sethash(calls_record, sid, val)
}

is_outlier <- function(calls_record, sid, threshold) {
    f <- function(s) {
        val <- utils::gethash(calls_record, s)
        if (is.null(val)) {
            FALSE
        } else {
            val[["count"]] > threshold
        }
    }

    vapply(sid, f, logical(1))
}
