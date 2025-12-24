#' Forage for genomic disorders.
#'
#' A convenience function to run the genomic disorder calling workflow on a set
#' of inputs.
#'
#' @param gds_path Path to the genomic disorders table file. See
#'   [`read_gdtable`][read_gdtable()]
#' @param sds_path Path to a BEDX file of hg38 segmental duplication regions.
#' @param bincov_path Path to the binned coverage matrix.
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
#' @param min_shift Minimum amount by which the read-depth ratio of a sample
#'   must be shifted from expectation to be considered a GD carrier.
#' @param pad Proportion of the GD region that should be included on either
#'   side when plotting.
#' @param max_calls_per_sample Maximum number of GD calls a sample can have
#'   before it is marked as an outlier and all of its calls are treated as false.
#' @param outliers_log_path Path to a file where the IDs of any outlier samples
#'   will be written. If `NULL`, outlier samples are not recorded.
#' @export
forage <- function(
    gds_path,
    sds_path,
    bincov_path,
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
    sds <- read_segdups(sds_path)
    medians <- read_median_coverages(medians_path)
    bf <- bincov_file(bincov_path, medians)
    target_samples <- readLines(samples_path)
    if (length(target_samples) == 0) {
        stop("at least one sample to check for GDs must be given")
    }

    if (!all(target_samples %in% bf$header)) {
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
        bincov_mat <- query(bf, gd$chr, gd$qstart, gd$qend)
        caller <- new_gdcaller(gd, bincov_mat, sds, sex_ploidy)
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
    plotter <- as_gdplotter(caller)
    for (i in seq_along(caller$carriers)) {
        plot_name <- sprintf(
            "%s_%d-%d_%s_%s_%s.jpg",
            caller$gd$chr,
            caller$gd$start_GRCh38,
            caller$gd$end_GRCh38,
            caller$gd$GD_ID,
            caller$gd$svtype,
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
