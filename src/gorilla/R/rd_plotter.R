# minimum size of region extracted from RD matrix required to use smoothing of
# the RD values
MIN_SMOOTH_REGION <- 10000

# number of bins to use for smoothing the RD values before plotting
SMOOTHING_WINDOW <- 31

# fraction of bins of a region to plot
# lower values speed up plotting, but limit the visible CNV resolution
BIN_PLOT_FRACTION <- 0.05

# minimum number of intervals to plot
MIN_PLOT_INTERVALS <- 100

#' Create a new `rd_plotter` object
#'
#' A `rd_plotter` can create a visualization of read depth.
#'
#' @param x A [`rd_mat`][query()] object.
#' @returns A `rd_plotter` object.
#' @export
#'
#' @examples
#' rd_path <- system.file(
#'     "extdata",
#'     "example.RD.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' medians_path <- system.file(
#'     "extdata",
#'     "example_medianCov.txt",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#'
#' rd_medians <- read_median_coverages(medians_path)
#' rd <- rd_file(rd_path, rd_medians)
#' mat <- query(rd, "chr16", 28743149, 28745149)
#' plotter <- rd_plotter(mat)
rd_plotter <- function(x) {
    new_rd_plotter(x)
}

#' Plot a `rd_plotter` object
#'
#' The `plot()` method for `rd_plotter`.
#'
#' The only additional parameter that is accepted is `col`, which should be a
#' string or a character vector equal to the length of `samples` giving the
#' colors to use for each sample.
#'
#' @param x A [`rd_plotter`][rd_plotter()] object.
#' @param y Ignored.
#' @param samples Samples to plot.
#' @param bg_samples Samples to use as background samples. If `NULL`,
#'   background samples are not included.
#' @param main Title for the plot. Use `NULL` to omit.
#' @param col Color(s) to use for `samples`. Either a string or a character
#'   vector of hexadecimal colors equal to the length of `samples`
#' @param ... Ignored.
#' @name plot.rd_plotter
#' @export
#'
#' @examples
#' rd_path <- system.file(
#'     "extdata",
#'     "example.RD.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' medians_path <- system.file(
#'     "extdata",
#'     "example_medianCov.txt",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#'
#' rd_medians <- read_median_coverages(medians_path)
#'
#' rd <- rd_file(rd_path, rd_medians)
#' mat <- query(rd, "chr16", 28743149, 28745149)
#' plotter <- rd_plotter(mat)
#' plot(plotter, samples = c("gorilla0000", "gorilla0001"))
plot.rd_plotter <- function(
    x,
    y,
    samples,
    bg_samples = NULL,
    main = "",
    col = "#000000",
    ...
) {
    stopifnot(is.character(samples) && length(samples) > 0)
    stopifnot(
        is.null(bg_samples) ||
            (is.character(bg_samples) && length(bg_samples) > 0)
    )
    stopifnot(is_string(main))
    stopifnot(
        is.character(col) &&
            (length(col) == 1 || length(col) == length(samples))
    )

    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(0, 2.5),
        ylab = "RD Ratio",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        las = 2,
        main = main
    )

    draw_bg_samples_rd(x, bg_samples)
    draw_samples_rd(x, samples, col)
}

draw_samples_rd <- function(x, samples, color) {
    if (is.null(samples)) {
        return()
    }

    for (i in seq_along(samples)) {
        graphics::lines(
            x$bin_mids,
            x$mat[x$bins_to_plot, samples[[i]], drop = FALSE],
            col = if (length(color) == 1) color else color[[i]],
            lwd = 2
        )
    }
}

draw_bg_samples_rd <- function(x, bg_samples) {
    if (is.null(bg_samples)) {
        return()
    }

    for (s in bg_samples) {
        graphics::lines(
            x$bin_mids,
            x$mat[x$bins_to_plot, s, drop = FALSE],
            col = "grey",
            lwd = 0.5
        )
    }
}

new_rd_plotter <- function(x) {
    stopifnot(inherits(x, "rd_mat"))

    if (x$region$end - x$region$start + 1 >= MIN_SMOOTH_REGION) {
        x <- smooth_rd(x, SMOOTHING_WINDOW)
    }
    bins_to_plot <- spaced_intervals(length(x$ranges))
    bin_mids <- GenomicRanges::start(x$ranges[bins_to_plot]) +
        (GenomicRanges::width(x$ranges[bins_to_plot]) %/% 2)

    structure(
        list(
            bins_to_plot = bins_to_plot,
            bin_mids = bin_mids,
            mat = x$mat,
            region = x$region
        ),
        class = "rd_plotter"
    )
}

# Generate a vector of best-effort equally spaced integers from 1 to n.
spaced_intervals <- function(n) {
    nintervals <- ceiling(n * BIN_PLOT_FRACTION)
    if (nintervals < MIN_PLOT_INTERVALS) {
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
