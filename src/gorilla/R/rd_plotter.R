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

#' Create a new `rd_plotter` object.
#'
#' A `rd_plotter` can create a visualization of read depth.
#'
#' @param x A [`rd_mat`][query()] object.
#' @returns A `rd_plotter` object.
#' @export
rd_plotter <- function(x) {
    new_rd_plotter(x)
}

#' Plot a `rd_plotter` object.
#'
#' @param x A [`rd_plotter`] object.
#' @param y Kept for compatability with generic, but ignored here.
#' @param samples Samples to plot.
#' @param bg_samples Samples to use as background samples. If `NULL`,
#'   background samples are not included.
#' @param ... Other plotting parameters. Currently the only accepted parameter
#'   is `col` which should either be a character vector of length 1 or length
#'   equal to `samples` giving the colors to use for each sample.
#' @name plot.rd_plotter
#' @export
plot.rd_plotter <- function(x, y, samples, bg_samples = NULL, ...) {
    stopifnot(is.character(samples) && length(samples) > 0)
    stopifnot(is.null(bg_samples) || (is.character(bg_samples) && length(bg_samples) > 0))

    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(0, 2.5),
        ylab = "RD Ratio",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        las = 2
    )

    col <- tryCatch(
        color_from_dots(length(samples), ...),
        value_error = function(e) stop("length of `col` must be 1 or equal to length of `samples`")
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
            col = color[[i]],
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
        x <- smooth(x, SMOOTHING_WINDOW)
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
