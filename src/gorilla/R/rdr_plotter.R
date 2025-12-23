# Functions for plotting binned coverage matrices.

# number of bins to use for smoothing the binned coverage values before
# plotting
SMOOTHING_WINDOW <- 31
# fraction of bins to plot
# lower values speed up plotting, but limit the visible CNV resolution
BIN_PLOT_FRACTION <- 0.05
# maximum number of background samples to plot
MAX_BACKGROUND <- 200L

#' Create a new `rdr_plotter` object.
#'
#' An `rdr_plotter` object represents a region of a binned coverage matrix to
#' visualize.
#'
#' @param x A `bincov_mat` object.
#' @param segdups A [`GRanges`][GenomicRanges::GRanges] object with ranges for
#'   segmental duplications.
#' @returns A `rdr_plotter` object.
#' @export
rdr_plotter <- function(x, segdups) {
    if (!inherits(x, "bincov_mat")) {
        stop("`x` must be a `bincov_mat` object")
    }

    if (!methods::is(segdups, "GRanges")) {
        stop("`segdups` must be a `GRanges` object")
    }

    validate_rdr_plotter(new_rdr_plotter(x, segdups))
}

new_rdr_plotter <- function(x, segdups) {
    mat <- smooth(x, SMOOTHING_WINDOW)
    bins_to_plot <- spaced_intervals(length(mat$ranges))
    bin_mids <- GenomicRanges::start(mat$ranges[bins_to_plot]) +
        (GenomicRanges::width(mat$ranges[bins_to_plot]) %/% 2)
    bg_samples <- sample(colnames(mat$rd), min(MAX_BACKGROUND, ncol(mat$rd)))
    span <- utils::head(mat$ranges, 1)
    GenomicRanges::end(span) <- utils::tail(GenomicRanges::end(mat$ranges), 1)
    ovp_sd <- segdups[S4Vectors::queryHits(GenomicRanges::findOverlaps(
        segdups,
        span
    ))]

    structure(
        list(
            bincov_mat = mat,
            bins_to_plot = bins_to_plot,
            bin_mids = bin_mids,
            bg_samples = bg_samples,
            segdups = ovp_sd
        ),
        class = "rdr_plotter"
    )
}

validate_rdr_plotter <- function(x) {
    bincov_contigs <- unique(methods::as(
        GenomicRanges::seqnames(x$segdups),
        "character"
    ))
    segdup_contigs <- unique(methods::as(
        GenomicRanges::seqnames(x$bincov_mat$ranges),
        "character"
    ))

    if (!all(segdup_contigs %in% bincov_contigs)) {
        stop("all segmental duplication contigs must be in the bincov matrix")
    }

    x
}

#' Plot an `rdr_plotter`.
#'
#' @param x An `rdr_plotter` object.
#' @param y Kept for compatibility with the `plot` generic. Ignored here.
#' @param main An optional string to use as the title of the plot.
#' @param carrier Which carriers to plot, if `x` has carriers. If 0, no
#'   carriers are plotting. If `NULL`, all carriers are plotted. If an integer
#'   or character vector, it will be used as an index into the carriers vector.
#' @export
plot.rdr_plotter <- function(x, y = NULL, ...) {
    dots <- list(...)
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(3.1, 4.1, 4.1, 2.1))
    rdr_plotter_plot_base(x, dots[["main"]])
    rdr_plotter_plot_bg(x)
    rdr_plotter_plot_carriers(x, dots[["carrier"]])
    rdr_plotter_plot_segdups(x)
    rdr_plotter_plot_pad(x)
    rdr_plotter_plot_axis(x)
    rdr_plotter_plot_border(x)
    graphics::par(old_par)
}

#' Create a new `gdplotter` object.
#'
#' A `gdplotter` is an object that can plot a genomic disorder (GD) region.
#'
#' @param gd A list-like object that represents a GD region.
#'   It must contain the following named components:
#'   * chr: Contig of the GD region
#'   * start_GRCh38: start of the GD region
#'   * end_GRCh38: end of the GD region
#'   * svtype: SV type of the GD (either 'DEL' or 'DUP')
#' @param bincov_mat A `bincov_mat` object over the GD region.
#' @param segdups A [`GRanges`][GenomicRanges::GRanges()] object of hg38
#'   segmental duplications.
#' @param carriers A character vector of sample IDs of the GD carriers.
#' @returns A `gdplotter` object.
#' @export
gdplotter <- function(gd, bincov_mat, segdups, carriers) {
    new_gdplotter(gd, bincov_mat, segdups, carriers)
}

new_gdplotter <- function(gd, bincov_mat, segdups, carriers) {
    carrier_col <- rep(
        if (gd$svtype == "DEL") "red" else "blue",
        length(carriers)
    )
    plotter <- new_rdr_plotter(bincov_mat, segdups) |>
        rdr_plotter_set_pad(gd$start_GRCh38, gd$end_GRCh38) |>
        rdr_plotter_set_carriers(stats::setNames(carrier_col, carriers))
    size_pretty <- format_size(gd$end_GRCh38 - gd$start_GRCh38 + 1)

    structure(
        list(
            gd = gd,
            carriers = carriers,
            plotter = plotter,
            size_pretty = size_pretty
        ),
        class = "gdplotter"
    )
}

#' @export
plot.gdplotter <- function(x, y = NULL, carrier = 1, ...) {
    carrier_pretty <- format_sample_id(x$carriers[[carrier]])
    main <- sprintf(
        "%s (hg38)\n%s (%s)",
        x$gd$GD_ID,
        carrier_pretty,
        x$size_pretty
    )

    plot(x$plotter, main = main, carrier = carrier)
}

#' Set the carriers for an `rdr_plotter`.
#'
#' @param x An `rdr_plotter` object.
#' @param carriers A named character vector of which the values are colors to
#'   use for plotting the carriers and the names are the IDs of the carriers.
#' @returns An `rdr_plotter` object.
#' @export
rdr_plotter_set_carriers <- function(x, carriers) {
    x$carriers <- carriers

    x
}

#' Set the padding for an `rdr_plotter`.
#'
#' The padding represents the genomic region outside the region of interest.
#' This function does not increase the size of the region that can be plotted
#' by `x`. It merely indicates which part of the existing region should be
#' considered the focal region.
#'
#' @param x An `rdr_plotter` object.
#' @param left The genomic coordinate of the left edge of the focal region.
#' @param right The genomic coordinate of the right egde of the focal region.
#' @returns An `rdr_plotter` object.
#' @export
rdr_plotter_set_pad <- function(x, left = NULL, right = NULL) {
    x$left_pad_border <- left
    x$right_pad_border <- right

    x
}

rdr_plotter_plot_base <- function(x, main = NULL) {
    plot(
        NULL,
        main = main,
        xlim = range(x$bin_mids),
        ylim = c(0, 3),
        ylab = "Normalized Read Depth Ratio",
        xlab = "",
        xaxs = "i",
        xaxt = "n"
    )
}

rdr_plotter_plot_bg <- function(x) {
    for (bg in x$bg_samples) {
        graphics::lines(
            x$bin_mids,
            x$bincov_mat$rd[x$bins_to_plot, bg],
            col = "grey",
            lwd = 0.5
        )
    }
}

rdr_plotter_plot_segdups <- function(x) {
    if (length(x$segdups) > 0) {
        graphics::rect(
            GenomicRanges::start(x$segdups),
            0.1,
            GenomicRanges::end(x$segdups),
            0.2,
            col = "brown4",
            border = NA
        )
    }
}

rdr_plotter_plot_carriers <- function(x, carrier) {
    if (
        is.null(x$carriers) ||
            !all(names(x$carriers) %in% colnames(x$bincov_mat$rd))
    ) {
        return()
    }

    if (is.numeric(carrier) && length(carrier) == 1 && carrier == 0) {
        return()
    }

    if (is.null(carrier)) {
        carriers_to_plot <- x$carriers
    } else {
        carriers_to_plot <- x$carriers[carrier]
    }

    for (id in names(carriers_to_plot)) {
        graphics::lines(
            x$bin_mids,
            x$bincov_mat$rd[x$bins_to_plot, id],
            col = carriers_to_plot[[id]],
            lwd = 0.5
        )
    }
}

rdr_plotter_plot_pad <- function(x) {
    axlims <- graphics::par("usr")
    if ("left_pad_border" %in% names(x)) {
        graphics::rect(
            axlims[[1]],
            axlims[[3]],
            x$left_pad_border,
            axlims[[4]],
            col = "#FFAF0044",
            border = NA
        )
    }

    if ("right_pad_border" %in% names(x)) {
        graphics::rect(
            x$right_pad_border,
            axlims[[3]],
            axlims[[2]],
            axlims[[4]],
            col = "#FFAF0044",
            border = NA
        )
    }
}

rdr_plotter_plot_axis <- function(x) {
    axlims <- graphics::par("usr")
    xticks <- grDevices::axisTicks(axlims[1:2], log = FALSE, nint = 10)
    xlabs <- formatC(xticks, big.mark = ",", format = "d")
    graphics::axis(
        1,
        at = xticks,
        labels = FALSE,
        las = 2,
        cex.axis = 0.8,
        tcl = -1
    )
    graphics::mtext(xlabs, 1, at = xticks, adj = 1.05, cex = 0.6)
    graphics::title(
        xlab = sprintf("%s Position (bp)", x$gdcaller$gd$chr),
        line = 2
    )
}

rdr_plotter_plot_border <- function(x) {
    graphics::box(lwd = 2)
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
