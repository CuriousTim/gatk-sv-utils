PE_LANES_PER_SAMPLE <- 15

#' Create a new `pe_plotter` object.
#'
#' A `pe_plotter` can create a visualization of discordant paired-end reads
#' alignments.
#'
#' @param x A [`pe_mat`][query()] object.
#' @returns A `pe_plotter` object.
#' @export
pe_plotter <- function(x) {
    new_pe_plotter(x)
}

#' Plot a `pe_plotter` object.
#'
#' @param x A [`pe_plotter`] object.
#' @param y Kept for compatability with generic, but ignored here.
#' @param samples Samples to plot. The samples will be plotted from top to
#'   bottom in the order they are given.
#' @param ... Other plotting parameters. Currently the only accepted parameter
#'   is `col` which should either be a character vector of length 1 or length
#'   equal to `samples` giving the colors to use for each sample.
#' @name plot.pe_plotter
#' @export
plot.pe_plotter <- function(x, y, samples, ...) {
    stopifnot(is.character(samples) && length(samples) > 0)

    lanes <- length(samples) * PE_LANES_PER_SAMPLE
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(1, lanes),
        ylab = "PE Evidence",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        yaxt = "n"
    )

    col <- tryCatch(
        color_from_dots(length(samples), ...),
        value_error = function(e) stop("length of `col` must be 1 or equal to length of `samples`")
    )

    for (i in seq_along(samples)) {
        target <- x$mat[samples[[i]], nomatch = NULL]
        rcontig <- NULL
        rstart <- NULL
        data.table::setkey(target, rcontig, rstart)
        if (nrow(target) > 0) {
            draw_sample_pe(
                target,
                lanes - (PE_LANES_PER_SAMPLE * (i - 1)),
                col[[i]]
            )
        }
    }
}

new_pe_plotter <- function(x) {
    stopifnot(inherits(x, "pe_mat"))

    rcontig <- NULL
    mcontig <- NULL
    mat <- x$mat[rcontig == mcontig & rstart >= x$region$start & mstart <= x$region$end, ]
    set_pe_arrow_widths(mat, x$region)

    structure(list(mat = mat, region = x$region), class = "pe_plotter")
}


draw_sample_pe <- function(x, y, color) {
    rcontig <- NULL
    mcontig <- NULL
    same_contig <- x[rcontig == mcontig, ]
    ends <- integer(PE_LANES_PER_SAMPLE)
    # assumes pairs is sorted by rstart
    for (i in seq_len(nrow(same_contig))) {
        pair <- same_contig[i, ]
        j <- 0
        for (k in seq_along(ends)) {
            if (pair$pe_start >= ends[[k]]) {
                j <- k
                break
            }
        }

        # skip drawing read pairs that would overlap an existing one
        # could be terrible
        if (j == 0) {
            next
        }

        ends[[j]] <- pair$pe_end
        draw_pe_pair(pair, y - (j - 1), color)
    }
}

set_pe_arrow_widths <- function(x, region) {
    arrow_len <- (region$end - region$start + 1) * 0.01
    rstart <- NULL
    mstart <- NULL
    rstrand <- NULL
    mstrand <- NULL
    x[, c("rarrstart", "marrstart") := list(
        rstart + ifelse(rstrand == "+", -..arrow_len, ..arrow_len),
        mstart + ifelse(mstrand == "+", -..arrow_len, ..arrow_len)
    )]
    x[rcontig == mcontig, c("pe_start", "pe_end") := list(
        pmin(rstart, rarrstart),
        pmax(mstart, marrstart)
    )]

    x
}

draw_pe_pair <- function(pair, y, col) {
    graphics::arrows(
        x0 = pair$rarrstart,
        y0 = y,
        x1 = pair$rstart,
        length = 0.05,
        col = col
    )
    graphics::arrows(
        x0 = pair$marrstart,
        y0 = y,
        x1 = pair$mstart,
        length = 0.05,
        col = col
    )
    graphics::segments(pair$pe_start, y, pair$pe_end, y, lty = 3, col = col)
}
