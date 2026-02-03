PE_LANES_PER_SAMPLE <- 16

#' Create a new `pe_plotter` object
#'
#' A `pe_plotter` can create a visualization of discordant paired-end reads
#' alignments.
#'
#' @param x A [`pe_mat`][query()] object.
#' @returns A `pe_plotter` object.
#' @export
#'
#' @examples
#' pe_path <- system.file("extdata", "example.PE.txt.gz", package = "gorilla", mustWork = TRUE)
#' pe <- pe_file(pe_path)
#' mat <- query(pe, "chr16", 28743149, 28745149)
#' plotter <- pe_plotter(mat)
pe_plotter <- function(x) {
    new_pe_plotter(x)
}

#' Plot a `pe_plotter` object
#'
#' The `plot()` method for `pe_plotter`.
#'
#' The only additional parameter that is accepted is `col`, which should be a
#' string or a character vector equal to the length of `samples` giving the
#' colors to use for each sample.
#'
#' @param x A [`pe_plotter`][pe_plotter()] object.
#' @param y Ignored.
#' @param samples Samples to plot. The samples will be plotted from top to
#'   bottom in the order they are given.
#' @param main Title for the plot.
#' @param col Color(s) to use for `samples`. Either a string or a character
#'   vector of hexadecimal colors equal to the length of `samples`
#' @param ... Ignored.
#' @name plot.pe_plotter
#' @export
#'
#' @examples
#' pe_path <- system.file("extdata", "example.PE.txt.gz", package = "gorilla", mustWork = TRUE)
#' pe <- pe_file(pe_path)
#' mat <- query(pe, "chr16", 28743149, 28745149)
#' plotter <- pe_plotter(mat)
#' plot(plotter, samples = "gorilla0003")
plot.pe_plotter <- function(x, y, samples, main = "", col = "#000000", ...) {
    stopifnot(is.character(samples) && length(samples) > 0)
    stopifnot(is_string(main))
    stopifnot(
        is.character(col) &&
            (length(col) == 1 || length(col) == length(samples))
    )

    lanes <- length(samples) * PE_LANES_PER_SAMPLE
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(1, lanes),
        ylab = "PE Evidence",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        yaxt = "n",
        main = main
    )

    col <- if (length(col) == 1) rep(col, length(sample)) else col
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
    rstart <- NULL
    mstart <- NULL
    mat <- x$mat[
        rcontig == mcontig & rstart >= x$region$start & mstart <= x$region$end,
    ]
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
    rcontig <- NULL
    mcontig <- NULL
    rstart <- NULL
    mstart <- NULL
    rstrand <- NULL
    mstrand <- NULL
    rarrstart <- NULL
    marrstart <- NULL
    ..arrow_len <- NULL

    x[,
        c("rarrstart", "marrstart") := list(
            rstart + ifelse(rstrand == "+", -..arrow_len, ..arrow_len),
            mstart + ifelse(mstrand == "+", -..arrow_len, ..arrow_len)
        )
    ]
    x[
        rcontig == mcontig,
        c("pe_start", "pe_end") := list(
            pmin(rstart, rarrstart),
            pmax(mstart, marrstart)
        )
    ]

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
