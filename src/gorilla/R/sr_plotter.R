#' Create a new `sr_plotter` object
#'
#' A `sr_plotter` can create a visualization of split read alignments.
#'
#' @param x A [`sr_mat`][query()] object.
#' @returns A `sr_plotter` object.
#' @export
#'
#' @examples
#' sr_path <- system.file("extdata", "example.SR.txt.gz", package = "gorilla", mustWork = TRUE)
#' sr <- sr_file(sr_path)
#' mat <- query(sr, "chr16", 28743149, 28745149)
#' plotter <- sr_plotter(mat)
sr_plotter <- function(x) {
    new_sr_plotter(x)
}

#' Plot a `sr_plotter` object
#'
#' The `plot()` method for `sr_plotter`.
#'
#' The only additional parameter that is accepted is `col`, which should be a
#' string or a character vector equal to the length of `samples` giving the
#' colors to use for each sample.
#'
#' @param x A [`sr_plotter`] object.
#' @param y Ignored.
#' @param samples Samples to plot. The samples will be plotted from top to
#'   bottom in the order they are given.
#' @param main Title for the plot.
#' @param col Color(s) to use for `samples`. Either a string or a character
#'   vector of hexadecimal colors equal to the length of `samples`
#' @param ... Ignored.
#' @name plot.sr_plotter
#' @export
#'
#' @examples
#' sr_path <- system.file("extdata", "example.SR.txt.gz", package = "gorilla", mustWork = TRUE)
#' sr <- sr_file(sr_path)
#' mat <- query(sr, "chr16", 28743149, 28745149)
#' plotter <- sr_plotter(mat)
#' plot(plotter, samples = "gorilla0003")
plot.sr_plotter <- function(x, y, samples, main = "", col = "#000000", ...) {
    stopifnot(is.character(samples) && length(samples) > 0)
    stopifnot(is_string(main))
    stopifnot(
        is.character(col) &&
            (length(col) == 1 || length(col) == length(samples))
    )

    lanes <- length(samples) * 2
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(0.9, lanes + 0.1),
        ylab = "SR Evidence",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        yaxt = "n",
        main = main
    )

    col <- if (length(col) == 1) rep(col, length(samples)) else col
    col <- validate_sr_colors(col)

    for (i in seq_along(samples)) {
        target <- x$mat[samples[[i]], nomatch = NULL]
        contig <- NULL
        pos <- NULL
        data.table::setkey(target, contig, pos)
        if (nrow(target) > 0) {
            draw_sample_sr(
                target,
                lanes - (2 * (i - 1)),
                col[[i]]
            )
        }
    }
}

new_sr_plotter <- function(x) {
    stopifnot(inherits(x, "sr_mat"))

    count <- NULL
    mat <- x$mat[count > 1, ]
    mat <- set_sr_alpha(mat)

    structure(list(mat = mat, region = x$region), class = "sr_plotter")
}

draw_sample_sr <- function(x, y, color) {
    side <- NULL
    left_clip <- x[side == "left", ]
    if (nrow(left_clip) > 0) {
        graphics::points(
            x = left_clip$pos,
            y = rep(y - 1, nrow(left_clip)),
            col = paste0(color, left_clip$alpha),
            pch = "\\",
        )
    }
    right_clip <- x[side == "right", ]
    if (nrow(right_clip) > 0) {
        graphics::points(
            x = right_clip$pos,
            y = rep(y, nrow(right_clip)),
            col = paste0(color, right_clip$alpha),
            pch = "/"
        )
    }
}

set_sr_alpha <- function(x) {
    alpha <- NULL
    count <- NULL

    # alpha decimal values converted to hex
    # 0.4 = 66
    # 0.6 = 99
    # 0.8 = CC
    x[,
        alpha := data.table::fcase(
            count == 2 , "66" , count == 3 , "99" , count >= 4 , "cc"
        )
    ]

    x
}

validate_sr_colors <- function(x) {
    if (!all(grepl("#[[:xdigit:]]{6}", x))) {
        stop("colors must be of the form #rrggbb")
    }

    x
}
