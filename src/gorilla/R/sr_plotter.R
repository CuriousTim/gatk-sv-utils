#' Create a new `sr_plotter` object.
#'
#' A `sr_plotter` can create a visualization of split read alignments.
#'
#' @param x A [`sr_mat`][query()] object.
#' @returns A `sr_plotter` object.
#' @export
sr_plotter <- function(x) {
    new_sr_plotter(x)
}

#' Plot a `sr_plotter` object.
#'
#' @param x A [`sr_plotter`] object.
#' @param y Kept for compatability with generic, but ignored here.
#' @param samples Samples to plot. The samples will be plotted from top to
#'   bottom in the order they are given.
#' @param ... Other plotting parameters. Currently the only accepted parameter
#'   is `col` which should either be a character vector of length 1 or length
#'   equal to `samples` giving the colors to use for each sample.
#' @name plot.sr_plotter
#' @export
plot.sr_plotter <- function(x, y, samples, ...) {
    stopifnot(is.character(samples) && length(samples) > 0)

    lanes <- length(samples) * 2
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(1, lanes),
        ylab = "SR Evidence",
        xlab = "",
        xaxs = "i",
        xaxt = "n",
        yaxt = "n"
    )

    col <- tryCatch(
        color_from_dots(length(samples), ...),
        value_error = function(e) stop("length of `col` must be 1 or equal to length of `samples`")
    )
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
            y = rep(y, nrow(left_clip)),
            col = paste0(color, left_clip$alpha),
            pch = "\\",
        )
    }
    right_clip <- x[side == "right", ]
    if (nrow(right_clip) > 0) {
        graphics::points(
            x = right_clip$pos,
            y = rep(y - 1, nrow(right_clip)),
            col = paste0(color, right_clip$alpha),
            pch = "/"
        )
    }
}

set_sr_alpha <- function(x) {
    # alpha decimal values converted to hex
    # 0.4 = 66
    # 0.6 = 99
    # 0.8 = CC
    x[, alpha := data.table::fcase(count == 2, "66", count == 3, "99", count >= 4, "cc")]

    x
}

validate_sr_colors <- function(x) {
    if (!all(grepl("#[[:xdigit:]]{6}", x))) {
        stop("colors must be of the form #rrggbb")
    }

    x
}
