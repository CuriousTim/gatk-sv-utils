is_string <- function(x) {
    is.character(x) && length(x) == 1
}

draw_padding <- function(inner_start, inner_end) {
    axlims <- graphics::par("usr")
    graphics::rect(
        axlims[[1]],
        axlims[[3]],
        inner_start,
        axlims[[4]],
        col = "#FFAF0044",
        border = NA
    )

    graphics::rect(
        inner_end,
        axlims[[3]],
        axlims[[2]],
        axlims[[4]],
        col = "#FFAF0044",
        border = NA
    )
}

draw_genomic_coords_axis <- function(contig) {
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
    graphics::title(xlab = sprintf("%s Position (bp)", contig), line = 2)
}

# Format a sample ID into a pretty string for plotting.
pretty_sample_id <- function(x) {
    if (nchar(x) > 30) {
        paste0(substr(x, 1, 30), "...")
    } else {
        x
    }
}

# Format a length in bases to a pretty string.
pretty_size <- function(size) {
    if (size <= 1000) {
        sprintf("%d b", size)
    } else if (size <= 100000) {
        sprintf("%.2f Kb", size / 1000)
    } else {
        sprintf("%.2f Mb", size / 1000000)
    }
}

color_from_dots <- function(n, ...) {
    dots <- list(...)
    col <- dots[["col"]]
    if (!is.null(col) && !is.character(col)) {
        stop(errorCondition("`col` must be `NULL` or a character vector", class = "type_error"))
    }

    if (is.null(col)) {
        rep("#000000", n)
    } else if (length(col) == 1) {
        rep(col, n)
    } else if (length(col) == n) {
        col
    } else {
        stop(errorCondition("", class = "value_error"))
    }
}
