is_string <- function(x) {
    is.character(x) && length(x) == 1
}

is_number <- function(x) {
    is.numeric(x) && length(x) == 1
}

# Add colored rectangles over the left and right sides of a plot where the
# x-axis is genomic coordinates. The rectangles represent regions of the genome
# flanking the region of interest.
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

# Draw an axis on the bottom of a plot showing genomic coordinates, including
# an axis label of the form "<contig> Position (bp)".
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
    graphics::title(xlab = sprintf("%s Position (bp)", contig), line = 1.5)
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
