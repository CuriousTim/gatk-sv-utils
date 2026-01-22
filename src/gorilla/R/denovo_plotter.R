# colors for child, father, mother
TRIO_COLORS <- c("#bc5d41", "#84a955", "#965da7")

#' Create a new `denovo_plotter` object.
#'
#' A `denovo_plotter` can create a visualization of de novo SV evidence for a
#' child, father, mother trio.
#'
#' @param x A [`svtrio`] object.
#' @returns A `denovo_plotter` object.
#' @export
denovo_plotter <- function(x) {
    new_denovo_plotter(x)
}

#' Plot a `denovo_plotter` object.
#'
#' @param x A [`denovo_plotter`] object.
#' @param y Kept for compatibility with generic, but ignored here.
#' @param ... Kept for compatibility with generic, but ignored here.
#' @name plot.denovo_plotter
#' @export
plot.denovo_plotter <- function(x, y, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(1:6, nrow = 6), heights = c(0.1, 0.75, 0.25, 0.25, 0.1, 1))
    # draw the plot title inside of a plot
    graphics::par(mar = c(0, 0, 0, 0))
    plot(
        NULL,
        xlim = c(0, 1),
        ylim = c(0, 1),
        ylab = "",
        xlab = "",
        xaxs = "i",
        yaxs = "i",
        xaxt = "n",
        yaxt = "n",
        frame.plot = FALSE
    )
    # TODO: add SV type and evidence types of the call
    main <- sprintf(
        "de novo SV (%s)\n%s (%s)",
        x$svtype,
        pretty_sample_id(x$trio$child),
        pretty_size(x$region$end - x$region$start + 1)
    )
    text(x = 0.5, y = 0.5, labels = main, adj = c(0.5, 0.5))

    trio <- unlist(x$trio)
    graphics::par(mar = c(0, 4.1, 0, 2.1), mgp = c(1, 1, 0))
    plot(x$pe, samples = trio, col = TRIO_COLORS)
    draw_padding(x$region$start, x$region$end)

    plot(x$sr, samples = trio, col = TRIO_COLORS)
    draw_padding(x$region$start, x$region$end)

    graphics::par(mgp = c(3, 1, 0))
    plot(x$rd, samples = trio, col = TRIO_COLORS)
    draw_padding(x$region$start, x$region$end)

    plot(x$segdups)

    graphics::par(mar = c(3.1, 4.1, 0, 2.1))
    plot(x$genes)
    draw_genomic_coords_axis(x$region$contig)
}

new_denovo_plotter <- function(x) {
    stopifnot(inherits(x, "svtrio"))

    pe <- new_pe_plotter(x$evidence$pe)
    sr <- new_sr_plotter(x$evidence$sr)
    rd <- new_rd_plotter(x$evidence$rd)

    genes <- new_genes_plotter(
        x$evidence$region$contig,
        x$evidence$region$qstart,
        x$evidence$region$qend
    )

    segdups <- new_segdups_plotter(
        x$evidence$region$contig,
        x$evidence$region$qstart,
        x$evidence$region$qend
    )

    structure(
        list(
            pe = pe,
            sr = sr,
            rd = rd,
            genes = genes,
            segdups = segdups,
            svtype = x$evidence$svtype,
            region = x$evidence$region,
            trio = x$trio
        ),
        class = "denovo_plotter"
    )
}
