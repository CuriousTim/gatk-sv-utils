# maximum number of background samples to plot
MAX_BACKGROUND <- 200L

new_gd_plotter <- function(x) {
    stopifnot(inherits(x, "gd_caller"))
    if (is.null(x[["carriers"]])) {
        stop("at least one GD carrier must be present")
    }

    rd <- new_rd_plotter(x$rd)
    genes <- new_genes_plotter(
        x$rd$region$contig,
        x$rd$region$start,
        x$rd$region$end
    )
    segdups <- new_segdups_plotter(
        x$rd$region$contig,
        x$rd$region$start,
        x$rd$region$end
    )
    bg_samples <- sample(
        colnames(x$rd$mat),
        min(MAX_BACKGROUND, ncol(x$rd$mat))
    )

    structure(
        list(gd = x$gd, rd = rd, genes = genes, segdups = segdups, bg_samples),
        class = "gd_plotter"
    )
}

#' @export
plot.gd_plotter <- function(x, y, carrier = 1, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(1:3, nrow = 3), width = 1, heights = c(6, 0.5, 4))

    carrier_pretty <- pretty_sample_id(x$gd$carriers[[carrier]])
    size_pretty <- pretty_size(x$gd$end_GRCh38 - x$gd$start_GRCh38 + 1)

    main <- sprintf(
        "%s (hg38)\n%s (%s)",
        x$gd$GD_ID,
        carrier_pretty,
        size_pretty
    )

    graphics::par(mar = c(0, 4.1, 4.1, 2.1))
    plot(
        x$rd,
        samples = x$gd$carriers[[1]],
        bg_samples = x$bg_samples,
        main = main,
        col = if (x$gd$svtype == "DEL") "red" else "blue"
    )
    draw_padding(x$gd$rd$region$start, x$gd$rd$region$end)
    graphics::box(lwd = 2)

    graphics::par(mar = c(0, 4.1, 0, 2.1))
    plot(x$segdups)

    graphics::par(mar = c(3.1, 4.1, 0, 2.1))
    plot(x$genes)
    draw_genomic_coords_axis(x$gd$rd$region$contig)
}
