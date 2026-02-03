# colors for child, father, mother
TRIO_COLORS <- c("#bc5d41", "#84a955", "#965da7")

#' Create a new `denovo_plotter` object
#'
#' A `denovo_plotter` can create a visualization of de novo SV evidence for a
#' child, father, mother trio.
#'
#' The plot will display evidence from discordant read pairs, split reads and
#' read depth, along with some annotations.
#'
#' @param x A [`svtrio`] object.
#' @returns A `denovo_plotter` object.
#' @export
#' @examples
#' pe_path <- system.file(
#'     "extdata",
#'     "example.PE.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' sr_path <- system.file(
#'     "extdata",
#'     "example.SR.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' rd_path <- system.file(
#'     "extdata",
#'     "example.RD.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' medians_path <- system.file(
#'     "extdata",
#'     "example_medianCov.txt",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#'
#' pe <- pe_file(pe_path)
#' sr <- sr_file(sr_path)
#' rd_medians <- read_median_coverages(medians_path)
#' rd <- rd_file(rd_path, rd_medians)
#'
#' sv <- svevidence("chr16", 28743149, 28745149, pe, sr, rd, "DUP")
#' trio <- svtrio(sv, "gorilla0000", "gorilla0001", "gorilla0002")
#' plotter <- denovo_plotter(trio)
denovo_plotter <- function(x) {
    new_denovo_plotter(x)
}

#' @export
plot.denovo_plotter <- function(x, y, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(1:5, nrow = 5), heights = c(1.2, 0.3, 0.4, 0.1, 1))
    # TODO: add evidence types of the call
    main <- sprintf(
        "de novo SV (%s)\n%s (%s)",
        x$svtype,
        pretty_sample_id(x$trio$child),
        pretty_size(x$region$end - x$region$start + 1)
    )

    trio <- unlist(x$trio)
    graphics::par(mar = c(0, 4.1, 4.1, 2.1), mgp = c(1, 1, 0))
    plot(x$pe, samples = trio, main = main, col = TRIO_COLORS)
    draw_padding(x$region$start, x$region$end)

    graphics::par(mar = c(0, 4.1, 0, 2.1), mgp = c(1, 1, 0))
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
