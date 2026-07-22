new_segdups_plotter <- function(contig, start, end) {
    query <- GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
    ovps <- GenomicRanges::findOverlaps(query, segdups_gr)
    segdups <- segdups_gr[S4Vectors::subjectHits(ovps)] |>
        GenomicRanges::reduce()

    structure(
        list(
            segdups = segdups,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "segdups_plotter"
    )
}

#' @export
plot.segdups_plotter <- function(x, y, ...) {
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(0, 1),
        xaxs = "i",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = ""
    )

    if (length(x$segdups) > 0) {
        graphics::rect(
            GenomicRanges::start(x$segdups),
            0.1,
            GenomicRanges::end(x$segdups),
            0.9,
            col = "brown4",
            border = NA
        )
    }
}
