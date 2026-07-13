new_repeats_plotter <- function(contig, start, end) {
    query <- GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
    # supress warnings about not having sequence levels in common
    sd_ovps <- suppressWarnings(GenomicRanges::findOverlaps(query, segdups_gr))
    sr_ovps <- suppressWarnings(GenomicRanges::findOverlaps(
        query,
        simple_repeats_gr
    ))
    rm_ovps <- suppressWarnings(GenomicRanges::findOverlaps(
        query,
        repeat_mask_gr
    ))

    segdups <- segdups_gr[S4Vectors::subjectHits(sd_ovps)] |>
        GenomicRanges::reduce()
    simple_repeats <- simple_repeats_gr[S4Vectors::subjectHits(sr_ovps)] |>
        GenomicRanges::reduce()
    repeat_mask <- repeat_mask_gr[S4Vectors::subjectHits(rm_ovps)] |>
        GenomicRanges::reduce()

    structure(
        list(
            segdups = segdups,
            simple_repeats = simple_repeats,
            repeat_mask = repeat_mask,
            region = list(conitg = contig, start = start, end = end)
        ),
        class = "repeats_plotter"
    )
}

#' @export
plot.repeats_plotter <- function(x, y, ...) {
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(0, 11),
        xaxs = "i",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = ""
    )

    if (length(x$segdups) > 0) {
        graphics::rect(
            GenomicRanges::start(x$segdups),
            1,
            GenomicRanges::end(x$segdups),
            4,
            col = "#8B2323",
            border = NA
        )
    }

    if (length(x$repeat_mask) > 0) {
        graphics::rect(
            GenomicRanges::start(x$repeat_mask),
            4,
            GenomicRanges::end(x$repeat_mask),
            7,
            col = "#238b23",
            border = NA
        )
    }

    if (length(x$simple_repeats) > 0) {
        graphics::rect(
            GenomicRanges::start(x$simple_repeats),
            7,
            GenomicRanges::end(x$simple_repeats),
            10,
            col = "#23238b",
            border = NA
        )
    }

    graphics::text(
        rep(x$region$end - (x$region$end - x$region$start + 1) * 0.001, 3),
        c(2.5, 5.5, 8.5),
        labels = c("SD", "RM", "SR"),
        adj = c(1, 0.5),
        cex = 0.7,
    )
}
