# Functions for plotting binned coverage matrices.

# number of bins to use for smoothing the binned coverage values before
# plotting
SMOOTHING_WINDOW <- 31
# fraction of bins to plot
# lower values speed up plotting, but limit the visible CNV resolution
BIN_PLOT_FRACTION <- 0.05
# maximum number of background samples to plot
MAX_BACKGROUND <- 200L

GENE_LANES <- 16

#' Create a new `rdr_plotter` object.
#'
#' An `rdr_plotter` object represents a region of a binned coverage matrix to
#' visualize.
#'
#' @param x A `bincov_mat` object.
#' @param segdups A [`GRanges`][GenomicRanges::GRanges] object with ranges for
#'   segmental duplications.
#' @returns A `rdr_plotter` object.
#' @export
rdr_plotter <- function(x, segdups) {
    if (!inherits(x, "bincov_mat")) {
        stop("`x` must be a `bincov_mat` object")
    }

    if (!methods::is(segdups, "GRanges")) {
        stop("`segdups` must be a `GRanges` object")
    }

    validate_rdr_plotter(new_rdr_plotter(x, segdups))
}

new_rdr_plotter <- function(x, segdups) {
    mat <- smooth(x, SMOOTHING_WINDOW)
    bins_to_plot <- spaced_intervals(length(mat$ranges))
    bin_mids <- GenomicRanges::start(mat$ranges[bins_to_plot]) +
        (GenomicRanges::width(mat$ranges[bins_to_plot]) %/% 2)
    bg_samples <- sample(colnames(mat$rd), min(MAX_BACKGROUND, ncol(mat$rd)))
    span <- utils::head(mat$ranges, 1)
    GenomicRanges::end(span) <- utils::tail(GenomicRanges::end(mat$ranges), 1)
    ovp_sd <- segdups[S4Vectors::subjectHits(GenomicRanges::findOverlaps(
        span,
        segdups
    ))]
    ovp_genes <- genes_gr[S4Vectors::subjectHits(GenomicRanges::findOverlaps(
        span,
        genes_gr
    ))]
    ovp_genes <- split_gene_table(ovp_genes)

    structure(
        list(
            bincov_mat = mat,
            bins_to_plot = bins_to_plot,
            bin_mids = bin_mids,
            bg_samples = bg_samples,
            segdups = ovp_sd,
            genes = ovp_genes
        ),
        class = "rdr_plotter"
    )
}

validate_rdr_plotter <- function(x) {
    bincov_contigs <- unique(methods::as(
        GenomicRanges::seqnames(x$segdups),
        "character"
    ))
    segdup_contigs <- unique(methods::as(
        GenomicRanges::seqnames(x$bincov_mat$ranges),
        "character"
    ))

    if (!all(segdup_contigs %in% bincov_contigs)) {
        stop("all segmental duplication contigs must be in the bincov matrix")
    }

    x
}

#' Plot an `rdr_plotter`.
#'
#' @param x An `rdr_plotter` object.
#' @param y Kept for compatibility with the `plot` generic. Ignored here.
#' @param main An optional string to use as the title of the plot.
#' @param carrier Which carriers to plot, if `x` has carriers. If 0, no
#'   carriers are plotting. If `NULL`, all carriers are plotted. If an integer
#'   or character vector, it will be used as an index into the carriers vector.
#' @param ... Kept for compatibility with the `plot` generic. Ignored here.
#' @export
plot.rdr_plotter <- function(x, y = NULL, main = NULL, carrier = NULL, ...) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    rdr_plotter_plot_base(x, main)
    rdr_plotter_plot_bg(x)
    rdr_plotter_plot_carriers(x, carrier)
    rdr_plotter_plot_pad(x)
    rdr_plotter_plot_border(x)
    rdr_plotter_plot_segdups(x)
    rdr_plotter_plot_genes(x)
    rdr_plotter_plot_axis(x)
}

#' Create a new `gdplotter` object.
#'
#' A `gdplotter` is an object that can plot a genomic disorder (GD) region.
#'
#' @param gd A list-like object that represents a GD region.
#'   It must contain the following named components:
#'   * chr: Contig of the GD region
#'   * start_GRCh38: start of the GD region
#'   * end_GRCh38: end of the GD region
#'   * svtype: SV type of the GD (either 'DEL' or 'DUP')
#' @param bincov_mat A `bincov_mat` object over the GD region.
#' @param segdups A [`GRanges`][GenomicRanges::GRanges()] object of hg38
#'   segmental duplications.
#' @param carriers A character vector of sample IDs of the GD carriers.
#' @returns A `gdplotter` object.
#' @export
gdplotter <- function(gd, bincov_mat, segdups, carriers) {
    new_gdplotter(gd, bincov_mat, segdups, carriers)
}

new_gdplotter <- function(gd, bincov_mat, segdups, carriers) {
    carrier_col <- rep(
        if (gd$svtype == "DEL") "red" else "blue",
        length(carriers)
    )
    plotter <- new_rdr_plotter(bincov_mat, segdups) |>
        rdr_plotter_set_pad(gd$start_GRCh38, gd$end_GRCh38) |>
        rdr_plotter_set_carriers(stats::setNames(carrier_col, carriers))
    size_pretty <- format_size(gd$end_GRCh38 - gd$start_GRCh38 + 1)

    structure(
        list(
            gd = gd,
            carriers = carriers,
            plotter = plotter,
            size_pretty = size_pretty
        ),
        class = "gdplotter"
    )
}

#' @export
plot.gdplotter <- function(x, y = NULL, carrier = 1, ...) {
    carrier_pretty <- format_sample_id(x$carriers[[carrier]])
    main <- sprintf(
        "%s (hg38)\n%s (%s)",
        x$gd$GD_ID,
        carrier_pretty,
        x$size_pretty
    )

    plot(x$plotter, main = main, carrier = carrier)
}

#' Set the carriers for an `rdr_plotter`.
#'
#' @param x An `rdr_plotter` object.
#' @param carriers A named character vector of which the values are colors to
#'   use for plotting the carriers and the names are the IDs of the carriers.
#' @returns An `rdr_plotter` object.
#' @export
rdr_plotter_set_carriers <- function(x, carriers) {
    x$carriers <- carriers

    x
}

#' Set the padding for an `rdr_plotter`.
#'
#' The padding represents the genomic region outside the region of interest.
#' This function does not increase the size of the region that can be plotted
#' by `x`. It merely indicates which part of the existing region should be
#' considered the focal region.
#'
#' @param x An `rdr_plotter` object.
#' @param left The genomic coordinate of the left edge of the focal region.
#' @param right The genomic coordinate of the right egde of the focal region.
#' @returns An `rdr_plotter` object.
#' @export
rdr_plotter_set_pad <- function(x, left = NULL, right = NULL) {
    x$left_pad_border <- left
    x$right_pad_border <- right

    x
}

rdr_plotter_plot_base <- function(x, main = NULL) {
    graphics::layout(matrix(1:3, nrow = 3), width = 1, height = c(6, 1.5, 4))
    graphics::par(mar = c(0, 4.1, 4.1, 2.1))
    plot(
        NULL,
        main = main,
        xlim = range(x$bin_mids),
        ylim = c(0, 3),
        ylab = "Normalized Read Depth Ratio",
        xlab = "",
        xaxs = "i",
        xaxt = "n"
    )
}

rdr_plotter_plot_bg <- function(x) {
    for (bg in x$bg_samples) {
        graphics::lines(
            x$bin_mids,
            x$bincov_mat$rd[x$bins_to_plot, bg],
            col = "grey",
            lwd = 0.5
        )
    }
}

rdr_plotter_plot_segdups <- function(x) {
    graphics::par(mar = c(0, 4.1, 0, 2.1))
    plot(
        NULL,
        xlim = range(x$bin_mids),
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

rdr_plotter_plot_genes <- function(x) {
    graphics::par(mar = c(3.1, 4.1, 0, 2.1))
    plot(
        NULL,
        xlim = range(x$bin_mids),
        ylim = c(-1, GENE_LANES),
        xaxs = "i",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xlab = "",
        ylab = ""
    )
    if (length(x$genes) == 0) {
        return()
    }

    lane_ends <- integer(GENE_LANES)
    lane <- GENE_LANES
    lane_ends[[lane]] <- x$genes[[1]][1, ]$gene_end
    plot_gene(x$genes[[1]], lane)
    for (i in seq.int(2, length(x$genes))) {
        gene <- x$genes[[i]]
        gene_name_width <- graphics::strwidth(gene[1, ]$gene_name, cex = 0.8)
        display_start <- gene[1, ]$gene_start - gene_name_width * 1.1
        lane <- 0
        for (j in seq.int(GENE_LANES, 1)) {
            if (display_start > lane_ends[[j]]) {
                lane <- j
                break
            }
        }
        # maybe the gene should be plotted at bottom if it would overlap
        lane <- if (lane == 0) GENE_LANES else lane
        lane_ends[[lane]] <- max(lane_ends[[lane]], gene[1, ]$gene_end)
        plot_gene(gene, lane)
    }
}

rdr_plotter_plot_carriers <- function(x, carrier) {
    if (
        is.null(x$carriers) ||
            !all(names(x$carriers) %in% colnames(x$bincov_mat$rd))
    ) {
        return()
    }

    if (is.numeric(carrier) && length(carrier) == 1 && carrier == 0) {
        return()
    }

    if (is.null(carrier)) {
        carriers_to_plot <- x$carriers
    } else {
        carriers_to_plot <- x$carriers[carrier]
    }

    for (id in names(carriers_to_plot)) {
        graphics::lines(
            x$bin_mids,
            x$bincov_mat$rd[x$bins_to_plot, id],
            col = carriers_to_plot[[id]],
            lwd = 0.5
        )
    }
}

rdr_plotter_plot_pad <- function(x) {
    axlims <- graphics::par("usr")
    if ("left_pad_border" %in% names(x)) {
        graphics::rect(
            axlims[[1]],
            axlims[[3]],
            x$left_pad_border,
            axlims[[4]],
            col = "#FFAF0044",
            border = NA
        )
    }

    if ("right_pad_border" %in% names(x)) {
        graphics::rect(
            x$right_pad_border,
            axlims[[3]],
            axlims[[2]],
            axlims[[4]],
            col = "#FFAF0044",
            border = NA
        )
    }
}

rdr_plotter_plot_axis <- function(x) {
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
    graphics::title(
        xlab = sprintf(
            "%s Position (bp)",
            as.character(GenomicRanges::seqnames(x$bincov_mat$ranges[1]))
        ),
        line = 2
    )
}

rdr_plotter_plot_border <- function(x) {
    graphics::box(lwd = 2)
}

plot_gene <- function(x, lane) {
    lane <- lane - 1
    graphics::lines(
        c(x[1, ]$gene_start, x[1, ]$gene_end),
        c(lane, lane),
        col = x[1, ]$color
    )
    graphics::rect(
        x$block_starts,
        ifelse(x$thick, lane - 0.4, lane - 0.3),
        x$block_ends,
        ifelse(x$thick, lane + 0.4, lane + 0.3),
        col = x$color,
        border = NA
    )

    # don't print ENSG gene names to save plot space
    if (!grepl("^ENSG[0-9]+", x[1, ]$gene_name)) {
        graphics::text(
            x[1, ]$gene_start,
            lane,
            labels = x[1, ]$gene_name,
            pos = 2,
            offset = 0.1,
            cex = 0.8
        )
    }
}

# Generate a vector of best-effort equally spaced integers from 1 to n.
spaced_intervals <- function(n) {
    nintervals <- ceiling(n * BIN_PLOT_FRACTION)
    if (nintervals < 2) {
        return(seq_len(n))
    }
    gaps <- nintervals - 1L
    total_gap_size <- n - nintervals
    remaining_pad <- total_gap_size %% gaps
    pads <- rep(floor(total_gap_size / gaps), gaps)
    i <- seq_len(remaining_pad)
    pads[i] <- pads[i] + 1L
    steps <- c(1, pads + 1L)

    cumsum(steps)
}

# Format a length in bases to a pretty string.
format_size <- function(size) {
    if (size <= 1000) {
        size_pretty <- sprintf("%d b", size)
    } else if (size <= 100000) {
        size_pretty <- sprintf("%.2f Kb", size / 1000)
    } else {
        size_pretty <- sprintf("%.2f Mb", size / 1000000)
    }

    size_pretty
}

# Format a sample ID into a pretty string for plotting.
format_sample_id <- function(x) {
    if (nchar(x) > 30) {
        paste0(substr(x, 1, 30), "...")
    } else {
        x
    }
}

split_gene_table <- function(x) {
    if (length(x) == 0) {
        return(list())
    }

    mapply(
        split_gene_blocks,
        as.character(GenomicRanges::seqnames(x)),
        S4Vectors::start(x),
        S4Vectors::end(x),
        S4Vectors::mcols(x)$thickStart,
        S4Vectors::mcols(x)$thickEnd,
        S4Vectors::mcols(x)$chromStarts,
        S4Vectors::mcols(x)$blockSizes,
        S4Vectors::mcols(x)$geneName,
        S4Vectors::mcols(x)$color,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )
}

split_gene_blocks <- function(
    chrom,
    chrom_start,
    chrom_end,
    thick_start,
    thick_end,
    offsets,
    blocks,
    gene_name,
    color
) {
    offsets <- as.integer(strsplit(
        sub(",$", "", offsets),
        split = ",",
        fixed = TRUE
    )[[1]])
    blocks <- as.integer(strsplit(
        sub(",$", "", blocks),
        split = ",",
        fixed = TRUE
    )[[1]])
    starts <- chrom_start + offsets
    ends <- starts + blocks - 1
    if (thick_start <= thick_end) {
        i <- utils::tail(which(starts <= thick_start), 1)
        j <- utils::tail(which(starts <= thick_end), 1)
        starts <- c(
            if (i > 1) starts[1:(i - 1)] else integer(),
            if (starts[[i]] < thick_start) starts[[i]] else integer(),
            thick_start,
            if (j > i) starts[(i + 1):j] else integer(),
            if (thick_end < ends[[j]]) thick_end + 1 else integer(),
            if (j < length(starts)) {
                starts[(j + 1):length(starts)]
            } else {
                integer()
            }
        )
        ends <- c(
            if (i > 1) ends[1:(i - 1)] else integer(),
            if (starts[[i]] < thick_start) thick_start - 1 else integer(),
            if (i == j) thick_end else ends[[i]],
            if (j > i + 1) ends[(i + 1):(j - 1)] else integer(),
            if (j > i) thick_end else integer(),
            if (thick_end < ends[[j]]) ends[[j]] else integer(),
            if (j < length(ends)) ends[(j + 1):length(ends)] else integer()
        )
        thick <- rep(FALSE, length(starts))
        thick[starts >= thick_start & ends <= thick_end] <- TRUE
    } else {
        thick <- rep(FALSE, length(starts))
    }

    data.table(
        chr = chrom,
        gene_start = chrom_start,
        gene_end = chrom_end,
        block_starts = starts,
        block_ends = ends,
        thick = thick,
        gene_name = gene_name,
        color = color
    )
}
