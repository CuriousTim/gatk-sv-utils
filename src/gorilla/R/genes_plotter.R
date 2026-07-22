#' @export
plot.genes_plotter <- function(x, y, gene_lanes = 16, ...) {
    # maybe it would be better to set the y-axis limits some fraction larger
    # than the number of gene lanes
    plot(
        NULL,
        xlim = c(x$region$start, x$region$end),
        ylim = c(-1, gene_lanes),
        xaxs = "i",
        xaxt = "n",
        yaxs = "i",
        yaxt = "n",
        xlab = "",
        ylab = ""
    )

    if (length(x$genes) > 0) {
        invisible(plot_genes(x$genes, gene_lanes))
    }
}

# construct an object that can plot genes
new_genes_plotter <- function(contig, start, end) {
    query <- GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
    ovps <- GenomicRanges::findOverlaps(query, genes_gr)
    genes <- genes_gr[S4Vectors::subjectHits(ovps)] |> split_gene_table()

    structure(
        list(
            genes = genes,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "genes_plotter"
    )
}

plot_genes <- function(x, nlanes) {
    lane_ends <- integer(nlanes)
    # y-axis of the plotting region ranges from -1 to nlanes, but the genes are
    # only drawn on lanes 0 to nlanes - 1, so that genes at the top and bottom
    # lanes don't get clipped (they might clipped at unusual plot heights and
    # number of gene lanes...)
    # offset is the number of lanes from the top lane + 1
    offset <- 1
    gene_name_width <- graphics::strwidth(x[[1]][1, ]$gene_name, cex = 0.8)
    display_start <- x[[1]][1, ]$gene_start - gene_name_width * 1.1
    lane_ends[[offset]] <- x[[1]][1, ]$gene_end
    draw_gene(x[[1]], nlanes - offset, display_start)
    if (length(x) == 1) {
        return()
    }
    for (i in seq.int(from = 2, to = length(x), by = 1)) {
        gene <- x[[i]]
        gene_name_width <- graphics::strwidth(gene[1, ]$gene_name, cex = 0.8)
        display_start <- gene[1, ]$gene_start - gene_name_width * 1.1
        offset <- 0
        for (j in seq_along(lane_ends)) {
            if (display_start > lane_ends[[j]]) {
                offset <- j
                break
            }
        }

        offset <- if (offset == 0) 1 else offset
        lane_ends[[offset]] <- max(lane_ends[[offset]], gene[1, ]$gene_end)
        draw_gene(gene, nlanes - offset, display_start)
    }
}

draw_gene <- function(x, lane, display_start) {
    graphics::lines(
        c(x[1, ]$gene_start, x[1, ]$gene_end),
        c(lane, lane),
        col = x[1, ]$color
    )
    graphics::rect(
        x$block_starts,
        ifelse(x$thick, lane - 0.4, lane - 0.25),
        x$block_ends,
        ifelse(x$thick, lane + 0.4, lane + 0.25),
        col = x$color,
        border = NA
    )

    xmin <- graphics::par("usr")[[1]]
    # don't print ENSG gene names to save space
    if (!grepl("^ENSG[0-9]+", x[1, ]$gene_name)) {
        # plotting gene name to left of gene can cause them to be truncated
        graphics::text(
            if (display_start <= xmin) xmin else x[1, ]$gene_start,
            lane,
            labels = x[1, ]$gene_name,
            pos = 2,
            offset = 0.1,
            cex = 0.8,
            xpd = NA
        )
    }
}

# Split a GRanges object of genes into one data.table per gene with segmented
# blocks.
split_gene_table <- function(x) {
    if (length(x) == 0) {
        return(list())
    }

    mapply(
        segment_gene_blocks,
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

# The gene annotation data retrieved from the UCSC Genome Browser is organized
# as start and end coordinates of the gene and blocks within that region
# corresponding to the exons of the gene. For protein-coding genes, there are
# additionally coordinates to mark the start and end of the protein-coding
# region of the gene (including intervening introns). The key point here is
# that for protein-coding genes, the start and end of the protein-coding
# region, referred to as 'thickStart' and 'thickEnd', do not neccessarily align
# to the block boundaries. Non-protein-coding regions of the gene are drawn
# thinner than protein-coding regions of the gene meaning that two sections of
# a block might need to be drawn with differing thicknesses. To facilitate
# drawing the blocks with the appropriate thickness, it is useful to segment
# the blocks so that each resulting segment is a single thickness that is what
# this function does.
segment_gene_blocks <- function(
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
