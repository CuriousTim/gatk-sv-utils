#' Read a BEDX file with segmental duplication regions.
#'
#' The first three columns must be contig, start, and end in 0-start, open-end
#' coordinates.
#'
#' @param path Path to the file.
#' @returns A `GRanges` object.
read_segdups <- function(path) {
    V1 <- NULL
    V2 <- NULL
    V3 <- NULL

    tmp <- fread(path, header = FALSE, select = 1:3, sep = "\t")
    tmp <- tmp[V1 %in% HG38_PRIMARY_CONTIGS, ]
    tmp[, c("V2", "V3") := list(as.integer(V2), as.integer(V3))]
    tmp[, V2 := V2 + 1L]

    tmp <- unique(tmp)

    GenomicRanges::reduce(GenomicRanges::GRanges(
        tmp$V1,
        IRanges::IRanges(tmp$V2, tmp$V3)
    ))
}
