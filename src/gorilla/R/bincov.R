#' Retrieve the header of a bincov matrix.
#'
#' @param con `TabixFile` Connection object to bincov matrix.
#' @returns `character` Header of the matrix.
#' @export
read_bincov_header <- function(con) {
    header <- Rsamtools::headerTabix(con)[["header"]]
    if (is.null(header)) {
        stop("bincov matrix is missing a header")
    }

    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) <= 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop("bincov matrix has an invalid header")
    }
    header[1:3] <- c("chr", "start", "end")

    header
}

#' Query ranges from a bincov matrix.
#'
#' @param con `TabixFile` Connection object to bincov matrix.
#' @param ranges `GRanges` Genomic ranges to use as the query.
#' @param header `character` The parsed header line from the bincov matrix.
#' @param select `integer` Which columns to keep.
#' @returns `list` A list of `data.table`s, one for each range in `ranges`. If
#'   a range does not overlap any intervals in the bincov matrix, the element
#'   corresponding to that range will be `NULL`.
#' @export
query_bincov <- function(con, ranges, header = NULL, select = NULL) {
    tabix <- Rsamtools::scanTabix(con, param = ranges)
    if (is.null(header)) {
        header <- read_bincov_header(con)
    }
    coltypes <- list("character" = 1L, "integer" = c(2L, 3L), "double" = seq.int(4, length(header)))
    if (!is.null(select)) {
        header <- header[select]
    } else {
        select <- seq_along(header)
    }

    mat <- mapply(\(w, x, y, z) parse_bincov_lines(w, x, y, z, header, select, coltypes),
            tabix, as.character(GenomicRanges::seqnames(ranges)),
            S4Vectors::start(ranges), S4Vectors::end(ranges), SIMPLIFY = FALSE,
            USE.NAMES = FALSE)
    mat <- rbindlist(mat, use.names = TRUE)

    chr <- NULL
    start <- NULL
    end <- NULL
    setkey(mat, chr, start, end)

    mat
}

#' Normalize coverage values from the bincov matrix.
#'
#' Normalization is done by dividing each sample's coverage values by the
#' genome-wide median coverage of the sample.
#'
#' @param x `data.table` Binned coverage values.
#' @param medians `double` Genome-wide median coverage for each sample in `x`.
#' @returns `data.table` Normalized coverage values.
#' @export
normalize_bincov <- function(x, medians) {
    samples <- colnames(x)[!colnames(x) %in% c("chr", "start", "end", "qchr", "qstart", "qend")]
    if (!all(samples %in% names(medians))) {
        stop("all samples in the bincov matrix must have a median coverage")
    }

    medians <- medians[samples]
    tmp <- data.table::copy(x)
    tmp[, names(.SD) := mapply(`/`, .SD, medians, SIMPLIFY = FALSE, USE.NAMES = FALSE), .SDcols = names(medians)]

    tmp
}

parse_bincov_lines <- function(lines, qchr, qstart, qend, header, select, coltypes) {
    if (length(lines) == 0) {
        return(NULL)
    }

    start <- NULL
    d <- fread(text = lines, sep = "\t", header = FALSE, colClasses = coltypes, select = select)
    colnames(d) <- header
    # bincov coordinates are 0-start, exclusive-end
    d[, start := start + 1L]
    d[, c("qchr", "qstart", "qend") := list(qchr, qstart, qend)]
    d
}
