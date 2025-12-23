#' Create a new `bincov_file` object.
#'
#' @param x Path to the binned coverage matrix file.
#' @param medians Named numeric vector of median coverage values for every
#'   sample in the binned coverage matrix file.
#' @name bincov_file
#' @export
bincov_file <- function(x, medians) {
    new_bincov_file(x, medians)
}

#' @export
print.bincov_file <- function(x, ...) {
    cat("binned coverage file\n")
    cat("path: ", x$handle$path, "\n", sep = "")
    cat("samples: ", length(x$header) - 3, "\n", sep = "")
}

#' Query an object with a genomic region.
#'
#' `query()` retrieves some genomic region from `x`.
#'
#' The method for [`bincov_file`][bincov_file()] retrieves all the bins
#' overlapping the requested region from the binned coverage matrix and returns
#' a `bincov_mat` object or `NULL` if the matrix does not have any overlappng
#' bins.
#'
#' @param x An object.
#' @param contig The contig containing the region.
#' @param start The start of the region (1-start).
#' @param end The end of the region (inclusive).
#' @returns An object representing the queried region.
#' @export
query <- function(x, contig, start, end) {
    UseMethod("query")
}

#' @rdname query
#' @export
query.bincov_file <- function(x, contig, start, end) {
    tmp <- Rsamtools::scanTabix(
        x$handle,
        param = GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
    )[[1]]
    if (length(tmp) == 0) {
        return(NULL)
    }

    coltypes <- c(
        list(character(), integer(), integer()),
        lapply(seq_len(length(x$header) - 3), double)
    )
    parsed <- scan(
        text = tmp,
        what = coltypes,
        nmax = length(tmp),
        sep = "\t",
        quiet = TRUE
    )
    names(parsed) <- x$header

    # bincov matrix coordinates are 0-start
    ranges <- GenomicRanges::GRanges(
        parsed$chr,
        IRanges::IRanges(parsed$start + 1L, parsed$end)
    )
    mat <- do.call(
        cbind,
        parsed[!names(parsed) %in% c("chr", "start", "end")]
    ) |>
        normalize_bincov(x$medians)

    structure(list(ranges = ranges, rd = mat), class = "bincov_mat")
}

#' @export
print.bincov_mat <- function(x, ...) {
    cat("binned coverage matrix\n")

    first_bin <- utils::head(x$ranges, 1)
    last_bin <- utils::tail(x$ranges, 1)
    cat(sprintf(
        "region: %s:%d-%d\n",
        levels(GenomicRanges::seqnames(first_bin))[[1]],
        S4Vectors::start(first_bin),
        S4Vectors::end(last_bin)
    ))
    cat(sprintf("bins: %d\n", length(x$ranges)))
    cat(sprintf("samples: %d\n", ncol(x$rd)))
}

# for smoothing of the bincov matrix
smooth <- function(x, n) {
    UseMethod("smooth")
}

#' @export
smooth.bincov_mat <- function(x, n) {
    if (nrow(x$rd) >= n) {
        x$rd <- apply(x$rd, 2, \(y) stats::runmed(y, n, na.action = "fail"))
    }

    x
}

# Normalize the raw binned coverage values by dividing by the median coverage
# across the genome.
normalize_bincov <- function(x, medians) {
    scale(x, FALSE, medians)
}

#' @export
`[.bincov_mat` <- function(x, i, j, ..., drop = FALSE) {
    x$ranges <- x$ranges[i]
    x$rd <- x$rd[i, j, drop = drop]

    x
}

new_bincov_file <- function(x, medians) {
    stopifnot(is.character(x) && length(x) == 1)
    handle <- Rsamtools::TabixFile(x)
    header <- read_bincov_header(handle)

    if (!all(header[-(1:3)] %in% names(medians))) {
        stop("every sample must have a median coverage value")
    }

    structure(
        list(
            handle = handle,
            header = header,
            medians = medians[header[-(1:3)]]
        ),
        class = "bincov_file"
    )
}

read_bincov_header <- function(con) {
    header <- Rsamtools::headerTabix(con)[["header"]]
    if (is.null(header) || length(header) == 0) {
        stop("bincov matrix is missing a header")
    }

    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) < 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop("bincov matrix has an invalid header")
    }
    header[1:3] <- c("chr", "start", "end")

    header
}
