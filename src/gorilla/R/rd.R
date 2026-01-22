#' Create a new `rd_file` object.
#'
#' @param x Path to the RD matrix file.
#' @param medians Named numeric vector of median coverage values of samples in
#'   the RD matrix file. See [read_median_coverages()].
#' @returns A `rd_file` object.
#' @export
rd_file <- function(x, medians) {
    new_rd_file(x, medians)
}

#' @export
print.rd_file <- function(x, ...) {
    cat("RD matrix file\n")
    cat("path: ", x$handle$path, "\n", sep = "")
    cat("samples: ", length(x$header) - 3, "\n", sep = "")
}

#' @rdname query
#' @export
query.rd_file <- function(x, contig, start, end) {
    coltypes <- c(
        list(character(), integer(), integer()),
        lapply(seq_len(length(x$header) - 3), double)
    )
    results <- tabix(x$handle, contig, start, end, coltypes)
    names(results) <- x$header

    mat <- do.call(
        cbind,
        results[!names(results) %in% c("chr", "start", "end")]
    )

    if (nrow(mat) > 0) {
        # RD matrix coordinates are 0-start
        ranges <- GenomicRanges::GRanges(
            results$chr,
            IRanges::IRanges(results$start + 1L, results$end)
        )
        mat <- normalize_bincov(mat, x$medians)
    } else {
        ranges <- GenomicRanges::GRanges()
    }

    structure(
        list(
            ranges = ranges,
            mat = mat,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "rd_mat"
    )
}

#' @rdname subset_samples
#' @export
subset_samples.rd_mat <- function(x, samples) {
    stopifnot(is.character(samples) && length(samples) > 0)

    samples <- unique(samples)
    mat <- x$mat[, samples, drop = FALSE]
    x$mat <- mat

    x
}

#' @export
print.rd_mat <- function(x, ...) {
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
    cat(sprintf("samples: %d\n", ncol(x$mat)))
}

smooth <- function(x, n) {
    UseMethod("smooth")
}

#' @export
smooth.rd_mat <- function(x, n) {
    mat <- x$mat
    if (nrow(x$mat) >= n) {
        mat <- apply(mat, 2, \(y) stats::runmed(y, n, na.action = "fail"))
    }

    structure(
        list(
            ranges = x$ranges,
            mat = mat,
            region = x$region,
            smooth_size = n
        ),
        class = "smooth_rd_mat"
    )
}

new_rd_file <- function(x, medians) {
    stopifnot(is_string(x))
    stopifnot(is.double(medians))

    handle <- Rsamtools::TabixFile(x)
    header <- read_rd_header(handle)
    check_rd_medians(header, medians)

    structure(
        list(
            handle = handle,
            header = header,
            medians = medians[header[-(1:3)]]
        ),
        class = "rd_file"
    )
}

read_rd_header <- function(con) {
    header <- Rsamtools::headerTabix(con)[["header"]]
    if (is.null(header) || length(header) == 0) {
        stop("RD matrix is missing a header")
    }

    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) < 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop("RD matrix has an invalid header")
    }
    header[1:3] <- c("chr", "start", "end")

    header
}

check_rd_medians <- function(header, medians) {
    if (!all(header[-(1:3)] %in% names(medians))) {
        stop("every sample in RD file must have a median coverage value")
    }
}

# Normalize the raw RD values by dividing by the median coverage
# across the genome.
normalize_bincov <- function(x, medians) {
    scale(x, FALSE, medians)
}
