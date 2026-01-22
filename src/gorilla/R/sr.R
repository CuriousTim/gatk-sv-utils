#' Create a new `sr_file` object.
#'
#' @param x Path to the SR matrix file.
#' @returns A `sr_file` object.
#' @export
sr_file <- function(x) {
    new_sr_file(x)
}

#' @rdname query
#' @export
query.sr_file <- function(x, contig, start, end) {
    coltypes <- list(character(), integer(), character(), integer(), character())
    results <- tabix(x$handle, contig, start, end, coltypes)

    names(results) <- c(
        "contig",
        "pos",
        "side",
        "count",
        "sample_id"
    )

    tmp <- data.table::as.data.table(results)
    sample_id <- NULL
    pos <- NULL
    data.table::setkey(tmp, sample_id, pos)

    structure(
        list(
            mat = tmp,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "sr_mat"
    )
}

#' @rdname subset_samples
#' @export
subset_samples.sr_mat <- function(x, samples) {
    stopifnot(is.character(samples) && length(samples) > 0)

    samples <- unique(samples)
    mat <- x$mat[samples, mult = "all", nomatch = NULL]
    data.table::setkey(mat, sample_id, pos)
    x$mat <- mat

    x
}

new_sr_file <- function(x) {
    stopifnot(is_string(x))

    handle <- Rsamtools::TabixFile(x)

    structure(list(handle = handle), class = "sr_file")
}
