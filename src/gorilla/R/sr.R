#' @rdname svevidencefiles
#' @export
sr_file <- function(path) {
    new_sr_file(path)
}

#' @rdname query
#' @export
query.sr_file <- function(x, contig, start, end) {
    coltypes <- list(
        character(),
        integer(),
        character(),
        integer(),
        character()
    )
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

    sample_id <- NULL
    pos <- NULL
    data.table::setkey(mat, sample_id, pos)
    x$mat <- mat

    x
}

new_sr_file <- function(path) {
    stopifnot(is_string(path))

    handle <- Rsamtools::TabixFile(path)

    structure(list(handle = handle), class = "sr_file")
}
