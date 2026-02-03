#' @rdname svevidencefiles
#' @export
pe_file <- function(path) {
    new_pe_file(path)
}

#' @rdname query
#' @export
query.pe_file <- function(x, contig, start, end) {
    coltypes <- c(list(
        character(),
        integer(),
        character(),
        character(),
        integer(),
        character(),
        character()
    ))
    results <- tabix(x$handle, contig, start, end, coltypes)

    names(results) <- c(
        "rcontig",
        "rstart",
        "rstrand",
        "mcontig",
        "mstart",
        "mstrand",
        "sample_id"
    )

    tmp <- data.table::as.data.table(results)
    sample_id <- NULL
    rstart <- NULL
    data.table::setkey(tmp, sample_id, rstart)

    structure(
        list(
            mat = tmp,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "pe_mat"
    )
}

#' @rdname subset_samples
#' @export
subset_samples.pe_mat <- function(x, samples) {
    stopifnot(is.character(samples))

    samples <- unique(samples)
    mat <- x$mat[samples, mult = "all", nomatch = NULL]

    sample_id <- NULL
    rstart <- NULL
    data.table::setkey(mat, sample_id, rstart)
    x$mat <- mat

    x
}

new_pe_file <- function(path) {
    stopifnot(is_string(path))

    handle <- Rsamtools::TabixFile(path)

    structure(list(handle = handle), class = "pe_file")
}
