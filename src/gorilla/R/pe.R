#' @rdname svevidencefiles
#' @export
pe_file <- function(path, cachedir = NULL) {
    new_pe_file(path, cachedir)
}

#' @rdname query
#' @export
query.pe_file <- function(x, contig, start, end) {
    results <- query(x$handle, contig, start, end)

    if (file.size(results) == 0) {
        tmp <- data.table::data.table(
            rcontig = character(),
            rstart = integer(),
            rstrand = character(),
            mcontig = character(),
            mstart = integer(),
            mstrand = character(),
            sample_id = character(),
            key = c("sample_id", "rstart")
        )
    } else {
        tmp <- data.table::fread(
            results,
            header = FALSE,
            colClasses = c(
                "character",
                "integer",
                "character",
                "character",
                "integer",
                "character",
                "character"
            ),
            col.names = c(
                "rcontig",
                "rstart",
                "rstrand",
                "mcontig",
                "mstart",
                "mstrand",
                "sample_id"
            ),
            key = c("sample_id", "rstart")
        )
    }

    file.remove(results)

    structure(
        list(
            mat = tmp,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "pe_mat"
    )
}

#' @export
c.pe_mat <- function(...) {
    dots <- list(...)

    if (length(dots) == 0) {
        return(NULL)
    }

    if (length(dots) < 1) {
        return(dots[[1]])
    }

    first <- dots[[1]]
    first_region <- first$region
    mats <- vector("list", length(dots))
    mats[[1]] <- first$mat
    for (i in seq(2, length(dots))) {
        if (!inherits(dots[[i]], "pe_mat")) {
            stop("all objects must be `pe_mat`")
        }

        if (!identical(dots[[i]]$region, first_region)) {
            stop("only objects from the same region can be merged")
        }

        mats[[i]] <- dots[[i]]$mat
    }

    mat <- data.table::rbindlist(mats)
    sample_id <- NULL
    rstart <- NULL
    data.table::setkey(mat, sample_id, rstart)

    structure(
        list(
            mat = mat,
            region = first_region
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

new_pe_file <- function(path, cachedir) {
    stopifnot(is_string(path))
    stopifnot(is.null(cachedir) || is_string(cachedir))

    handle <- new_tabix_handle(path, cachedir)

    structure(list(handle = handle), class = "pe_file")
}
