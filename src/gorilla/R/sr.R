#' @rdname svevidencefiles
#' @export
sr_file <- function(path, cachedir = NULL) {
    new_sr_file(path, cachedir)
}

#' @rdname query
#' @export
query.sr_file <- function(x, contig, start, end) {
    results <- query(x$handle, contig, start, end)

    if (file.size(results) == 0) {
        tmp <- data.table::data.table(
            contig = character(),
            pos = integer(),
            side = character(),
            count = integer(),
            sample_id = character(),
            key = c("sample_id", "pos")
        )
    } else {
        tmp <- data.table::fread(
            results,
            header = FALSE,
            colClasses = c(
                "character",
                "integer",
                "character",
                "integer",
                "character"
            ),
            col.names = c(
                "contig",
                "pos",
                "side",
                "count",
                "sample_id"
            ),
            key = c("sample_id", "pos")
        )
    }

    file.remove(results)

    structure(
        list(
            mat = tmp,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "sr_mat"
    )
}

#' @export
c.sr_mat <- function(...) {
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
        if (!inherits(dots[[i]], "sr_mat")) {
            stop("all objects must be `sr_mat`")
        }

        if (!identical(dots[[i]]$region, first_region)) {
            stop("only objects from the same region can be merged")
        }

        mats[[i]] <- dots[[i]]$mat
    }

    mat <- data.table::rbindlist(mats)
    sample_id <- NULL
    pos <- NULL
    data.table::setkey(mat, sample_id, pos)

    structure(
        list(
            mat = mat,
            region = first_region
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

new_sr_file <- function(path, cachedir) {
    stopifnot(is_string(path))
    stopifnot(is.null(cachedir) || is_string(cachedir))

    handle <- new_tabix_handle(path, cachedir)

    structure(list(handle = handle), class = "sr_file")
}
