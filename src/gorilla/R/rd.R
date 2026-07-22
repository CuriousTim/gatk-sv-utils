#' @rdname svevidencefiles
#' @export
rd_file <- function(path, medians, cachedir = NULL) {
    new_rd_file(path, medians, cachedir)
}

#' @export
print.rd_file <- function(x, ...) {
    cat("RD matrix file\n")
    cat("path: ", x$handle$path, "\n", sep = "")
}

#' @rdname query
#' @export
query.rd_file <- function(x, contig, start, end) {
    results <- query(x$handle, contig, start, end)

    tmp <- data.table::fread(results, header = TRUE)
    header <- colnames(tmp)
    if (length(header) < 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop("RD matrix has an invalid header")
    }

    check_rd_medians(header, x$medians)

    ccols <- c("chr", "start", "end")
    data.table::setnames(tmp, c("#Chr", "Start", "End"), ccols)
    tmp[, names(.SD) := lapply(.SD, as.integer), .SDcols = c("start", "end")]
    tmp[, names(.SD) := lapply(.SD, as.double), .SDcols = !ccols]

    mat <- as.matrix(tmp[, .SD, .SDcols = !ccols])

    if (nrow(mat) > 0) {
        # RD matrix coordinates are 0-start
        ranges <- GenomicRanges::GRanges(
            tmp$chr,
            IRanges::IRanges(tmp$start + 1L, tmp$end)
        )
        mat <- normalize_rd(mat, x$medians)
    } else {
        ranges <- GenomicRanges::GRanges()
    }

    file.remove(results)

    structure(
        list(
            ranges = ranges,
            mat = mat,
            region = list(contig = contig, start = start, end = end)
        ),
        class = "rd_mat"
    )
}

#' @export
c.rd_mat <- function(...) {
    dots <- list(...)

    if (length(dots) == 0) {
        return(NULL)
    }

    if (length(dots) < 1) {
        return(dots[[1]])
    }

    first <- dots[[1]]
    first_region <- first$region
    first_ranges <- first$ranges
    mats <- vector("list", length(dots))
    mats[[1]] <- first$mat
    for (i in seq(2, length(dots))) {
        if (!inherits(dots[[i]], "rd_mat")) {
            stop("all objects must be `rd_mat`")
        }

        if (!identical(dots[[i]]$region, first_region)) {
            stop("only objects from the same region can be merged")
        }

        ranges <- dots[[i]]$ranges
        if (!identical(ranges, first_ranges)) {
            stop("only objects with the same ranges can be merged")
        }

        mats[[i]] <- dots[[i]]$mat
    }

    mat <- do.call(cbind, mats)
    if (anyDuplicated(colnames(mat)) != 0) {
        stop("column names of all matrices must be unique")
    }

    structure(
        list(
            ranges = first_ranges,
            mat = mat,
            region = first_region
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

#' @importFrom stats median
#' @export
median.rd_mat <- function(x, na.rm = FALSE, region = NULL, ...) {
    if (!is.null(region)) {
        gr <- GenomicRanges::GRanges(
            region$contig,
            IRanges::IRanges(region$start, region$end)
        )
        rows <- suppressWarnings(GenomicRanges::findOverlaps(gr, x$ranges)) |>
            S4Vectors::subjectHits()

        if (length(rows) == 0) {
            return(stats::setNames(rep(NA_real_, ncol(x$mat)), colnames(x$mat)))
        }
    } else {
        rows <- NULL
    }

    matrixStats::colMedians(x$mat, rows = rows, na.rm = na.rm)
}

smooth_rd <- function(x, n) {
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

new_rd_file <- function(x, medians, cachedir) {
    stopifnot(is_string(x))
    stopifnot(is_string(medians))
    stopifnot(is.null(cachedir) || is_string(cachedir))

    if (is.null(cachedir)) {
        cachedir <- tempfile("rd", tempdir(TRUE))
    }

    tryCatch(
        mkdir(cachedir),
        error = function(e) {
            stop("failed to create RD cache directory")
        }
    )

    if (file.exists(medians)) {
        medians_path <- medians
    } else {
        cached_medians <- file.path(cachedir, basename(medians))
        if (!file.exists(cached_medians)) {
            gcs_download_file(medians, cached_medians)
        }
        medians_path <- cached_medians
    }

    handle <- new_tabix_handle(x, cachedir)

    structure(
        list(
            handle = handle,
            medians = read_median_coverages(medians_path)
        ),
        class = "rd_file"
    )
}

check_rd_medians <- function(header, medians) {
    if (!all(header[-(1:3)] %in% names(medians))) {
        stop("every sample in RD file must have a median coverage value")
    }
}

# Normalize the raw RD values by dividing by the median coverage
# across the genome.
normalize_rd <- function(x, medians) {
    # the order of the samples in the medians vector may not match the order of
    # the samples in the binned read-depth matrix
    scale(x, FALSE, medians[colnames(x)])
}
