#' Read the median coverages file.
#'
#' @param path Path to the file.
#' @returns A named vector in which the names are the sample IDs and the values
#'   are the median coverages.
#' @export
read_median_coverages <- function(path) {
    tmp <- readLines(path, n = 2L, ok = FALSE)
    ids <- strsplit(tmp[[1]], split = "\t", fixed = TRUE)[[1]]
    if (anyDuplicated(ids) != 0) {
        stop("samples in median coverages file must be unique")
    }
    covs <- as.double(strsplit(tmp[[2]], split = "\t", fixed = TRUE)[[1]])

    if (any(covs <= 0)) {
        stop("median coverages must be positive")
    }

    stats::setNames(covs, ids)
}
