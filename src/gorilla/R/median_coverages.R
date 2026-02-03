#' Median sequencing coverage file
#'
#' Read a file with genome-wide median sequencing coverage for samples.
#'
#' The file should contain two lines. The first line must be tab-separated
#' sample IDs and the second line must be tab-separated genome-wide median
#' sequencing coverage for each sample.
#'
#' @param path Path to the file.
#' @returns A named vector in which the names are the sample IDs and the values
#'   are the median coverages.
#' @export
#'
#' @examples
#' medians_path <- system.file(
#'     "extdata",
#'     "example_medianCov.txt",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' rd_medians <- read_median_coverages(medians_path)
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
