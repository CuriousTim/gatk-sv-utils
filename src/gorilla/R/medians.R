#' Read coverage medians file into a named vector.
#'
#' @param path `character(1)` Path to the coverage medians file.
#' @returns `double` The coverage medians. The names of the vector are the
#'   sample IDs.
#' @export
read_medians_file <- function(path) {
    stopifnot("`path` should be a character string" = is.character(path) && length(path) == 1)
    lines <- readLines(path, n = 2L, ok = FALSE)
    ids <- strsplit(lines[[1]], split = "\t", fixed = TRUE)[[1]]
    if (anyDuplicated(ids) != 0) {
        stop("duplicate samples in medians file")
    }
    medians <- as.double(strsplit(lines[[2]], split = "\t", fixed = TRUE)[[1]])

    names(medians) <- ids

    medians
}
