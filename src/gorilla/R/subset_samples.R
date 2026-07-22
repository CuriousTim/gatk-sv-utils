#' Subset an object to some samples
#'
#' `subset_samples()` subsets `x` to the samples in `samples`.
#'
#' Here, `samples` refers to entities in an experiment rather than the a random
#' subset. Methods should ignore samples in `samples` that are not present in
#' `x`.
#'
#' There are methods for `pe_mat`, `sr_mat`, `rd_mat`, and `svevidence`.
#'
#' @param x Object to subset.
#' @param samples String identifiers for the samples to keep.
#' @returns `x` subset to `samples`.
#' @export
#'
#' @examples
#' pe_path <- system.file("extdata", "example.PE.txt.gz", package = "gorilla", mustWork = TRUE)
#' pe <- pe_file(pe_path)
#' mat <- query(pe, "chr16", 28743149, 28747149)
#' mat_sub <- subset_samples(mat, c("gorilla0001", "gorilla0002"))
subset_samples <- function(x, samples) {
    UseMethod("subset_samples")
}
