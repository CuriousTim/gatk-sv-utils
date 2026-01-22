#' Create a new `svtrio` object.
#'
#' The `svtrio` object represents SV evidence over a genomic region for a
#' child, father, mother trio of samples.
#'
#' @param evidence A `sveidence` object.
#' @param child The child sample.
#' @param father The father sample.
#' @param mother The mother sample.
#' @returns A `svtrio` object.
svtrio <- function(evidence, child, father, mother) {
    new_svtrio(evidence, child, father, mother)
}

new_svtrio <- function(evidence, child, father, mother) {
    evidence <- subset_samples(evidence, c(child, father, mother))
    trio <- list(child = child, father = father, mother = mother)

    structure(list(evidence = evidence, trio = trio), class = "svtrio")
}
