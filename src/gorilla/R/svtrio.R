#' Create a new `svtrio` object
#'
#' The `svtrio` object represents SV evidence over a genomic region for a
#' child, father, mother trio of samples.
#'
#' @param evidence A [`svevidence`] object.
#' @param child The child sample ID.
#' @param father The father sample ID.
#' @param mother The mother sample ID.
#' @returns A `svtrio` object.
#' @export
#'
#' @examples
#' pe_path <- system.file(
#'     "extdata",
#'     "example.PE.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' sr_path <- system.file(
#'     "extdata",
#'     "example.SR.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' rd_path <- system.file(
#'     "extdata",
#'     "example.RD.txt.gz",
#'     package = "gorilla",
#'     mustWork = TRUE
#' )
#' medians_path <- system.file(
#'    "extdata",
#'    "example_medianCov.txt",
#'    package = "gorilla",
#'    mustWork = TRUE
#' )
#'
#' pe <- pe_file(pe_path)
#' sr <- sr_file(sr_path)
#' rd_medians <- read_median_coverages(medians_path)
#' rd <- rd_file(rd_path, rd_medians)
#'
#' sv <- svevidence("chr16", 28743149, 28745149, pe, sr, rd, "DUP")
#' trio <- svtrio(sv, "gorilla0000", "gorilla0001", "gorilla0002")
svtrio <- function(evidence, child, father, mother) {
    new_svtrio(evidence, child, father, mother)
}

new_svtrio <- function(evidence, child, father, mother) {
    stopifnot(inherits(evidence, "svevidence"))
    stopifnot(is_string(child))
    stopifnot(is_string(father))
    stopifnot(is_string(mother))

    evidence <- subset_samples(evidence, c(child, father, mother))
    trio <- list(child = child, father = father, mother = mother)

    structure(list(evidence = evidence, trio = trio), class = "svtrio")
}
