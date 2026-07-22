#' Query an object with a genomic region
#'
#' `query()` retrieves some genomic region from `x`. Usually `x` is an object
#' representing some on-disk file indexed by genomic coordinates.
#'
#' The method for [`rd_file`][rd_file()] retrieves all the bins
#' overlapping a region from the RD matrix and returns an `rd_mat` object.
#'
#' The method for [`sr_file`][sr_file()] retrieves all split read alignments
#' overlapping a region from the SR matrix and returns an `sr_mat` object.
#'
#' The method for [`pe_file`][pe_file()] retrieves all discordant paired-end
#' reads alignments overlapping a region from the PE matrix and returns a
#' `pe_mat` object. Only the upstream read is checked for overlap.
#'
#' @param x An object.
#' @param contig The contig containing the region.
#' @param start The start of the region (1-start).
#' @param end The end of the region (inclusive).
#' @returns An object representing the queried region.
#' @export
#'
#' @examples
#' pe_path <- system.file("extdata", "example.PE.txt.gz", package = "gorilla", mustWork = TRUE)
#' pe <- pe_file(pe_path)
#' mat <- query(pe, "chr16", 28743149, 28745149)
query <- function(x, contig, start, end) {
    UseMethod("query")
}
