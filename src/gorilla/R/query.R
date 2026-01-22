#' Query an object with a genomic region.
#'
#' `query()` retrieves some genomic region from `x`.
#'
#' The method for [`rd_file`][rd_file()] retrieves all the bins
#' overlapping a region from the RD matrix and returns an `rd_mat` object.
#'
#' The method for [`sr_file`][sr_file()] retrieves all split read alignments
#' overlapping a region from the SR matrix and returns an `sr_mat` object.
#'
#' The method for [`pe_file`][pe_file()] retrieves all aberrant paired-end
#' reads alignments overlapping a region from the PE matrix and returns a
#' `pe_mat` object. Only the upstream read is checked for overlap.
#'
#' @param x An object.
#' @param contig The contig containing the region.
#' @param start The start of the region (1-start).
#' @param end The end of the region (inclusive).
#' @returns An object representing the queried region or `NULL` if the region
#'   does not overlap any intervals in `x`.
#' @export
query <- function(x, contig, start, end) {
    UseMethod("query")
}
