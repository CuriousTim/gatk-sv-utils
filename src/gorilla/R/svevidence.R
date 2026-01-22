#' Create a new `svevidence` object.
#'
#' This class represents SV evidence over a genomic region. It is assumed that
#' the evidence sources all have the same set of samples.
#'
#' @param contig The contig containing the region.
#' @param start The start of the region (1-start).
#' @param end The end of the region (inclusive).
#' @param pe A [`pe_file`] object for the discordant read pairs evidence.
#' @param sr A [`sr_file`] object for the split read evidence.
#' @param rd A [`rd_file`] object for the read depth evidence.
#' @param svtype SV type of the variant.
#' @param pad Fraction of region to add as padding. The minimum padded start is
#'   1 and the maximum padded end is 536870912.
#' @returns A `svevidence` object.
svevidence <- function(contig, start, end, pe, sr, rd, svtype, pad = 0) {
    if (pad < 0) {
        stop("`pad` must be non-negative")
    }
    new_svevidence(contig, start, end, pe, sr, rd, svtype, pad)
}

#' @rdname subset_samples
#' @export
subset_samples.svevidence <- function(x, samples) {
    sample_id <- NULL
    rstart <- NULL
    pos <- NULL

    pe <- subset_samples(x$pe, samples)
    sr <- subset_samples(x$sr, samples)
    rd <- subset_samples(x$rd, samples)

    structure(list(
        pe = pe,
        sr = sr,
        rd = rd,
        svtype = x$svtype,
        region = x$region
    ), class = "svevidence")
}

new_svevidence <- function(contig, start, end, pe, sr, rd, svtype, pad) {
    pad_size <- ceiling((end - start + 1) * pad)
    qstart <- max(1, start - pad_size)
    qend <- min(TABIX_MAX_SEQLEN, end + pad_size)
    pe_mat <- query(pe, contig, qstart, qend)
    sr_mat <- query(sr, contig, qstart, qend)
    rd_mat <- query(rd, contig, qstart, qend)
    region <- list(contig = contig, start = start, end = end, qstart = qstart, qend = qend)

    structure(list(pe = pe_mat, sr = sr_mat, rd = rd_mat, svtype = svtype, region = region), class = "svevidence")
}
