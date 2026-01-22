# query a file indexed with tabix for a genomic region
# `coltypes` is a list of column types like `list(double(), integer())`
tabix <- function(handle, contig, start, end, coltypes) {
    tmp <- Rsamtools::scanTabix(
        handle,
        param = GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
    )[[1]]

    if (length(tmp) == 0) {
        return(coltypes)
    }

    scan(
        text = tmp,
        what = coltypes,
        nmax = length(tmp),
        sep = "\t",
        quiet = TRUE
    )
}
