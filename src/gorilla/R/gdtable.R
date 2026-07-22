GDTABLE_REQ_HEADER <- c(
    "chr",
    "start_GRCh38",
    "end_GRCh38",
    "GD_ID",
    "svtype",
    "NAHR",
    "terminal",
    "cluster"
)

#' Read a genomic disorders table file
#'
#' Read a TSV file describing genomic disorder regions.
#'
#' The file must have the following columns:
#' 1. chr: contig of the region
#' 2. start_GRCh38: start of the region
#' 3. end_GRCh38: end of the region
#' 4. GD_ID: ID of the region
#' 5. svtype: either 'DEL' or 'DUP'
#' 6. NAHR: either 'yes' or 'no' indicating NAHR-mediated or not
#' 7. terminal: one of 'p', 'q', or 'no' indicating the p-arm, q-arm, or
#'    non-terminal region
#' 8. cluster: either a cluster ID or an empty string
#' The file must have column headers as given.
#'
#' @param path Path to the file.
#' @returns A [`data.table`][data.table::data.table].
#' @export
#'
#' @examples
#' tmp <- tempfile()
#' fp <- file(tmp, open = "wt")
#' writeLines(
#'     "chr\tstart_GRCh38\tend_GRCh38\tGD_ID\tsvtype\tNAHR\tterminal\tcluster",
#'     fp
#' )
#' writeLines("chr1\t898703\t6229913\tGD_0000\tDUP\tno\tp\t", fp)
#' writeLines("chr15\t31727418\t32153204\tGD_0001\tDEL\tyes\tno\tcluster00", fp)
#' close(fp)
#' gds <- read_gdtable(tmp)
read_gdtable <- function(path) {
    tmp <- data.table::fread(
        path,
        header = TRUE,
        sep = "\t",
        strip.white = FALSE,
        na.strings = "",
        colClasses = c(
            "character",
            "integer",
            "integer",
            "character",
            "character",
            "character",
            "character",
            "character"
        )
    )

    if (!identical(colnames(tmp), GDTABLE_REQ_HEADER)) {
        stop(sprintf(
            "genomic disorder regions table must have header: '%s'",
            paste0(GDTABLE_REQ_HEADER, collapse = ", ")
        ))
    }

    assert_hg38_contigs(tmp$chr)
    assert_valid_contig_pos(tmp$start_GRCh3)
    assert_valid_contig_pos(tmp$end_GRCh38)
    assert_positive_genomic_range(tmp$start_GRCh38, tmp$end_GRCh38)
    assert_valid_cnvtype(tmp$svtype)
    assert_valid_nahr(tmp$NAHR)
    assert_valid_gd_terminal_type(tmp$terminal)

    NAHR <- NULL
    tmp[, NAHR := NAHR == "yes"]

    tmp
}

#' Expand genomic disorder regions
#'
#' Expand the regions in a `data.table` of genomic disorder regions.
#'
#' The start of the expanded region will not be less than 1 and the end of the
#' region will not be greater than 5536870912. The object is modified in place.
#'
#' @param x A `data.table` produced by [read_gdtable()].
#' @param prop The fraction of each region by which to expand it.
#' @returns A `data.table` with columns `qstart` and `qend` with the expanded
#'   start and end, invisibly.
#' @export
#'
#' @examples
#' tmp <- tempfile()
#' fp <- file(tmp, open = "wt")
#' writeLines(
#'     "chr\tstart_GRCh38\tend_GRCh38\tGD_ID\tsvtype\tNAHR\tterminal\tcluster",
#'     fp
#' )
#' writeLines("chr1\t898703\t6229913\tGD_0000\tDUP\tno\tp\t", fp)
#' writeLines("chr15\t31727418\t32153204\tGD_0001\tDEL\tyes\tno\tcluster00", fp)
#' close(fp)
#'
#' gds <- read_gdtable(tmp)
#' expanded_gds <- expand_gdtable(gds, 0.2)
expand_gdtable <- function(x, prop) {
    start_GRCh38 <- NULL
    end_GRCh38 <- NULL

    pad_size <- ceiling((x[["end_GRCh38"]] - x[["start_GRCh38"]] + 1L) * prop)
    x[,
        c("qstart", "qend") := list(
            pmax(1L, start_GRCh38 - pad_size),
            pmin(TABIX_MAX_SEQLEN, end_GRCh38 + pad_size)
        )
    ]

    invisible(x)
}
