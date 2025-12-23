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

# Add a `qstart` and `qend` columns to the `data.table` with the new
# coordinates of the regions after expanding. Modification is in place.
expand_gdtable <- function(x, prop) {
    start_GRCh38 <- NULL
    end_GRCh38 <- NULL

    pad_size <- ceiling((x$end_GRCh38 - x$start_GRCh38 + 1L) * prop)
    x[,
        c("qstart", "qend") := list(
            pmax(1L, start_GRCh38 - pad_size),
            pmin(TABIX_MAX_SEQLEN, end_GRCh38 + pad_size)
        )
    ]

    invisible(x)
}

# Read the genomic disorders table file.
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
