assert_hg38_contigs <- function(x) {
    if (any(!x %in% HG38_PRIMARY_CONTIGS)) {
        stop("only chromosomes chr1-22,X,Y are permitted")
    }
}

assert_valid_contig_pos <- function(x) {
    if (any(x <= 0L | x > TABIX_MAX_SEQLEN)) {
        stop(sprintf(
            "only chromosome positions between 1 and %d inclusive are permitted",
            TABIX_MAX_SEQLEN
        ))
    }
}

assert_positive_genomic_range <- function(start, end) {
    if (any(end < start)) {
        stop("only positive genomic ranges are permitted")
    }
}

assert_valid_cnvtype <- function(x) {
    if (any(!x %in% c("DUP", "DEL"))) {
        stop("only DUP and DEL CNV types are permitted")
    }
}

assert_valid_nahr <- function(x) {
    if (any(!x %in% c("yes", "no"))) {
        stop("NAHR values must be 'yes' or 'no'")
    }
}

assert_valid_gd_terminal_type <- function(x) {
    if (any(!x %in% c("p", "q", "no"))) {
        stop("GD terminal type must be 'p', 'q', or 'no'")
    }
}
