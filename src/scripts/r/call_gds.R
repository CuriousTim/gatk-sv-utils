# Call Genomic Disorder CNVs
#
# Usage:
# Rscript call_gds.R [options] <gd_regions> <sd_regions> <bincov> <medians> \
#   <samples> <sex_ploidy> <outdir>
# gd_regions      genomic disorder regions to visualize
# sd_regions      segmental duplications
# bincov          binned coverage matrix
# medians         coverage medians
# samples         list of samples to check for CNVs
# sex_ploidy      sex choromosome ploidy table, 0 for unknown, 1 for male, 2 for female
# outdir          output directory
#
# options
# --min-shift             minimum amount a sample's read depth ratio must shifted from 1 to
#                         considered a CNV carrier
# --pad                   fraction by which the genomic disorder region should be expanded
#                         for plotting
# --max-calls-per-sample  maximum number of calls per sample
# --outliers              samples that had more than the max number of calls

validate_args <- function(x) {
    if (!is.finite(x$min_shift) || x$min_shift < 0) {
        stop("minimum shift must be a non-negative number")
    }

    if (!is.finite(x$pad) || x$pad < 0) {
        stop("padding must be a non-negative number")
    }

    if (x$max_calls_per_sample <= 0) {
        stop("max calls per sample must be a positive number")
    }

    x
}

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    pos_args <- vector("list", 7)
    opts <- list(min_shift = 0.3, pad = 0.5,
                 max_calls_per_sample = 3, outliers = NULL)
    i <- 1
    j <- 1
    repeat {
        if (i > length(args)) {
            break
        }

        if (args[[i]] == "--min-shift") {
            i <- i + 1
            opts$min_shift <- as.double(args[[i]])
        } else if (args[[i]] == "--pad") {
            i <- i + 1
            opts$pad <- as.double(args[[i]])
        } else if (args[[i]] == "--max-calls-per-sample") {
            i <- i + 1
            opts$max_calls_per_sample <- trunc(as.double(args[[i]]))
        } else if (args[[i]] == "--outliers") {
            i <- i + 1
            opts$outliers <- args[[i]]
        } else {
            pos_args[[j]] <- args[[i]]
            j <- j + 1
        }

        i <- i + 1
    }

    if (any(sapply(pos_args, is.null))) {
        stop("incorrect number of arguments", call. = FALSE)
    }
    names(pos_args) <- c("gd_regions", "sd_regions", "bincov", "medians",
                         "samples", "sex_ploidy", "outdir")

    validate_args(append(pos_args, opts))
}

# Main ------------------------------------------------------------------------

argv <- parse_args()

suppressPackageStartupMessages(library(gorilla))

forage(argv$gd_regions, argv$sd_regions, argv$bincov, argv$medians, argv$samples, argv$sex_ploidy, argv$outdir, argv$min_shift, argv$pad, argv$max_calls_per_sample, argv$outliers)
