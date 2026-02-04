# Visualize potential de novo SV
#
# The script only works for a single batch and all members of the trio must be
# in the same batch.
#
# Usage:
# Rscript visualize_denovos.R <variants> <pedigree> <pe> <sr> <rd> <median_cov> <outdir> <exclusions>

nzchar2 <- function(x) {
    !is.na(x) && !nzchar(x)
}

argv <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(gorilla))

dir.create(argv[[7]])

variants <- fread(
    argv[[1]],
    header = FALSE,
    sep = "\t",
    colClasses = c(
        "character",
        "integer",
        "integer",
        "integer",
        "character",
        "character",
        "character"
    ),
    col.names = c(
        "chr",
        "start",
        "end",
        "svlen",
        "vid",
        "svtype",
        "sample_id"
    )
)

pedigree <- fread(
    argv[[2]],
    header = FALSE,
    sep = "\t",
    colClasses = c(
        "character",
        "character",
        "character",
        "character",
        "integer",
        "character"
    ),
    col.names = c(
        "family_id",
        "sample_id",
        "paternal_id",
        "maternal_id",
        "sex",
        "phenotype"
    ),
    key = "sample_id"
)

pe <- pe_file(argv[[3]])
sr <- sr_file(argv[[4]])
rd_medians <- read_median_coverages(argv[[6]])
rd <- rd_file(argv[[5]], rd_medians)
exclusions_fp <- file(argv[[8]], open = "wt")
writeLines("sample_id\tvid\texclusion_reason", exclusions_fp)

for (i in seq_len(nrow(variants))) {
    v <- variants[i, ]
    fam <- pedigree[v$sample_id, ]
    if (!nzchar2(fam$paternal_id) || !nzchar2(fam$maternal_id)) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "parents missing from pedigree"),
            exclusions_fp
        )
        next
    }

    if (!v$sample_id %in% names(rd_medians)) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "sample missing from evidence"),
            exclusions_fp
        )
        next
    }

    if (!fam$paternal_id %in% names(rd_medians)) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "father missing from evidence"),
            exclusions_fp
        )
        next
    }

    if (!fam$maternal_id %in% names(rd_medians)) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "mother missing from evidence"),
            exclusions_fp
        )
        next
    }

    sve <- tryCatch(
        svevidence(v$chr, v$start, v$end, pe, sr, rd, v$svtype, pad = 0.3),
        scanTabix_io = function(e) {
            writeLines(
                sprintf("%s\t%s\t%s", v$sample_id, v$vid, "error querying evidence files"),
                exclusions_fp
            )
            NULL
        }
    )
    if (is.null(sve)) {
        next
    }
    trio <- svtrio(sve, fam$sample_id, fam$paternal_id, fam$maternal_id)
    plotter <- denovo_plotter(trio)
    plot_path <- file.path(argv[[7]], sprintf("%s~~%s.png", v$vid, v$sample_id))
    png(plot_path, width = 3840, height = 2160, res = 300)
    plot(plotter)
    dev.off()
}

close(exclusions_fp)
