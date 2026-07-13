# Visualize potential de novo SV
#
# The script only works for a single batch and all members of the trio must be
# in the same batch.
#
# Usage:
# Rscript visualize_denovos.R <variants> <pedigree> <evidence> <batches> <outdir> <exclusions>
#
# <variants>
# TSV of variants to visualize
#
# <pedigree>
# Pedigree in GATK PED format
#
# <evidence>
# TSV with paths to the SV evidence files with columns in the order
# batch ID, PE, SR, RD, median cov
#
# <batches>
# TSV with mapping between batch IDs and sample IDs in the order
# batch ID, sample ID
#
# <outdir>
# Directory to write the plots
#
# <exclusions>
# File to write variants that were skipped

REQ_VARIANT_FIELDS <- c("chr", "start", "end", "svlen", "vid", "svtype", "sample_id")

nzchar2 <- function(x) {
    !is.na(x) && nzchar(x)
}

argv <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(gorilla))

dir.create(argv[[5]])

# read variants to visualize
variants <- fread(argv[[1]], header = TRUE, sep = "\t")

if (!all(REQ_VARIANT_FIELDS %in% colnames(variants))) {
    stop(
        paste0(
            "all of these fields must be in the variant file: ",
            paste0(REQ_VARIANT_FIELDS, collapse = ", ")
        ),
        call. = FALSE
    )
}

# read pedigree
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

# read PE/SR/RD/medians paths
evidence_paths <- fread(
    argv[[3]], header = FALSE, sep = "\t",
    colClasses = c(
        "character", "character", "character", "character", "character"
    ),
    col.names = c(
        "batch_id", "pe_path", "sr_path", "rd_path", "medians_path"
    ),
    key = "batch_id"
)

# sample to batch mapping
samples <- fread(
    argv[[4]], header = FALSE, sep = "\t",
    colClasses = c("character", "character"),
    col.names = c("sample_id", "batch_id"),
    key = "sample_id"
)

exclusions_fp <- file(argv[[6]], open = "wt")
writeLines("sample_id\tvid\texclusion_reason", exclusions_fp)

tmpdir <- tempfile("temp", ".")
if (!dir.create(tmpdir, showWarnings = FALSE)) {
    stop("failed to create temporary directory in working directory")
}

for (i in seq_len(nrow(variants))) {
    v <- variants[i, ]
    fam <- pedigree[v$sample_id, nomatch = NULL]

    if (nrow(fam) == 0) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "sample missing from pedigree"),
            exclusions_fp
        )
        next
    }

    if (!nzchar2(fam$paternal_id) || !nzchar2(fam$maternal_id)) {
        writeLines(
            sprintf("%s\t%s\t%s", v$sample_id, v$vid, "parents missing from pedigree"),
            exclusions_fp
        )
        next
    }

    child_batch <- samples[fam$sample_id, batch_id, nomatch = NULL]
    paternal_batch <- samples[fam$paternal_id, batch_id, nomatch = NULL]
    maternal_batch <- samples[fam$maternal_id, batch_id, nomatch = NULL]

    child_evidence_paths <- evidence_paths[child_batch, nomatch = NULL]
    if (nrow(child_evidence_paths) == 0) {
        stop(sprintf("batch '%s' has no evidence files", child_batch))
    }
    paternal_evidence_paths <- evidence_paths[paternal_batch, nomatch = NULL]
    if (nrow(paternal_evidence_paths) == 0) {
        stop(sprintf("batch '%s' has no evidence files", paternal_batch))
    }
    maternal_evidence_paths <- evidence_paths[maternal_batch, nomatch = NULL]
    if (nrow(maternal_evidence_paths) == 0) {
        stop(sprintf("batch '%s' has no evidence files", maternal_batch))
    }

    child_cachedir <- file.path(tmpdir, child_batch)
    child_pe <- pe_file(child_evidence_paths$pe_path, child_cachedir)
    child_sr <- sr_file(child_evidence_paths$sr_path, child_cachedir)
    child_rd <- rd_file(
        child_evidence_paths$rd_path,
        child_evidence_paths$medians_path,
        child_cachedir
    )
    # INS start to end length is always 1b which is not long enough to
    # visualize anything so we add 1Kb padding.
    child_batch_svevidence <- svevidence(
        v$chr, v$start, v$end, child_pe, child_sr, child_rd, v$svtype,
        pad = if (v$svtype == "INS") 1000 else 0.3,
        sr_pad = if (v$svtype == "INS") 300 else NULL
    )

    child_svevidence <- subset_samples(child_batch_svevidence, fam$sample_id)

    if (paternal_batch == child_batch) {
        paternal_svevidence <- subset_samples(child_batch_svevidence, fam$paternal_id)
    } else {
        paternal_cachedir <- file.path(tmpdir, paternal_batch)
        paternal_pe <- pe_file(paternal_evidence_paths$pe_path, paternal_cachedir)
        paternal_sr <- sr_file(paternal_evidence_paths$sr_path, paternal_cachedir)
        paternal_rd <- rd_file(
            paternal_evidence_paths$rd_path,
            paternal_evidence_paths$medians_path,
            paternal_cachedir
        )
        paternal_batch_svevidence <- svevidence(
            v$chr, v$start, v$end, paternal_pe, paternal_sr, paternal_rd, v$svtype,
            pad = if (v$svtype == "INS") 1000 else 0.3,
            sr_pad = if (v$svtype == "INS") 300 else NULL
        )
        paternal_svevidence <- subset_samples(paternal_batch_svevidence, fam$paternal_id)
    }

    if (maternal_batch == child_batch) {
        maternal_svevidence <- subset_samples(child_batch_svevidence, fam$maternal_id)
    } else if (maternal_batch == paternal_batch) {
        maternal_svevidence <- subset_samples(paternal_batch_svevidence, fam$maternal_id)
    } else {
        maternal_cachedir <- file.path(tmpdir, maternal_batch)
        maternal_pe <- pe_file(maternal_evidence_paths$pe_path, maternal_cachedir)
        maternal_sr <- sr_file(maternal_evidence_paths$sr_path, maternal_cachedir)
        maternal_rd <- rd_file(
            maternal_evidence_paths$rd_path,
            maternal_evidence_paths$medians_path,
            maternal_cachedir
        )
        maternal_batch_svevidence <- svevidence(
            v$chr, v$start, v$end, maternal_pe, maternal_sr, maternal_rd, v$svtype,
            pad = if (v$svtype == "INS") 1000 else 0.3,
            sr_pad = if (v$svtype == "INS") 300 else NULL
        )
        maternal_svevidence <- subset_samples(maternal_batch_svevidence, fam$maternal_id)
    }

    sve <- c(child_svevidence, paternal_svevidence, maternal_svevidence)

    trio <- svtrio(sve, fam$sample_id, fam$paternal_id, fam$maternal_id)
    plotter <- denovo_plotter(trio)
    plot_path <- file.path(
        argv[[5]],
        sprintf("%s~~%s~~%s_%d-%d.png", v$vid, v$sample_id, v$chr, v$start, v$end)
    )
    # numbers based on the relative track heights in denovo_plotter
    height_scale <- (1.2 * plotter$pe_track_scale + 1.8) / 3
    png(plot_path, width = 3840, height = 2160 * height_scale, res = 300)
    plot_call <- call("plot", x = quote(plotter))
    if ("o_ev" %in% colnames(v)) {
        plot_call[["evtype"]] <- v$o_ev
    }
    if ("ac" %in% colnames(v)) {
        plot_call[["ac"]] <- v$ac
    }
    eval(plot_call)
    dev.off()
}

close(exclusions_fp)
