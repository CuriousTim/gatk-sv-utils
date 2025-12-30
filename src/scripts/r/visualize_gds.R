# Visualize VCF matched to genomic disorder regions
#
# Usage:
# Rscript visualize_gd.R <samples> <gd_regions> <sd_regions> <bincov> <medians> <outdir>
#
# <samples>
# Tab-separated file of VCF variants matched to GD regions
# 1. VCF ID
# 2. GD ID
# 3. sample ID
# No header
#
# <gd_regions>
# The GD regions
#
# <sd_regions>
# BEDX file of hg38 segmental duplication regions.
#
# <bincov>
# Binned coverage matrix.
#
# <medians>
# Median coverages file.
#
# <outdir>
# Output directory.

argv <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(gorilla))

dir.create(argv[[6]])

samples <- fread(argv[[1]], header = FALSE, sep = "\t", colClasses = rep("character", 3), col.names = c("vid", "gdid", "sid"))
samples <- samples[, list(samples = paste0(sid, collapse = ",")), by = c("vid", "gdid")]
setkey(samples, "gdid")
gds <- expand_gdtable(read_gdtable(argv[[2]]), 0.5)
setkey(gds, "GD_ID")

goodies <- gds[samples, nomatch = NULL]

segdups <- read_segdups(argv[[3]])
medians <- read_median_coverages(argv[[5]])
bincov_handle <- bincov_file(argv[[4]], medians)

for (i in seq_len(nrow(goodies))) {
    gd <- goodies[i, ]
    mat <- query(bincov_handle, gd$chr, gd$qstart, gd$qend)
    carriers <- strsplit(gd$samples, ",", fixed = TRUE)[[1]]
    plotter <- gdplotter(gd, mat, segdups, carriers)
    for (j in seq_along(carriers)) {
        if (!carriers[[j]] %in% bincov_handle$header) {
            stop(sprintf("%s is not a sample in bincov at '%s'", carriers[[j]], argv[[4]]))
        }
        plot_name <- sprintf("%s~~%s~~%s.jpg", gd$vid, gd$GD_ID, carriers[[j]])
        plot_path <- file.path(argv[[6]], plot_name)
        jpeg(plot_path, res = 100, width = 960, height = 540)
        plot(plotter, carrier = j)
        dev.off()
    }
}
