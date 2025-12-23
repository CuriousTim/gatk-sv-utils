# Visualize Genomic Disorder Regions
#
# Usage:
# Rscript visualize_gd.R <gd_regions> <sd_regions> <bincov> <medians> <outdir>
#
# <gd_regions>
# Tab-separated file of genomic disorder regions and samples to visualize.
# 1. chr: contig of the region
# 2. start_GRCh38: start of the region
# 3. end_GRCh38: end of the region
# 4. GD_ID: ID of the region
# 5. svtype: either 'DEL' or 'DUP'
# 6. samples: comma-separated samples to visualize
#    The file must have column headers as given.
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

dir.create(argv[[5]])

gds <- fread(argv[[1]], header = TRUE, sep = "\t", colClasses = c("character", "integer", "integer", "character", "character", "character"))
pad_size <- ceiling((gds$end_GRCh38 - gds$start_GRCh38 + 1L) * 0.5)
gds[, c("qstart", "qend") := list(pmax(1L, start_GRCh38 - pad_size), pmin(536870912L, end_GRCh38 + pad_size))]

medians <- read_median_coverages(argv[[4]])
bincov_handle <- bincov_file(argv[[3]], medians)

for (i in seq_len(nrow(gds))) {
    gd <- gds[i, ]
    mat <- query(bincov_handle, gd$chr, gd$qstart, gd$qend)
    carriers <- strsplit(gd$samples, ",", fixed = TRUE)[[1]]
    plotter <- gdplotter(gd, mat, segdups, carriers)
    for (j in carriers) {
        plot_name <- sprintf("%s_%d-%d_%s_%s_%s.jpg",
                             gd$chr, gd$start_GRCh38, gd$end_GRCh38, gd$GD_ID, gd$svtype, carriers[[j]])
        plot_path <- file.path(argv[[5]], plot_name)
        jpeg(plot_path, res = 100, width = 960, height = 540)
        plot(plotter, carrier = i)
        dev.off()
    }
}
