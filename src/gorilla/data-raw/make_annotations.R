# Retrieve gene and segmental duplication annotations from UCSC

library(GenomicRanges)
library(UCSC.utils)

genes <- fetch_UCSC_track_data("hg38", "knownGene")
genes <- as.data.table(genes)
genes <- genes[chrom %in% paste0("chr", c(1:22, "X", "Y")), ]
# only keep genes in the basic set
genes <- genes[grepl("\\<basic\\>", tier), ]
# pick highest rank transcript for each gene
genes <- genes[, .SD[which.min(rank)], by = "geneName"]
genes[, color := fcase(transcriptClass == "pseudo", "#ff33ff",
                       transcriptClass == "nonCoding", "#006400",
                       transcriptClass == "coding", "#0c0c78")]
setkey(genes, chrom, chromStart)

# thickStart: where to start drawing the thick block (start codon)
# thickEnd: where to stop drawing the thick block (stop codon)
# chromStarts: start positions of blocks as offset from chromStart (exons)
# blockSizes: sizes of the blocks
genes_gr <- GRanges(
    genes$chrom,
    IRanges(genes$chromStart + 1L, genes$chromEnd),
    strand = genes$strand,
    thickStart = genes$thickStart + 1L,
    thickEnd = genes$thickEnd,
    chromStarts = genes$chromStarts,
    blockSizes = genes$blockSizes,
    geneName = genes$geneName,
    color = genes$color
)

# fracMatch: fraction of matching bases (similarity)
segdups <- fetch_UCSC_track_data("hg38", "genomicSuperDups")
segdups <- segdups[segdups$chrom %in% paste0("chr", c(1:22, "X", "Y")), ]
segdups_gr <- GRanges(
    segdups$chrom,
    IRanges(segdups$chromStart + 1L, segdups$chromEnd),
    strand = segdups$strand,
    fracMatch = segdups$fracMatch
)

usethis::use_data(genes_gr, segdups_gr, internal = TRUE, overwrite = TRUE)
