# Retrieve gene and segmental duplication annotations from UCSC

library(data.table)
library(GenomicRanges)
library(UCSC.utils)

PRIMARY_CHROMS <- paste0("chr", c(1:22, "X", "Y"))

genes <- fetch_UCSC_track_data("hg38", "knownGene")
genes <- as.data.table(genes)
genes <- genes[chrom %in% paste0("chr", c(1:22, "X", "Y")), ]
genes <- genes[grepl("\\<MANE_(Select|Plus_Clinical)\\>", tag), ]
# pick highest rank transcript for each gene
genes <- genes[, .SD[which.min(rank)], by = "geneName"]
genes[, color := fcase(transcriptClass == "nonCoding", "#006400",
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
segdups <- segdups[segdups$chrom %in% PRIMARY_CHROMS, ]
segdups_gr <- GRanges(
    segdups$chrom,
    IRanges(segdups$chromStart + 1L, segdups$chromEnd),
    strand = segdups$strand,
    fracMatch = segdups$fracMatch
)

simple_repeats <- fetch_UCSC_track_data("hg38", "simpleRepeat")
simple_repeats <- simple_repeats[simple_repeats$chrom %in% PRIMARY_CHROMS, ]
simple_repeats_gr <- GRanges(
    simple_repeats$chrom,
    IRanges(simple_repeats$chromStart + 1L, simple_repeats$chromEnd)
)

repeat_mask <- fetch_UCSC_track_data("hg38", "rmsk")
repeat_mask <- repeat_mask[repeat_mask$genoName %in% PRIMARY_CHROMS, ]
repeat_mask_gr <- GRanges(
    repeat_mask$genoName,
    IRanges(repeat_mask$genoStart + 1L, repeat_mask$genoEnd)
)

usethis::use_data(
    genes_gr, segdups_gr, simple_repeats_gr, repeat_mask_gr,
    internal = TRUE, overwrite = TRUE, compress = "xz"
)
