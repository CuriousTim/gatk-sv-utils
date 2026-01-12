# Prepare files to run ManualCnvRevisionAcrossContigs in GATK-SV to modify VCFs
# with the changes found from genomic disorder review.
#
# Paths to input/output files that need to be read/written are marked with
# 'CHANGEME'.

########################################################
# Directory
########################################################
# CHANGEME: directory to the manual review 2 folder, i.e. review of calls in
# VCF that were not found the GD calling workflow
setwd("/path/to/directory")


########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(stringr)
library(purrr)
library(openxlsx)


########################################################
# STEP 0: UPDATE gd_missing
########################################################

# CHANGEME: file from comparing GD calls with VCF calls
gd_missing <- as.data.frame(fread("/path/to/file.tsv",
                                  col.names = c("SVID", "chr", "start", "end", "type", "GD",
                                                "cat1", "cat2", "cat3", "cat4", "cat5")))

##################################
# Move samples in 00_correct (i.e., true missed GDs) from cat5 to cat1
# CHANGEME: path to directory with true GD plots from manual review 2, i.e.
# review of VCF calls that overlapped a GD region, but did not trigger a call
# in the GD calling workflow
correct_gds <- list.files("./00_correct/")
correct_gds <- data.frame(SVID = gsub("~~.*", "", correct_gds),
                          GD = gsub("~~.*", "", gsub(".*~~GD", "GD", correct_gds)),
                          sample = gsub(".jpg", "", gsub(".*~~", "", correct_gds)))
correct_gds <- correct_gds %>%
                group_by(across(-sample)) %>%
                summarise(samples = paste(unique(sample), collapse = ","),
                          .groups = "drop")
correct_gds <- as.data.frame(correct_gds)

for(i in 1:nrow(correct_gds)) {

  # Get the samples to correct (cat5 --> cat1)
  samples_to_correct <-  unlist(strsplit(correct_gds[i, "samples"], ","))

  # Identify the relevant row & extract samples from it
  row <- which(gd_missing$SVID == correct_gds[i, "SVID"] & gd_missing$GD == correct_gds[i, "GD"])
  cat1_samples <- unlist(strsplit(gd_missing[row, "cat1"], ","))
  cat5_samples <- unlist(strsplit(gd_missing[row, "cat5"], ","))

  # Correct the file
  gd_missing[row, "cat1"] <- paste0(unique(c(cat1_samples, samples_to_correct)), collapse = ",")
  if (length(setdiff(cat5_samples, samples_to_correct)) == 0) {
    gd_missing[row, "cat5"] <- ""
  } else {
    gd_missing[row, "cat5"] <- paste0(unique(setdiff(cat5_samples, samples_to_correct)), collapse = ",")
  }
}
rm(i, samples_to_correct, row, cat1_samples, cat5_samples)


########################################################
# STEP 1: CREATE new-cnv-table
########################################################
# Table of (1) chrom, (2) pos, (3) end, (4) unique ID possibly corresponding to IDs in --gd-table, and (5) comma-delimited list of carrier samples, for new DEL/DUP records

# CHANGEME: genomic disorder regions file
new_cnv_table <- as.data.frame(read.xlsx("/path/to/file.xlsx",
                     cols = c(1:4)))
colnames(new_cnv_table) <- c("chr", "start", "end", "GD")

# CHANGEME: path to directory of true positive GD calls from manual review 1,
# i.e. review of calls from the GD calling workflow
gd_mr1 <- list.files("/path/to/00_correct", recursive = T)
gd_mr1 <- data.frame(file_name = gsub(".*\\/chr", "chr", gd_mr1))
gd_mr1 <- data.frame(GD = gsub("_DEL___.*", "", (gsub("_DUP___.*", "", gsub(".*_GD", "GD", gd_mr1$file_name)))),
                    sample = gsub(".jpg", "", gsub(".*___", "__", gd_mr1$file_name)))
gd_mr1 <- gd_mr1 %>%
  group_by(across(-sample)) %>%
  summarise(sample_mr1 = paste(unique(sample), collapse = ","),
            .groups = "drop")

new_cnv_table <- as.data.frame(left_join(new_cnv_table, gd_mr1, by = "GD"))
rm(gd_mr1)

# Add carriers from manual review 2
colnames(correct_gds)[3] <- "sample_mr2"
new_cnv_table <- as.data.frame(left_join(new_cnv_table, correct_gds[, c(2,3)], by = "GD"))

# Merge samples from manual review 1 and 2
new_cnv_table$samples <- apply(new_cnv_table[, c("sample_mr1", "sample_mr2")], 1, function(x) {paste(na.omit(x), collapse = ",")})

# Delete GDs with no carriers
new_cnv_table <- new_cnv_table[which(new_cnv_table$samples != ""), c(1:4, 7)]



########################################################
# STEP 2: CREATE remove-call-table
########################################################
# Table of (1) variant ID and (2) sample ID to change to hom-ref genotypes

# Create a dataframe
remove_call_table <- data.frame()

##################################
# Samples that are in cat1/2/4 of gd_missing
with_gd <- data.frame()
for(i in 1:nrow(gd_missing)) {

  # Get the SVID
  svid <- gd_missing[i, "SVID"]

  # Identify samples to remove
  samples_with_gd <- unique(c(unlist(strsplit(gd_missing[i, "cat1"], ",")),
                              unlist(strsplit(gd_missing[i, "cat2"], ",")),
                              unlist(strsplit(gd_missing[i, "cat4"], ","))))
  if (length(samples_with_gd) == 0) {next}

  # Add them to the list
  df_temp <- data.frame(SVID = svid, sample = samples_with_gd)
  with_gd <- rbind(with_gd, df_temp)
}
rm(i, svid, samples_with_gd, df_temp)

# Add to remove_call_table
remove_call_table <- rbind(remove_call_table, with_gd)
rm(with_gd)


##################################
# Manual review 2: cat5 samples that are classified as "02_no_cnv"
# These are variants that are not supported by RD.

# Read in samples
# CHANGEME: path to directory with RD plots of calls in the VCF that were
# determined to be a false positive, i.e. the sample has an alternate genotype
# at a CNV site, but it is homref based on read depth
no_cnv <- list.files("./02_no_cnv/")
no_cnv <- data.frame(SVID = gsub("~~.*", "", no_cnv),
                     sample = gsub(".jpg", "", gsub(".*~~", "", no_cnv)))

# Add to remove_call_table
remove_call_table <- rbind(remove_call_table, no_cnv)
rm(no_cnv)



########################################################
# STEP 3: CREATE remove-vids-list
########################################################
# List of variant IDs to remove

# VIDs that are on the exclusion list because they only contain samples with manually validated GDs
# CHANGEME: path to the "remove_vids" file from MakeGDRevisionTable workflow.
remove_vids_list <- as.data.frame(fread("/path/to/gd_remove_vids.txt",
                                        header = F, col.names = "SVID"))

# Variants that after Manual review 2 only have samples in cat 1/2/4
remove_vids_list <- rbind(remove_vids_list, gd_missing[which(gd_missing$cat3 == "" & gd_missing$cat5 == ""), "SVID", drop = F])


########################################################
# STEP 4: SAVE
########################################################

# new-cnv-table - Needs to be saved 1 file/chr
for (c in paste0("chr", c(seq(1:22), "X", "Y"))) {
  new_cnv_table_chr <- new_cnv_table[which(new_cnv_table$chr == c), ]
  # CHANGEME: where to write the per-contig new CNV tables
  fwrite(new_cnv_table_chr, paste0("/path/to/new_cnv_table_", c,"_suffix.txt"),
         col.names = F, row.names = F, sep = " ", quote = F)
}
# CHANGEME: where to write the entire new CNV table
fwrite(new_cnv_table, "/path/to/new_cnv_table_suffix.txt",
       col.names = F, row.names = F, sep = " ", quote = F)

# remove-call-table
# CHANGEME
fwrite(remove_call_table, "/path/to/remove_calls_table_suffix.txt",
       col.names = F, row.names = F, sep = " ", quote = F)

# remove-vids-list
# CHANGEME
fwrite(remove_vids_list, "/path/to/remove_vids_list_suffix.txt",
       col.names = F, row.names = F, sep = " ", quote = F)

# Save a correct version of gd_missing
# CHANGEME
fwrite(gd_missing[!gd_missing$SVID %in% remove_vids_list$SVID, ], "/path/to/gd_missing-correct.txt",
       col.names = T, row.names = F, sep = "\t", quote = F)
