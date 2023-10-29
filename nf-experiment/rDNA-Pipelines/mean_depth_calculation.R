# Load required packages
library(tidyverse)
library(GenomicRanges)

# Set the directory path
directory <- "/g/data/te53/ontsv/hrp561/cnv"

# List files with .1kbp.bedcov and .1kbp.bedcov.done extensions
covfiles <- list.files(directory, pattern = "\\.1kbp\\.bedcov$")
donefiles <- list.files(directory, pattern = "\\.1kbp\\.bedcov.done$")

# Extract the base file names from the covfiles
fbase <- str_split(covfiles, "\\.1kbp\\.bedcov$", simplify = TRUE)[, 1]

# Extract sample names from the base file names
snames <- str_split(fbase, "_pass$", simplify = TRUE)[, 1]

# Initialize an empty vector to store names of .bedcov files with corresponding .bedcov.done files
matching_bedcov_files <- character()

# Loop through each sample name and find matching .bedcov files with a .bedcov.done file
for (sample in snames) {
  done_filename <- paste0(sample, "_pass.1kbp.bedcov.done")
  if (done_filename %in% donefiles) {
    matching_bedcov_files <- c(matching_bedcov_files, 
                               file.path(directory, paste0(sample, "_pass.1kbp.bedcov")))
  }
}

# Display the names of the matching .bedcov files
cat("Matching .bedcov file names:\n")
cat(matching_bedcov_files, sep = "\n")
matching_bedcov_files

# Count and display the number of matched files
num_files <- length(matching_bedcov_files)
cat("\nTotal number of files:", num_files, "\n")

# Read in centromere coordinates and convert to GRanges format
centromeres <- read_tsv("/g/data/te53/ontsv/hrp561/scripts/chm13-t2t-censat.ucsc.20230611.bed", 
                        col_names = c("chromosome", "start", "end", "name", "score", 
                                      "strand", "thickStart", "thickEnd", "itemRgb"))
centromere_ranges <- GRanges(seqnames = centromeres$chromosome, 
                             ranges = IRanges(start = centromeres$start, end = centromeres$end))

# Read in rDNA coordinates and convert to GRanges format
rDNA <- read_tsv("/g/data/te53/ontsv/hrp561/scripts/chm13-t2t-rdnamodels.ucsc.20230611.bed", 
                 col_names = c("chromosome", "start", "end", "name"))
rDNA_ranges <- GRanges(seqnames = rDNA$chromosome, 
                       ranges = IRanges(start = rDNA$start, end = rDNA$end))

# Define a function to calculate the mean read depth after filtering 
get_mean_depth <- function(file) {
  df <- read_tsv(file, col_names = c("chromosome", "start", "end", "read_depth"))
  
  # Remove chrM, chrY, and chrX
  df <- df %>% filter(!(chromosome %in% c("chrM", "chrY", "chrX")))
  # Remove rows with zero read depth
  df <- df %>% filter(read_depth != 0)
  
  # Convert data frame to a GRanges object
  df_ranges <- GRanges(seqnames = df$chromosome, 
                       ranges = IRanges(start = df$start, end = df$end))
  
  # Identify overlaps with centromeres and rDNA regions
  centromere_overlaps <- findOverlaps(df_ranges, centromere_ranges)
  rDNA_overlaps <- findOverlaps(df_ranges, rDNA_ranges)
  
  # Exclude overlapping regions from the dataset
  overlaps_to_remove <- unique(c(queryHits(centromere_overlaps), queryHits(rDNA_overlaps)))
  df <- df[-overlaps_to_remove, ]
  
  # Calculate the mean read depth and return
  mean_depth <- mean(df$read_depth, na.rm = TRUE) / 1000
  return(data.frame(file_ID = file, mean_read_depth = mean_depth))
}

# Apply the function to the first matching .bedcov file and store the result
results <- bind_rows(lapply(matching_bedcov_files[1], get_mean_depth))
results

# Calculate mean read depth across all matched files
results_list <- list()

for (file in matching_bedcov_files) {
  df <- read_tsv(file, col_names = c("chromosome", "start", "end", "read_depth"))
  
  # Remove chrM, chrY, and chrX
  df <- df %>% filter(!(chromosome %in% c("chrM", "chrY", "chrX")))
  # Remove rows with zero read depth
  df <- df %>% filter(read_depth != 0)
  
  # Convert data frame to a GRanges object
  df_ranges <- GRanges(seqnames = df$chromosome, 
                       ranges = IRanges(start = df$start, end = df$end))
  
  # Identify overlaps with centromeres and rDNA regions
  centromere_overlaps <- findOverlaps(df_ranges, centromere_ranges)
  rDNA_overlaps <- findOverlaps(df_ranges, rDNA_ranges)
  
  # Exclude overlapping regions from the dataset
  overlaps_to_remove <- unique(c(queryHits(centromere_overlaps), queryHits(rDNA_overlaps)))
  df <- df[-overlaps_to_remove, ]
  
  # Calculate and store the mean read depth for each file
  mean_depth <- mean(df$read_depth, na.rm = TRUE) / 1000
  results_list[[file]] <- mean_depth
}

# Combine results into a data frame and display
results_df <- data.frame(file_ID = names(results_list), mean_read_depth = unlist(results_list))
print(results_df)
view(results_df)
