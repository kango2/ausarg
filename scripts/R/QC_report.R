#' Perform raw data QC analysis and generate plots and tables
#'
#' @param projectdir Base project directory path
#' @param plotdir Directory name for plots (relative to projectdir)
#' @param tablesdir Directory name for tables (relative to projectdir)
#' @param rawdataqcdir Directory name for raw data QC files (relative to projectdir)
#' @param technology Technology to use for genome size calculation (e.g., "illumina", "PB", "ont")
#' @param kmer_table_filename Name of the k-mer peaks table file (default: "kmerpeaks.tsv")
#' @param plot_width Width of output plots in inches (default: 8)
#' @param plot_height Height of output plots in inches (default: 5)
#' @param plot.out Logical, whether to display plots in plot window (default: TRUE)
#' @return invisible(NULL)
#' @example
#' rawdata_qc_analysis(projectdir="/g/data/xl04/genomeprojects/TILRUE", technology="illumina")

rawdata_qc_analysis <- function(projectdir,
                                plotdir = "figures",
                                tablesdir = "tables", 
                                rawdataqcdir = "rawdataqc",
                                technology = "illumina",
                                kmer_table_filename = "kmerpeaks.tsv",
                                plot_width = 8,
                                plot_height = 5,
                                plot.out = TRUE) {
  
  # Load required libraries
  require(dplyr)
  require(readr)
  require(ggplot2)
  require(stringr)
  require(tidyr)
  
  # Store original working directory
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  
  # Validate directories exist
  if (!dir.exists(projectdir)) {
    stop("Project directory does not exist: ", projectdir)
  }
  
  rawdata_path <- file.path(projectdir, rawdataqcdir)
  if (!dir.exists(rawdata_path)) {
    stop("Raw data QC directory does not exist: ", rawdata_path)
  }
  
  plot_path <- file.path(projectdir, plotdir)
  table_path <- file.path(projectdir, tablesdir)
  
  # Create output directories if they don't exist
  if (!dir.exists(plot_path)) {
    dir.create(plot_path, recursive = TRUE)
    message("Created plot directory: ", plot_path)
  }
  
  if (!dir.exists(table_path)) {
    dir.create(table_path, recursive = TRUE) 
    message("Created table directory: ", table_path)
  }
  
  # Read k-mer peaks table for genome size calculation
  kmer_table_path <- file.path(table_path, kmer_table_filename)
  if (!file.exists(kmer_table_path)) {
    warning("K-mer peaks table not found: ", kmer_table_path, ". Genome size will not be calculated.")
    kmer_peak <- NA
  } else {
    kmer_data <- read_delim(kmer_table_path, delim = "\t", show_col_types = FALSE)
    kmer_match <- kmer_data %>% filter(tolower(Technology) == tolower(technology))
    if (nrow(kmer_match) == 0) {
      warning("Technology '", technology, "' not found in k-mer peaks table. Genome size will not be calculated.")
      kmer_peak <- NA
    } else {
      kmer_peak <- kmer_match$`K-mer Peak`[1]  # Take first match if multiple
      message("Using k-mer peak of ", kmer_peak, " for technology '", technology, "'")
    }
  }
  
  ########################################################################################################
  ################Raw data QC metrics
  setwd(rawdata_path)
  
  # Read quality frequency files
  files <- list.files(".", pattern = "\\_quality_freq.csv$")
  
  if (length(files) == 0) {
    stop("No _quality_freq.csv files found in: ", rawdata_path)
  }
  
  message("Found ", length(files), " quality frequency files")
  
  qvlen <- list()
  for (i in 1:length(files)) {
    tryCatch({
      x <- read_delim(files[i], delim = ",", show_col_types = FALSE)
      qvlen[[i]] <- x 
    }, error = function(e) {
      warning("Failed to read file: ", files[i], ". Error: ", e$message)
    })
  }
  
  # Remove NULL entries and combine
  qvlen <- qvlen[!sapply(qvlen, is.null)]
  qvlen <- bind_rows(qvlen)
  
  if (nrow(qvlen) == 0) {
    stop("No valid quality frequency data found")
  }
  
  qvlen <- qvlen %>% 
    separate(sample, into = c("samplename", "platform", "flowcell"), sep = "\\.", remove = FALSE)
  
  ### Contour plot for QV vs Read Lengths vs Read numbers
  qvcountourreads <- qvlen %>% 
    group_by(platform, read_length, qv) %>% 
    summarise(read_numbers = sum(read_numbers), .groups = "drop") %>%
    group_by(platform) %>%
    mutate(read_length = read_length + 50, readfraction = read_numbers/sum(read_numbers)) %>% 
    ggplot(aes(x=read_length, y = qv, z = readfraction)) + 
    geom_contour_filled(alpha=1) +
    facet_wrap(~platform, scales = "free") +
    theme_bw() + 
    guides(fill = guide_legend(title = "Read Fraction")) +
    ylab("Average Read QV") + xlab("Read Length (bp)") +
    theme(axis.text.y = element_text(size = 16), 
          axis.text.x = element_text(hjust=1, angle=45, size = 16), 
          legend.position = "bottom")
  
  ## Contour plot for QV vs read length vs base counts
  qvcontourbases <- qvlen %>% 
    group_by(platform, read_length, qv) %>% 
    summarise(basenumbers = sum(read_numbers*(read_length+50)), .groups = "drop") %>%
    group_by(platform) %>%
    mutate(read_length = read_length + 50, basefraction = basenumbers/sum(basenumbers)) %>%
    ggplot(aes(x=read_length, y = qv, z = basefraction)) + 
    geom_contour_filled(alpha=1) +
    facet_wrap(~platform, scales = "free") +
    theme_bw() + 
    guides(fill = guide_legend(title = "Base Fraction")) +
    ylab("Average Read QV") + xlab("Read Length (bp)") +
    theme(axis.text.y = element_text(size = 16), 
          axis.text.x = element_text(hjust=1, angle=45, size = 16), 
          legend.position = "bottom")
  
  ## Plot QV readnumbers
  ## For ONT, only few reads have >30QV. messes up the plot so we mutate >30 QV to 30 QV for ONT
  qvlenplot <- qvlen %>% 
    group_by(platform, read_length, qv) %>% 
    summarise(read_numbers = sum(read_numbers), .groups = "drop") %>% 
    mutate(read_length = read_length + 50) %>% 
    mutate(basecounts = read_length * read_numbers) %>% 
    mutate(qv = case_when(qv > 30 & platform == "ont" ~ 30, TRUE ~ qv)) %>%
    ggplot(aes(x=qv, y=read_numbers)) + 
    geom_bar(stat = "identity") +
    xlab("Read average QV") + ylab("No. of reads") +
    facet_wrap(~platform, scales = "free") +
    theme_bw()
  
  # Display plots if requested
  if (plot.out) {
    print(qvcontourbases)
    print(qvcountourreads)
    print(qvlenplot)
  }
  
  #### Plot QV figures for ONT and PacBio data
  setwd(plot_path)
  ggsave("qvvslength.basescountscontour.pdf", qvcontourbases, width = plot_width, height = plot_height)
  ggsave("qvvslength.readcountscontour.pdf", qvcountourreads, width = plot_width, height = plot_height)
  ggsave("qvvsreadcounts.histogram.pdf", qvlenplot, width = plot_width, height = plot_height)
  
  ### Read in the stats.csv file to get the number of reads etc and generate supplementary and main report table
  setwd(rawdata_path)
  
  files <- list.files(".", "*_stats.csv", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No *_stats.csv files found in: ", rawdata_path)
  }
  
  message("Found ", length(files), " stats files")
  
  datametrics <- list()
  for (i in 1:length(files)) {
    tryCatch({
      x <- read_delim(files[i], delim = ",", show_col_types = FALSE)
      datametrics[[i]] <- x  
    }, error = function(e) {
      warning("Failed to read file: ", files[i], ". Error: ", e$message)
    })
  }
  
  # Remove NULL entries and combine
  datametrics <- datametrics[!sapply(datametrics, is.null)]
  datametrics <- bind_rows(datametrics)
  
  if (nrow(datametrics) == 0) {
    stop("No valid stats data found")
  }
  
  datametrics <- datametrics %>% 
    separate(sample, c("samplename", "tech", "runid"), sep = "\\.", remove = FALSE)
  
  # Calculate genome sequence metrics with median read length
  genomeseqmetrics <- datametrics %>% 
    group_by(tech) %>% 
    summarise(
      total_bases = sum(total_bases), 
      total_reads = sum(total_reads),
      .groups = "drop"
    ) %>% 
    mutate(
      average_read_length = total_bases / total_reads,
      median_read_length = NA  # Placeholder - would need individual read lengths for true median
    )
  
  # Add genome size calculation if k-mer peak is available
  estimated_genome_size <- NA_real_
  if (!is.na(kmer_peak)) {
    tech_bases <- genomeseqmetrics %>% 
      filter(tolower(tech) == tolower(technology)) %>% 
      pull(total_bases)
    
    if (length(tech_bases) > 0) {
      estimated_genome_size <- tech_bases[1] / kmer_peak
      genomeseqmetrics <- genomeseqmetrics %>%
        mutate(estimated_genome_size = case_when(
          tolower(tech) == tolower(technology) ~ estimated_genome_size,
          TRUE ~ NA_real_
        ))
      message("Estimated genome size: ", round(estimated_genome_size, 0), " bp using ", technology, " data")
    } else {
      warning("Technology '", technology, "' not found in stats data for genome size calculation")
      genomeseqmetrics <- genomeseqmetrics %>% mutate(estimated_genome_size = NA_real_)
    }
  } else {
    genomeseqmetrics <- genomeseqmetrics %>% mutate(estimated_genome_size = NA_real_)
  }
  
  ### Write sequence data metrics for supplementary table and main reporting
  setwd(table_path)
  write_delim(datametrics, "sequencedatametrics_supplementary.tsv", delim = "\t")
  write_delim(genomeseqmetrics, "sequencedatametrics.tsv", delim = "\t")
  
  # Output summary to screen
  cat("\n=== Raw Data QC Analysis Complete ===\n")
  cat("Generated outputs:\n")
  cat("- Plots: qvvslength.basescountscontour.pdf, qvvslength.readcountscontour.pdf, qvvsreadcounts.histogram.pdf\n")
  cat("- Tables: sequencedatametrics_supplementary.tsv, sequencedatametrics.tsv\n")
  cat("- Technologies analyzed:", paste(unique(genomeseqmetrics$tech), collapse = ", "), "\n")
  if (!is.na(kmer_peak)) {
    cat("- Used k-mer peak of", kmer_peak, "for genome size estimation\n")
  }
  
  # Print genome size estimate
  if (!is.na(estimated_genome_size)) {
    cat("\n=== GENOME SIZE ESTIMATE ===\n")
    cat("Estimated genome size:", format(round(estimated_genome_size, 0), big.mark = ","), "bp\n")
    cat("Based on technology:", technology, "\n")
    cat("K-mer peak:", kmer_peak, "\n")
  } else {
    cat("\nGenome size could not be estimated (missing k-mer peak data or technology data)\n")
  }
  
  return(invisible(NULL))
}

# Example usage:
# rawdata_qc_analysis(projectdir="/g/data/xl04/genomeprojects/TILRUE", technology="illumina")
# 
# # With custom parameters:
# rawdata_qc_analysis(
#   projectdir = "/path/to/my/project",
#   technology = "PB",
#   plot_width = 10,
#   plot_height = 6
# )
