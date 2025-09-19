#' Perform k-mer histogram analysis and generate plots
#'
#' @param projectdir Base project directory path
#' @param plotdir Directory name for plots (relative to projectdir)
#' @param tablesdir Directory name for tables (relative to projectdir)
#' @param rawdataqcdir Directory name for raw data QC files (relative to projectdir)
#' @param plot_filename Name for the output PDF plot (default: "kmerhisto.pdf")
#' @param table_filename Name for the output table (default: "kmerpeaks.tsv")
#' @param exclude_ont Logical, whether to exclude ONT data from analysis (default: TRUE)
#' @param kcount_min Minimum k-count for peak detection (default: 10)
#' @param kcount_max Maximum k-count for peak detection (default: 100)
#' @param plot_width Width of output plot in inches (default: 8)
#' @param plot_height Height of output plot in inches (default: 6)
#' @param plot.out Logical, whether to display plot in plot window (default: TRUE)
#' @return invisible(NULL)
#' @example
#' kmer_histograms(projectdir="/g/data/xl04/genomeprojects/TILRUE")

kmer_histograms <- function(projectdir,
                            plotdir = "figures",
                            tablesdir = "tables", 
                            rawdataqcdir = "rawdataqc",
                            plot_filename = "kmerhisto.pdf",
                            table_filename = "kmerpeaks.tsv",
                            exclude_ont = TRUE,
                            kcount_min = 10,
                            kcount_max = 100,
                            plot_width = 8,
                            plot_height = 6,
                            plot.out = TRUE) {
  
  # Load required libraries
  require(dplyr)
  require(readr)
  require(ggplot2)
  require(stringr)
  require(scales)
  
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
  
  # Create output directories if they don't exist
  plot_path <- file.path(projectdir, plotdir)
  table_path <- file.path(projectdir, tablesdir)
  
  if (!dir.exists(plot_path)) {
    dir.create(plot_path, recursive = TRUE)
    message("Created plot directory: ", plot_path)
  }
  
  if (!dir.exists(table_path)) {
    dir.create(table_path, recursive = TRUE) 
    message("Created table directory: ", table_path)
  }
  
  # Read k-mer histogram files
  setwd(rawdata_path)
  files <- list.files(".", pattern = "\\.histo$")
  
  if (length(files) == 0) {
    stop("No .histo files found in: ", rawdata_path)
  }
  
  message("Found ", length(files), " histogram files")
  
  # Process histogram files
  khisto <- list()
  for (i in 1:length(files)) {
    f <- files[i]
    
    # Parse filename to extract tech and klen
    x <- str_split(f, "\\.", simplify = TRUE)
    if (ncol(x) < 3) {
      warning("Unexpected filename format: ", f, ". Skipping.")
      next
    }
    
    # Read histogram data
    tryCatch({
      a <- read_delim(files[i], delim = " ", col_names = FALSE, show_col_types = FALSE)
      colnames(a) <- c("kcounts", "kfrequency")
      a$tech <- x[2]
      a$klen <- x[3]
      khisto[[i]] <- a
    }, error = function(e) {
      warning("Failed to read file: ", f, ". Error: ", e$message)
    })
  }
  
  # Remove NULL entries from failed reads
  khisto <- khisto[!sapply(khisto, is.null)]
  
  # Combine all histogram data
  khisto <- bind_rows(khisto)
  
  if (nrow(khisto) == 0) {
    stop("No valid histogram data found")
  }
  
  # Filter out ONT data if requested (recommended due to analysis issues)
  if (exclude_ont) {
    khisto <- khisto %>% filter(tolower(tech) != "ont")
    message("Excluded ONT data from analysis")
  }
  
  # Find peaks (homozygous peaks)
  peaks <- khisto %>% 
    filter(kcounts < kcount_max & kcounts > kcount_min) %>% 
    group_by(tech, klen) %>% 
    filter(kfrequency == max(kfrequency)) %>%
    ungroup()
  
  message("Detected ", nrow(peaks), " k-mer peaks")
  
  # Create subtitle for plot
  unique_techs <- paste(unique(khisto$tech), collapse = ", ")
  unique_klens <- paste(unique(khisto$klen), collapse = ", ")
  sub_title <- paste0("A plot of k-mer frequency histograms for technologies (", unique_techs,
                     ")\nand k-mer lengths (", unique_klens, 
                     ")\nRed lines indicate detected homozygous peaks. \nBoxed numbers indicate k-mer count at the peaks.")
  
  # Create k-mer plot
  kmerplot <- khisto %>% 
    mutate(kcounts = case_when(kcounts > kcount_max ~ kcount_max, TRUE ~ kcounts)) %>% 
    group_by(tech, klen, kcounts) %>% 
    summarise(kfrequency = sum(kfrequency), .groups = "drop") %>%
    filter(
      (kcounts > 1 & kcounts <= kcount_max & kfrequency < 5e7 & tolower(tech) == "illumina") | 
        (kcounts > 2 & kcounts <= kcount_max & tolower(tech) == "pb") | 
        (kcounts > 3 & kcounts <= kcount_max & kfrequency < 5e7 & tolower(tech) == "ont")
    ) %>%
    ggplot(aes(x = kcounts, y = kfrequency)) + 
    geom_area(stat = "identity", fill = "lightblue", alpha = 0.7) + 
    geom_vline(data = peaks, aes(xintercept = kcounts), color = "red", linetype = "dashed") +
    geom_label(data = peaks, aes(x = kcounts, y = kfrequency * 1.5, label = kcounts), 
               fill = "white", alpha = 0.8) + 
    theme_bw() +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6, sep = "")) +
    facet_grid(tech ~ klen, scales = "free_y") +
    labs(x = "K-mer counts", y = "Frequency", 
         title = "K-mer Frequency Histograms",
         subtitle = sub_title) +
    theme(text = element_text(size = 16), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 12))
  
  # Save the k-mer plot
  plot_output_path <- file.path(plot_path, plot_filename)
  ggsave(plot_output_path, kmerplot, width = plot_width, height = plot_height)
  
  # Display plot in plot window if requested
  if (plot.out) {
    print(kmerplot)
  }
  
  # Prepare and save peak data table
  peaks_table <- peaks %>% 
    select(Technology = tech, `K-mer Length` = klen, `K-mer Peak` = kcounts, Frequency = kfrequency) %>% 
    arrange(`K-mer Length`, Technology)
  
  table_output_path <- file.path(table_path, table_filename)
  write_delim(peaks_table, table_output_path, delim = "\t")
  
  # Output summary to screen
  cat("\n=== K-mer Analysis Complete ===\n")
  cat("Generated outputs:\n")
  cat("- Figure:", plot_output_path, "\n")
  cat("- Table: ", table_output_path, "\n")
  cat("- Analyzed", nrow(peaks_table), "technology/k-length combinations\n")
  cat("- Peak k-counts range:", min(peaks_table$`K-mer Peak`), "to", max(peaks_table$`K-mer Peak`), "\n")
  
  return(invisible(NULL))
}
