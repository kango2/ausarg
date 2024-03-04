# Load required libraries
library(tidyverse)
library(optparse)

# Create an option parser
option_parser <- OptionParser(usage = "usage: %prog [options]",
                              description = "Script for processing and saving PDF files.")

# Define the options
option_parser <- add_option(option_parser, c("-t", "--top_folder"), 
                            type = "character", 
                            default = "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/longread_qc",
                            help = "Top folder path for reading files")

option_parser <- add_option(option_parser, c("-o", "--output_folder"), 
                            type = "character", 
                            default = "/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental",
                            help = "Output folder path for saving the PDF")

option_parser <- add_option(option_parser, c("-p", "--pacbio_only"), 
                            action = "store_true", 
                            default = FALSE,
                            help = "Process only the PacBio part")

option_parser <- add_option(option_parser, c("-u", "--ont_only"), 
                            action = "store_true", 
                            default = FALSE,
                            help = "Process only the ONT part")

# Parse the options
options <- parse_args(option_parser)

# Assign parsed options to variables
top_folder <- options$top_folder
output_folder <- options$output_folder
pacbio_only <- options$pacbio_only
ont_only <- options$ont_only

process_data <- function(file_path, dataset_name) {
  data <- read.csv(file_path)
  data$Dataset <- dataset_name
  data$Number_of_Bases <- data$Read_Length * data$Summed_Read_Numbers
  return(data)
}


if (pacbio_only || (!pacbio_only && !ont_only)) {

    # Detect all files with _ont_length_freq in their names within the top folder
    file_paths <- list.files(path = top_folder, pattern = "_pacbio_length_freq\\.csv$", full.names = TRUE)

    # Extract dataset names from file paths
    dataset_names <- str_replace(basename(file_paths), "_pacbio_length_freq\\.csv$", "")

    # Apply the function to each dataset and combine them
    combined_data <- map2_df(file_paths, dataset_names, process_data)

    # Create the stacked histogram
    plot <- ggplot(combined_data, aes(x = Read_Length, y = Summed_Read_Numbers, fill = Dataset)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "PacBio Read Numbers vs Read Length",
        x = "Read Length", y = "Read Numbers") +
    xlim(0, 50000) +
    scale_fill_viridis_d()

    # Concatenate all dataset names to form the PDF file name
    sample_name <- str_extract(dataset_names[1], "^[^_]+")
    pdf_file_name <- paste(sample_name, collapse = "_")
    pdf_file_name <- paste0(pdf_file_name, "_pacbio_length_freq.pdf")

    # Check if the output folder exists, if not, create it
    if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    }

    # Full path for the PDF file
    pdf_full_path <- file.path(output_folder, pdf_file_name)

    # Save the plot as a PDF
    ggsave(filename = pdf_full_path, plot = plot, device = "pdf", width = 11, height = 8.5) 

        # List all CSV files in the top folder that contain "pacbio_quality_freq"
    file_paths <- list.files(path = top_folder, pattern = "pacbio_quality_freq", full.names = TRUE)

    # Initialize an empty list to store data frames
    dfs <- list()

    # Loop through the CSV files, read the data, and extract the date part as 'Legend'
    for (file_path in file_paths) {
    data <- read.csv(file_path)
    file_name <- basename(file_path)
    legend <- gsub("_pacbio_quality_freq.*", "", file_name)  # Extract date part as legend
    data$Legend <- legend  # Add a 'Legend' column with the extracted date
    dfs[[length(dfs) + 1]] <- data
    }

    # Concatenate the data frames into one
    concatenated_data <- do.call(rbind, dfs)

    # Data preparation
    # Grouping by quality value (QV) and 'Legend', summing up the read numbers for each file
    grouped_data <- concatenated_data %>%
    group_by(Legend, QV) %>%
    summarise(Read_Numbers = sum(Read_Numbers))

    # Creating the stacked bar plot
    p <- ggplot(grouped_data, aes(x = QV, y = Read_Numbers, fill = Legend)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "PacBio Quality Value vs Read Numbers",
        x = "Quality Value (QV)",
        y = "Read Numbers"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    pdf_file_name <- paste(sample_name, collapse = "_")
    pdf_file_name <- paste0(pdf_file_name, "_pacbio_quality_freq.pdf")

    # Check if the output folder exists, if not, create it
    if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    }
    print(p)

    # Display the plot
    # Full path for the PDF file
    pdf_full_path <- file.path(output_folder, pdf_file_name)

    # Save the plot as a PDF
    ggsave(filename = pdf_full_path, plot = p, device = "pdf", width = 11, height = 8.5)


}

if (ont_only || (!pacbio_only && !ont_only)) {

        # Detect all files with _ont_length_freq in their names within the top folder
    file_paths <- list.files(path = top_folder, pattern = "_ont_length_freq\\.csv$", full.names = TRUE)

    # Extract dataset names from file paths
    dataset_names <- str_replace(basename(file_paths), "_ont_length_freq\\.csv$", "")

    # Apply the function to each dataset and combine them
    combined_data <- map2_df(file_paths, dataset_names, process_data)

    # Create the stacked histogram
    plot <- ggplot(combined_data, aes(x = Read_Length, y = Summed_Read_Numbers * Read_Length, fill = Dataset)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "ONT Number of Bases vs Read Length",
        x = "Read Length", y = "Number of Bases") +
    scale_fill_viridis_d()

    # Concatenate all dataset names to form the PDF file name
    sample_name <- str_extract(dataset_names[1], "^[^_]+")
    pdf_file_name <- paste(sample_name, collapse = "_")
    pdf_file_name <- paste0(pdf_file_name, "_ont_length_freq.pdf")

    # Check if the output folder exists, if not, create it
    if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    }

    # Full path for the PDF file
    pdf_full_path <- file.path(output_folder, pdf_file_name)

    # Save the plot as a PDF
    ggsave(filename = pdf_full_path, plot = plot, device = "pdf", width = 11, height = 8.5)

        # List all CSV files in the top folder that contain "pacbio_quality_freq"
    file_paths <- list.files(path = top_folder, pattern = "ont_quality_freq", full.names = TRUE)

    # Initialize an empty list to store data frames
    dfs <- list()

    # Loop through the CSV files, read the data, and extract the date part as 'Legend'
    for (file_path in file_paths) {
    data <- read.csv(file_path)
    file_name <- basename(file_path)
    legend <- gsub("_ont_quality_freq.*", "", file_name)  # Extract date part as legend
    data$Legend <- legend  # Add a 'Legend' column with the extracted date
    dfs[[length(dfs) + 1]] <- data
    }

    # Concatenate the data frames into one
    concatenated_data <- do.call(rbind, dfs)

    # Data preparation
    # Grouping by quality value (QV) and 'Legend', summing up the read numbers for each file
    grouped_data <- concatenated_data %>%
    group_by(Legend, QV) %>%
    summarise(Read_Numbers = sum(Read_Numbers))

    # Creating the stacked bar plot
    p <- ggplot(grouped_data, aes(x = QV, y = Read_Numbers, fill = Legend)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "ONT Quality Value vs Read Numbers",
        x = "Quality Value (QV)",
        y = "Read Numbers"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    pdf_file_name <- paste(sample_name, collapse = "_")
    pdf_file_name <- paste0(pdf_file_name, "_ont_quality_freq.pdf")

    # Check if the output folder exists, if not, create it
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }
    print(p)

    # Display the plot
    # Full path for the PDF file
    pdf_full_path <- file.path(output_folder, pdf_file_name)

    # Save the plot as a PDF
    ggsave(filename = pdf_full_path, plot = p, device = "pdf", width = 11, height = 8.5)

   
}