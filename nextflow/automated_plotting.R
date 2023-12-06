# Load required libraries
library(tidyverse)

# Specify the top folder path for reading files
top_folder <- "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/longread_qc"

# User-provided output folder path for saving the PDF
output_folder <- "/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental" # User should modify this path


# Function to process each dataset
process_data <- function(file_path, dataset_name) {
  data <- read.csv(file_path)
  data$Dataset <- dataset_name
  data$Number_of_Bases <- data$Read_Length * data$Summed_Read_Numbers
  return(data)
}


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
  xlim(0, 100000) +
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

dev.off()

# User-supplied top folder
top_folder <- "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/longread_qc/"
output_folder <- "/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental"

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


