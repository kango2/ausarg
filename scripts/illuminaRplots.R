library(tidyverse)
library(ggplot2)

generate_qc_plot_combined_side_by_side <- function(csv_paths, identifiers) {
  # Initialize an empty data frame for plotting data
  plot_data <- tibble(Read_Position = integer(), Mean_QV = numeric(), File = character(), Pair = character())
  
  # Loop through each CSV file
  for (i in seq_along(csv_paths)) {
    qc_metrics_df <- read_csv(csv_paths[i])
    
    # Splitting 'Mean_QV_at_read_position' and converting into numeric vectors
    qc_metrics_df <- qc_metrics_df %>%
      separate(Mean_QV_at_read_position, into = c("Mean_QV_at_read_position_1", "Mean_QV_at_read_position_2"), sep = ";")
    
    mean_qv_1 <- map_dbl(str_split(qc_metrics_df$Mean_QV_at_read_position_1[1], ",")[[1]], as.numeric)
    mean_qv_2 <- map_dbl(str_split(qc_metrics_df$Mean_QV_at_read_position_2[1], ",")[[1]], as.numeric)
    
    # Append to plot_data with file and pair identifiers
    plot_data <- plot_data %>%
      bind_rows(
        tibble(
          Read_Position = 1:length(mean_qv_1),
          Mean_QV = mean_qv_1,
          File = identifiers[i],
          Pair = "R1"
        ),
        tibble(
          Read_Position = 1:length(mean_qv_2),
          Mean_QV = mean_qv_2,
          File = identifiers[i],
          Pair = "R2"
        )
      )
  }
  
  # Plot with data from all files, side by side for R1 and R2
  p <- ggplot(plot_data, aes(x = Read_Position, y = Mean_QV)) +
    geom_line(aes(color = File)) +
    facet_wrap(~Pair, scales = "free", ncol = 2) +
    scale_color_viridis_d(begin = 0.2, end = 0.8, direction = 1, name = "File") +
    labs(x = 'Read Position', y = 'Mean QV', title = 'Mean Quality Value (QV) at Read Position for R1 and R2') +
    theme_minimal() +
    ylim(0, 50)
  
  print(p)
}

# Example usage, assuming you have multiple CSVs and corresponding identifiers
csv_paths <- c("/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc/350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R1_001.fastq_350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R2_001.fastq_QC.csv", "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc/350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L002_R1_001.fastq_350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L002_R2_001.fastq_QC.csv") # Update these paths
identifiers <- c("File1", "File2") # Update these identifiers

generate_qc_plot_combined_side_by_side(csv_paths, identifiers)



library(tidyverse)

parse_qv_data <- function(qv_data_str) {
  # Convert the comma-separated string to numeric values
  qv_frequencies <- as.numeric(strsplit(qv_data_str, ",")[[1]])
  
  # Create a data frame: 'Quality_Value' as a sequence, and 'Frequency'
  qv_df <- tibble(Quality_Value = seq_along(qv_frequencies), Frequency = qv_frequencies)
  return(qv_df)
}

generate_stacked_qv_histogram <- function(csv_paths) {
  all_qv_data <- map_dfr(csv_paths, function(csv_path) {
    # Read the CSV and separate 'Mean_QV_per_read' into 'R1_QV' and 'R2_QV'
    df <- read_csv(csv_path) %>%
      mutate(R1_QV = str_split(Mean_QV_per_read, ";", simplify = TRUE)[,1],
             R2_QV = str_split(Mean_QV_per_read, ";", simplify = TRUE)[,2]) %>%
      slice(1)
    
    # Parse QV data strings into data frames for plotting
    r1_qv_df <- parse_qv_data(df$R1_QV[1])
    r2_qv_df <- parse_qv_data(df$R2_QV[1])
    
    # Add a column to distinguish between R1 and R2, and another for the file source
    r1_qv_df <- r1_qv_df %>%
      mutate(Read = "R1", Source = basename(csv_path))
    r2_qv_df <- r2_qv_df %>%
      mutate(Read = "R2", Source = basename(csv_path))
    
    # Combine R1 and R2 data
    bind_rows(r1_qv_df, r2_qv_df)
  }, .id = "Source")
  
  # Plotting
  ggplot(all_qv_data, aes(x = Quality_Value, y = Frequency, fill = Read)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Source, scales = "free_x") +
    labs(title = "QV Frequencies by Read and Source", x = "Quality Value", y = "Frequency") +
    theme_minimal()
}


generate_stacked_qv_histogram(c("/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc/350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R1_001.fastq_350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R2_001.fastq_QC.csv", "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc/350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L002_R1_001.fastq_350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L002_R2_001.fastq_QC.csv"))
