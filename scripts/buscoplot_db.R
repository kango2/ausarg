# Load necessary libraries
library(DBI)
library(tidyverse)

# Connect to the SQLite database
db_path <- "/g/data/xl04/genomeprojects/database/ausarg.db"
conn <- dbConnect(RSQLite::SQLite(), db_path)

# Query the database
query <- "
SELECT 
  organismName, 
  single_copy_buscos, 
  fragmented_buscos, 
  missing_buscos, 
  multi_copy_buscos, 
  n_markers 
FROM busco_metrics
"
data <- dbGetQuery(conn, query)

# Close the database connection
dbDisconnect(conn)

# Calculate percentages for stacked bar chart
data <- data %>%
  mutate(
    single_copy_buscos = single_copy_buscos / n_markers * 100,
    fragmented_percentage = fragmented_buscos / n_markers * 100,
    missing_percentage = missing_buscos / n_markers * 100,
    multi_copy_percentage = multi_copy_buscos / n_markers * 100
  ) %>%
  select(organismName, single_copy_buscos, fragmented_percentage, missing_percentage, multi_copy_percentage)

# Reshape data to long format
data_long <- data %>%
  pivot_longer(
    cols = starts_with("single_copy_buscos"):starts_with("multi_copy_percentage"),
    names_to = "category",
    values_to = "percentage"
  )

# Define pastel colors for the categories
complementary_pastel_colors <- c(
  "single_copy_buscos" = "#6BAED6",  # Darker pastel blue
  "fragmented_percentage" = "#FC9272",  # Darker pastel pink
  "missing_percentage" = "#A1D99B",  # Darker pastel green
  "multi_copy_percentage" = "#9E9AC8"  # Darker pastel purple
)

# Plot the horizontal stacked bar chart
ggplot(data_long, aes(x = percentage, y = organismName, fill = category)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Horizontal Stacked Bar Chart of BUSCO Metrics",
    x = "Percentage",
    y = "Organism",
    fill = "Category"
  ) +
  scale_x_continuous(breaks = seq(80, 100, by = 5)) +  # Adjust x-axis ticks
  coord_cartesian(xlim = c(80, 100)) +  # Set Cartesian coordinates for the x-axis
  scale_fill_manual(values = complementary_pastel_colors) +  # Apply pastel color scheme
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




##ggtree
tree_path <- "/g/data/xl04/genomeprojects/referencedata/tmp/referencespecies.nwk"  # Replace with the actual path
tree <- read.tree(tree_path)

common_organisms <- intersect(data_long$organismName, tree$tip.label)

# Filter the data for matching organisms
data_long_filtered <- data_long %>%
  filter(organismName %in% common_organisms)




pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))
pruned_tree$edge.length <- NULL
p <- ggtree(pruned_tree) +
  geom_tiplab(size = 3) 

colnames(data_long_filtered)[colnames(data_long_filtered) == "organismName"] <- "label"

p <- p + 
  geom_facet(
    data = data_long_filtered,
    mapping = aes(x = percentage, y = label, fill = category),
    geom = geom_col,       # Horizontal bar chart
    panel = "BUSCO",
    width = 0.6
  ) +
  scale_fill_manual(values = complementary_pastel_colors) +  # Use custom colors
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") + xlim_expand(c(85, 100), "BUSCO") # Adjust theme for tree + data visualization

facet_widths(p, widths = c(1, 1))



