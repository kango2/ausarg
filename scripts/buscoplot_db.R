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


data <- data %>%
  mutate(
    # Store original counts (for labels)
    abs_single = single_copy_buscos,
    abs_fragmented = fragmented_buscos,
    abs_missing = missing_buscos,
    abs_multi = multi_copy_buscos,
    
    # Create percentage columns (for stacked bar plot)
    single_copy_buscos = abs_single / n_markers * 100,
    fragmented_percentage = abs_fragmented / n_markers * 100,
    missing_percentage = abs_missing / n_markers * 100,
    multi_copy_percentage = abs_multi / n_markers * 100
  ) %>%
  select(organismName, n_markers,
         single_copy_buscos, fragmented_percentage, missing_percentage, multi_copy_percentage,
         abs_single, abs_fragmented, abs_missing, abs_multi)

data <- data %>%
  mutate(
    complete_buscos = abs_single + abs_multi,
    tlabel = paste0(
      "C:", complete_buscos,
      " [S:", abs_single,
      ", D:", abs_multi,
      "], F:", abs_fragmented,
      ", M:", abs_missing
    )
  )


species_filter_path <- "/g/data/xl04/genomeprojects/referencedata/tmp/plot/secondpass/species.txt"
species_to_include <- readLines(species_filter_path)

data <- data %>%
  filter(organismName %in% species_to_include) %>%
  mutate(organismName = factor(organismName, levels = rev(species_to_include)))

# Reshape data to long format
data_long <- data %>%
  pivot_longer(
    cols = c("single_copy_buscos", "multi_copy_percentage", "fragmented_percentage", "missing_percentage"),
    names_to = "category",
    values_to = "percentage"
  ) 
# %>%
#   mutate(category = factor(category, levels = c(
#     "single_copy_buscos",
#     "multi_copy_percentage",
#     "fragmented_percentage",
#     "missing_percentage"
#   )))
# 

# Define pastel colors for the categories
complementary_pastel_colors <- c(
  "multi_copy_percentage" = "#9E9AC8",
  "fragmented_percentage" = "#FC9272",
  "single_copy_buscos" = "#6BAED6",
  "missing_percentage" = "#A1D99B"
)



# Plot the horizontal stacked bar chart
data_long %>%
  ggplot(aes(x = percentage, y = organismName, fill = category)) +
  geom_bar(stat = "identity",width = 0.6) +
  geom_text(data = data, aes(x = 80, y = organismName, label = tlabel),inherit.aes = FALSE, hjust = 0, size = 3.3) +
  labs(
    title = "BUSCO Metrics",
    x = "Percentage",
    y = "Organism",
    fill = "Category"
  ) +
  scale_x_continuous(breaks = seq(80, 100, by = 5)) +  # Adjust x-axis ticks
  coord_cartesian(xlim = c(80, 100)) +  # Set Cartesian coordinates for the x-axis
  scale_fill_manual(values = c("#6BAED6", "#FC9272", "#A1D99B", "#9E9AC8")) +  # Apply pastel color scheme
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("/g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/plots/busco_plot.pdf", width = 10, height = 6)



##ggtree
tree_path <- "/g/data/xl04/genomeprojects/referencedata/tmp/referencespecies.nwk"  
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
  theme_tree2() + xlim_expand(c(0, 20), "Tree") # Adjust theme for tree + data visualization

facet_widths(p, widths = c(1, 1))



