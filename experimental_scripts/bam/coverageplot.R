# Load the coverage.txt file into R
coverage <- read.table("/g/data/xl04/ka6418/sequence_alignment/bam_outputs/tiliqua_rugosa/pacbio/merged/assembly.haplotype1/coverage.txt", header=FALSE, col.names=c("chrom", "pos", "coverage"))

# Calculate the coverage distribution using the table function
coverage_dist <- table(coverage$coverage)

# Load the ggplot2 package
library(ggplot2)

# Create a histogram of the coverage distribution using ggplot2
ggplot(data.frame(cov=names(coverage_dist), count=as.numeric(coverage_dist)), aes(x=cov, y=count)) +
  geom_bar(stat="identity") +
  labs(title="Coverage Distribution", x="Coverage", y="Count")

# Save the plot as a PNG file
ggsave("coverage_plot.png")
