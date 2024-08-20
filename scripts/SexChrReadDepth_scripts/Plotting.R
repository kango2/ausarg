library(karyoploteR)
library(ggplot2)
library(tidyverse)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
chr_length <- args[1]
chr_filter <- args[2]
RD_male_input <- args[3]
RD_female_input <- args[4]
blast_input <- args[5]
density_input <- args[6]
output <- args[7]


custom.genome <- toGRanges(chr_length)
# manually set giestain in all regions of genome to gneg because it's white in karyoploter ideogram, this is purely for plotting nice white bars and not a real feature of genome
filter <- scan(chr_filter, character())

#Read depth
## male
df <- read.table(RD_male_input, header=FALSE)
df_median <- median(df$V4)
df$V5 <- df$V4 / df_median * 100
df <- df %>% select(V1,V2,V3,V5)
names(df) <- c("Chromosome", "Start", "End", "AverageDepth")
malegr <- toGRanges(data.frame(mutate(df, AverageDepth = case_when(AverageDepth > 150 ~ 150, TRUE ~ AverageDepth)) %>%
                                   select(chr = Chromosome, 
                                          start = Start, 
                                          end = End, 
                                          y = AverageDepth)
))

## female
df2 <- read.table(RD_female_input, header=FALSE)
df2_median <- median(df2$V4)
df2$V5 <- df2$V4 / df2_median * 100
df2 <- df2 %>% select(V1,V2,V3,V5)
names(df2) <- c("Chromosome", "Start", "End", "AverageDepth")
femalegr <- toGRanges(data.frame(mutate(df2, AverageDepth = case_when(AverageDepth > 150 ~ 150, TRUE ~ AverageDepth)) %>%
                                   select(chr = Chromosome, 
                                          start = Start, 
                                          end = End, 
                                          y = AverageDepth)
))

# Y kmer blast
kmer_df <- read.table(blast_input, header=FALSE, sep = "\t")
colnames(kmer_df) <- c("query", "V1", "pid", "length", "mismatch", "gapopen", "qstart", "qend", "start", "end", "evalue", "bitscore")
kmer_df$ystart <- 0
kmer_df$yend <- 1
kmer_df <- kmer_df %>%
        mutate(start = ifelse(end < start, end, start),
               end = ifelse(end < start, start, end))


blastgr <- toGRanges(data.frame(kmer_df %>%
                                   select(chr = V1,
                                          x0 = start,
                                          x1 = end,
                                          y0 = ystart,
                                          y1 = yend)
))

# Repeat density
density <- read.table(density_input, header=FALSE, sep = "\t")
colnames(density) <- c("chromosome", "start", "end", "fstart", "fend", "winsize", "density")
palette <- gradient_n_pal(
  colours = c("#1C5796", "#397AA8", "#579EB9", "#89C0C4", "#BCE2CF", "#FFFFE0", "#FAD4AC", "#F0A882", "#E47961", "#C65154", "#A52747"),
  values = c(0, 0.5, 1)
)
density$colour <- palette(density$density)



# PLOTTING
pp <- getDefaultPlotParams(plot.type=2)
pp$data1max <- 150
pp$data1min <- 0
pp$data2height <- 50
pp$data1inmargin <- 10
pp$data2inmargin <- 30
pp$leftmargin <- 0.1

pdf(width=32, height=64, file=output)
kp <- plotKaryotype(genome=custom.genome, chromosomes=filter, plot.type=2, plot.params=pp, cytobands = custom.genome)
kpDataBackground(kp)
kpAxis(kp, numticks=16, labels="", tick.len=0.004 * max(width(kp$plot.region)), col="gray")
kpAxis(kp, numticks=4, tick.len=0.008 * max(width(kp$plot.region)), col="black")
kpPoints(kp, data = malegr,
        ymin = 0, ymax = 150,
        col = "#6495ED",data.panel = 1, cex=0.1)
kpPlotLoess(kp, data = malegr, col = "#6495ED", span=0.1, conf.interval=NA)
kpPoints(kp, data = femalegr,
        ymin = 0, ymax = 150,
        col = "#ff0000",data.panel = 1, cex=0.1)
kpPlotLoess(kp, data = femalegr, col = "#ff0000", span=0.1, conf.interval=NA)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,
                 minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
kpRect(kp, data = blastgr, col = "black", border=NA, data.panel = "ideogram")
kpPlotRegions(kp, data=density, col=density$colour, border=, avoid.overlapping=FALSE, data.panel = 2)
kpAddLabels(kp, labels="Read depth\n(% of median)", data.panel = 1, label.margin = 0.03)
kpAddLabels(kp, labels="Repeat Density", data.panel = 2, label.margin = 0.03)
dev.off()

