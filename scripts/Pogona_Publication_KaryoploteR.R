myLibPath <- "/g/data/xl04/ka6418/temp/Rlibraries"
.libPaths(c(myLibPath, .libPaths()))
.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
.libPaths(c("/g/data/te53/software/Rpackages/are-4.2.2g/", .libPaths()))


library(karyoploteR)
library(tidyverse)
library(optparse)
library(tools) # For file_path_sans_ext to extract base name

# Define options for file paths and other parameters
option_list <- list(
  make_option(c("--genome_csv"), type = "character", default = NULL,
              help = "Path to the genome CSV file", metavar = "character"),
  make_option(c("--telomeres_csv"), type = "character", default = NULL,
              help = "Path to the telomeres CSV file", metavar = "character"),
  make_option(c("--illumina_bed"), type = "character", default = NULL,
              help = "Path to the Illumina depth BED file", metavar = "character"),
  make_option(c("--ont_bed"), type = "character", default = NULL,
              help = "Path to the ONT depth BED file", metavar = "character"),
  make_option(c("--hifi_bed"), type = "character", default = NULL,
              help = "Path to the HiFi depth BED file", metavar = "character"),
  make_option(c("--gc_csv"), type = "character", default = NULL,
              help = "Path to the GC content CSV file", metavar = "character"),
  make_option(c("--nregions_csv"), type = "character", default = NULL,
              help = "Path to the N regions CSV file", metavar = "character"),
  make_option(c("--min_seq_len"), type = "numeric", default = 4e6,
              help = "Minimum sequence length", metavar = "numeric"),
  make_option(c("--output_folder"), type = "character", default = ".",
              help = "Output folder for the PDF", metavar = "character")
)

# Parse the options
opt <- parse_args(OptionParser(option_list = option_list))

# Extract the base name of the genome CSV file and create the output PDF path
genome_base_name <- file_path_sans_ext(basename(opt$genome_csv))
output_pdf_path <- file.path(opt$output_folder, paste0(genome_base_name, ".pdf"))

# Load data using file paths provided as arguments
mygenome <- read_delim(opt$genome_csv, delim = ",")

telomeresdf <- read_delim(opt$telomeres_csv, delim = ",")

column_names <- c("Chromosome", "Start", "End", "AverageDepth")
ilmnrd <- read_delim(opt$illumina_bed, delim = " ", col_names = column_names)

ontrd <- read_delim(opt$ont_bed, delim = " ", col_names = column_names)

hifird <- read_delim(opt$hifi_bed, delim = " ", col_names = column_names)

gc <- read_delim(opt$gc_csv, delim = ",")

column_names <- c("chromosome", "start", "end", "gaplen")
Nregions <- read_delim(opt$nregions_csv, delim = ",", col_names = column_names)

# Process data
mygenomegr <- makeGRangesFromDataFrame(data.frame(
  mutate(mygenome, start = 1) %>%
    filter(Length >= opt$min_seq_len) %>%
    arrange(desc(Length)) %>%
    select(chr = `Sequence ID`, start, end = Length)
))

telomeresgr <- toGRanges(data.frame(select(telomeresdf, chr = Sequence_ID, start = Start, end = End)))

ilmnrdgr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 100 ~ 100, TRUE ~ AverageDepth)) %>%
                                   select(chr = Chromosome, start = Start, end = End, y = AverageDepth)))

ontgr <- toGRanges(data.frame(mutate(ontrd, AverageDepth = case_when(AverageDepth > 100 ~ 100, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, start = Start, end = End, y = AverageDepth)))

hifigr <- toGRanges(data.frame(mutate(hifird, AverageDepth = case_when(AverageDepth > 60 ~ 60, TRUE ~ AverageDepth)) %>%
                                 select(chr = Chromosome, start = Start, end = End, y = AverageDepth)))

gcgr <- toGRanges(data.frame(gc %>% 
                               mutate(y = `GC Count` * 100 / 10000, 
                                      y = case_when(y > 60 ~ 60, TRUE ~ y), 
                                      y = case_when(y < 30 ~ 30, TRUE ~ y)) %>%
                               select(chr = `Sequence Header`, start = `Position Start`, end = `Position End`, y)))

Nregions_gr <- toGRanges(data.frame(Nregions) %>% select(chr = chromosome, start = start, end = end))

# Create PDF output
pdf(output_pdf_path, width = 12, height = 8)

# Visualization
pp <- getDefaultPlotParams(plot.type = 1)
pp$data1inmargin <- 20
pp$data1outmargin <- 120
pp$leftmargin <- 0.15

kp <- plotKaryotype(genome = mygenomegr, plot.params = pp)
kpAddBaseNumbers(kp)
kpPlotRegions(kp, data = telomeresgr, col = "orange", r0 = -0.1, r1 = -0.35, clipping = TRUE, avoid.overlapping = FALSE)
kpPlotRegions(kp, data = Nregions_gr, col = "black", r0 = -0.1, r1 = -0.35, clipping = TRUE, avoid.overlapping = FALSE, lwd = 0.1)

kpDataBackground(kp, data.panel = 1, r0 = 0.0, r1 = 0.15)
kpAxis(kp, r0 = 0.0, r1 = 0.15, cex = 0.2, numticks = 2, ymin = 0, ymax = 100, data.panel = 1, mai = c(0, 0, 0, 1))
kpAddLabels(kp, labels = "Illumina", r0 = 0.0, r1 = 0.15, data.panel = 1, cex = 0.2, label.margin = 0.003, side = "right")

kpLines(kp, data = ilmnrdgr, r0 = 0.0, r1 = 0.15, ymin = 0, ymax = 100, col = "#6495ED", lwd = 0.1, data.panel = 1)
kpLines(kp, data = ontgr, r0 = 0.25, r1 = 0.40, ymin = 0, ymax = 100, col = "#FA8072", lwd = 0.1, data.panel = 1)
kpLines(kp, data = hifigr, r0 = 0.50, r1 = 0.65, ymin = 0, ymax = 60, col = "#3CB371", lwd = 0.1, data.panel = 1)
kpLines(kp, data = gcgr, r0 = 0.75, r1 = 0.90, ymin = 30, ymax = 60, col = "#6A4DA8", lwd = 0.1, data.panel = 1)

dev.off()
