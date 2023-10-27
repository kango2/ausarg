#!/usr/bin/env Rscript

library(tidyverse)
library(karyoploteR)
library(optparse)

# Define command-line arguments
option_list <- list(
    make_option(c("-s", "--seqtable"), type="character", default=NULL, 
                help="Path to Sequence Table file", metavar="character"),
    make_option(c("-t", "--telomeres"), type="character", default=NULL, 
                help="Path to Telomeres csv", metavar="character"),
    make_option(c("-r", "--repetitive"), type="character", default=NULL, 
                help="Path to Summary csv from Centromeric repeats", metavar="character"),
    make_option(c("-i", "--ilmnrd"), type="character", default=NULL, 
                help="Path to Illumina read depth csv", metavar="character"),
    make_option(c("-o", "--ontrd"), type="character", default=NULL, 
                help="Path to ONT read depth csv", metavar="character"),
    make_option(c("-p", "--pacbio"), type="character", default=NULL, 
                help="Path to PacBio read depth csv", metavar="character"),
    make_option(c("-g", "--gc"), type="character", default=NULL, 
                help="Path to GC content CSV", metavar="character"),
    make_option(c("-u", "--output"), type="character", default=NULL, 
                help="Path to the output PDF file", metavar="character")
)

# Parse command-line arguments
args <- parse_args(OptionParser(option_list=option_list))

# Determine default name for the output PDF
if (is.null(args$output)) {
    base_name <- sub("_Telomeres.csv", "", basename(args$telomeres))
    args$output <- paste0(base_name, ".pdf")
}



minseqlen <- 20e5

mygenome <- read_delim(args$seqtable, delim = ",")
mygenomegr <- makeGRangesFromDataFrame(data.frame(mutate(mygenome, start = 1) %>% 
                                                  filter(Length >= minseqlen) %>%
                                                  arrange(desc(Length)) %>%
                                                  select(chr = `Sequence ID`, start, end = Length)
                                                ))

telomeresdf <- read_delim(args$telomeres, delim = ",")
telomeresgr <- toGRanges(data.frame(select(telomeresdf, chr=Sequence_ID, start = Start, end = End)))

cendf <- read_delim(args$repetitive, delim = ",")
filtered_cendf <- filter(cendf, !is.na(class))

# Split data based on the class values
list_of_dataframes <- split(filtered_cendf, filtered_cendf$class)

# Convert each dataframe in the list to GRanges and store back in the list
list_of_granges <- lapply(list_of_dataframes, function(df) {
  toGRanges(data.frame(select(df, chr = name, start = start, end = end)))
})

# If you want to store them in separate variables
names(list_of_granges) <- paste0("cengr_", unique(filtered_cendf$class))
list2env(list_of_granges, envir = .GlobalEnv)

cengr <- toGRanges(data.frame(select(cendf, chr=name, start = start, end = end)))

column_names <- c("chromosome", "start", "end", "width")
Nregions <- read_delim("/g/data/xl04/ka6418/bassiana/all_assemblies/Ns/rBasDup_HifiASM_YAHS_Nregions.csv",delim = ",", col_names = column_names)
Nregions_gr <- toGRanges(data.frame(Nregions) %>% select(chr = chromosome, start = start, end = end) )

column_names <- c("Chromosome", "Start", "End", "AverageDepth")
ilmnrd <- read_delim(args$ilmnrd, delim = ",",col_names = column_names)
ilmnrdgr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, 
                                      start = Start, 
                                      end = End, 
                                      y = AverageDepth)
                             ))

column_names <- c("Chromosome", "Start", "End", "AverageDepth")
ontrd <- read_delim(args$ontrd, delim = ",",col_names = column_names)
ontgr <- toGRanges(data.frame(mutate(ontrd, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, 
                                      start = Start, 
                                      end = End, 
                                      y = AverageDepth)
                             ))

column_names <- c("Chromosome", "Start", "End", "AverageDepth")
hifird <- read_delim(args$pacbio, delim = ",",col_names = column_names)
hifigr <- toGRanges(data.frame(mutate(hifird, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, 
                                      start = Start, 
                                      end = End, 
                                      y = AverageDepth)
                             ))

gc <- read_delim(args$gc, delim = ",")
gcgr <- toGRanges(data.frame(gc %>% mutate(y = `GC Count`*100/10000) %>%
                                select(chr = `Sequence Header`, 
                                      start = `Position Start`, 
                                      end = `Position End`, 
                                      y)))


pp <- getDefaultPlotParams(plot.type=1)
pp$data1inmargin <- 20
pp$data1outmargin <- 120

pdf(args$output, width = 11.7, height = 8.3)
kp <- plotKaryotype(genome = mygenomegr,plot.params = pp)

kpAddBaseNumbers(kp)
#kpPlotRegions(kp, data=cengr, col="#FF000080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=cengr_CEN199, col="#FF000080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=cengr_CEN187, col="#80008080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=cengr_MIC187, col="green", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=cengr_MIC199, col="gold", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)

kpPlotRegions(kp, data=telomeresgr, col="orange", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)

kpPlotRegions(kp, data=Nregions_gr, col="black", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE,lwd = 0.1)

kpLines(kp, data = ilmnrdgr,
        r0=0, r1=0.25, ymin = 0, ymax = 200,
        col = "#6495ED", lwd = 0.1) 
kpLines(kp, data = ontgr,
        r0=0.25, r1=0.4, ymin = 0, ymax = 200,
        col = "#FA8072", lwd = 0.1)
kpLines(kp, data = hifigr,
        r0=0.45, r1=0.65, ymin = 0, ymax = 200,
        col = "#3CB371", lwd = 0.1)
kpLines(kp, data = gcgr,
        r0=0.65, r1=0.90, ymin = 40, ymax = 60,
        col = "#6A4DA8", lwd = 0.3)

dev.off()


