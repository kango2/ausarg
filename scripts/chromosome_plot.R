#TODO:Add user arguments

library(tidyverse)
library(karyoploteR)

custom.genome <- toGRanges("/Users/u7119530/Downloads/converted_data.txt")
pp <- getDefaultPlotParams(plot.type=1)
pp$data2height <- 150
pp$data1height <- 150


r <- read_delim("Downloads/reformatted_stalk_telomeres_data.txt", delim = "\t")
c <- read_delim("Downloads/trash_tilrug.csv", delim = ",")
c1 <- filter(c, width>10000 & !is.na(class) & class == "CEN171") %>% select(chr = name, start, end)
c2 <- filter(c, width>10000 & !is.na(class) & class == "CEN82") %>% select(chr = name, start, end)

telomeres <- toGRanges(data.frame(filter(r, str_starts(name,"TRF")) %>% select(chr, start, end)))
c1 <- toGRanges(data.frame(c1))
c2 <- toGRanges(data.frame(c2))

dev.off()
kp <- plotKaryotype(genome = custom.genome,plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 10e6)
kpPlotRegions(kp, telomeres, col = "red", r0 = 0, r1 = 0.3)
kpPlotRegions(kp, c1, col = "blue", r0 = 0.33, r1 = 0.63)
kpPlotRegions(kp, c2, col = "green", r0 = 0.67, r1 = 1)


