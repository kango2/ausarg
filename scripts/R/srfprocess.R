library(tidyverse)
library(tidygraph)
library(GenomicRanges)

#####merge the class together
srfpaf <- read_paf("analysis/SRF/POGVIT.v2.1.asm2srf.genome.paf", tibble = T)
srfpaf <- srfpaf %>% separate(tname, c("sampleid", "repid"), sep = "#", remove = F) %>% mutate(replen = as.numeric(str_extract(repid, "\\d+$")))

classpaf <- read_paf("analysis/SRF/POGVIT.v2.1.enlong1K.paf", tibble = T)
reppairs <- classpaf %>% select(from = qname, to = tname) %>% distinct()
repgraph <- tbl_graph(nodes = NULL, edges = reppairs, directed = FALSE)
repclass <- repgraph %>% mutate(component = group_components()) %>% as_tibble()
repclass <- repclass %>% mutate(repclass = paste("srfclass", component, sep = "-"))

srfpaf <- left_join(srfpaf, repclass, by = c("tname" = "name"))
repclasses <- srfpaf %>% pull(repclass) %>% unique()

##merge consecutive repeats of the same type into a window
repreduced <- NULL
repreduced <- list()
for (i in 1:length(repclasses)){
  repaln <- srfpaf %>% filter(repclass==repclasses[i]) %>% select(chr = qname, start = qstart, end = qend, de)
  repaln <- makeGRangesFromDataFrame(repaln, keep.extra.columns = T)
  ##merge two consecutive alignments as one if there is <=10bp gap between them, this can be increased to 100bp or corresponding to the gap width of scaffolds
  ##merge can be based on length of repeats as well, i.e. if two consecutive repeats
  repregions <- reduce(repaln, min.gapwidth = 10)
  ###get overlapping alignments to collect metrics of how many alignments overlap the region and the mean de (divergence) score for all alignments
  olap <- findOverlapPairs(repregions, repaln)
  olap <- data.frame(olap) %>% group_by(first.seqnames, first.start, first.end, first.width) %>% summarise(alncounts = n(), meande = mean(second.de))
  olap <- ungroup(olap)
  colnames(olap) <- c("chr","start", "end", "width", "alncounts", "meande")
  olap$repid <- repclasses[i]
  tmprid <- srfpaf %>% filter(repclass==repclasses[i]) %>% pull(tname) %>% unique() %>% sort()
  olap$replen <- as.numeric(str_extract(tmprid[1], "\\d+$"))
  olap <- olap %>% mutate(rangeid = paste(chr, paste(start,end,sep = "-"), sep = ":"))
  repreduced[[i]] <- olap
}
repreduced <- bind_rows(repreduced)
print(paste(nrow(repreduced), " ranges before filtering", sep = ""))

##find overlapping alignments of different repeats
repranges <- makeGRangesFromDataFrame(repreduced, keep.extra.columns = T)
repolap <- data.frame(findOverlapPairs(repranges, repranges, minoverlap=10))
repolap <- repolap %>% group_by(first.rangeid, first.repid) %>% mutate(olapcount = n())
repolap <- repolap %>% ungroup()
##if there is only one overlapping repeat alignments, they are good to go as they will be self overlap
selection1 <- repolap %>%
  filter(olapcount == 1) %>%
  select(chr = first.X.seqnames,
         start = first.X.start,
         end = first.X.end,
         width = first.X.width,
         repid = first.repid,
         replen = first.replen,
         alncounts=first.alncounts,
         meande = first.meande,
         rcopy = first.rcopy,
         rangeid = first.rangeid,
         olapcount)
print(paste(pull(selection1,rangeid) %>% unique() %>% length(), " ranges with no other overlapping ranges", sep = ""))
### for overlapping repeats, discard the one with shorter overall width
### logic is that smaller alignments overall within a larger alignment are likely to be spurious hits
### also discard logic says that there is a bigger alignment compared to the discarded one. The biggest will never be discarded
discard1 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width > second.X.width &
                                                           first.X.start <= second.X.start &
                                                           first.X.end >= second.X.end
) %>%
  pull(second.rangeid)
discard2 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width < second.X.width &
                                                           first.X.start >= second.X.start &
                                                           first.X.end <= second.X.end) %>%
  pull(first.rangeid)
discard <- c(discard1, discard2)
discard <- unique(discard)
print(paste(length(discard), " ranges discarded in favour of the larger ones overlapping them", sep = ""))
selection2 <- repolap %>% filter(olapcount > 1 & ! first.rangeid %in% discard) %>% select(chr = first.X.seqnames,
                                                                                          start = first.X.start,
                                                                                          end = first.X.end,
                                                                                          width = first.X.width,
                                                                                          repid = first.repid,
                                                                                          replen = first.replen,
                                                                                          alncounts=first.alncounts,
                                                                                          meande = first.meande,
                                                                                          rcopy = first.rcopy,
                                                                                          rangeid = first.rangeid,
                                                                                          olapcount)
# selection2 needs to be of unique ranges. redundancy removed by selecting alignments with maximum overlap count, minimum row number and smallest meande
selection2 <- selection2 %>% group_by(rangeid) %>% filter(olapcount == max(olapcount) & meande == min(meande)) %>% filter(row_number()==1)
print(paste(nrow(selection2), " ranges selected in favour of the shorter/less optimal ones overlapping them", sep = ""))
selection <- bind_rows(selection1, selection2) %>% arrange(chr, start)
