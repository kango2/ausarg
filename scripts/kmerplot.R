library(tidyverse)
library(scales)

setwd("/g/data/xl04/ka6418/bassiana/publication-v2/eval/kmers")
files <- list.files("./","*.histo")

khisto <- list()
for (i in seq_along(files)) {
  f <- files[i]
  x <- str_split(str_split(f, "\\.", simplify = T)[1], "_", simplify = T)
  a <- read_delim(files[i], delim = " ", col_names = F, col_types = cols(
    X1 = col_double(),
    X2 = col_double()
  ))
  colnames(a) <- c("kcounts", "kfrequency")
  a$tech <- x[2]
  a$klen <- x[3]
  khisto[[i]] <- a
}
khisto <- bind_rows(khisto)

# Calculate peaks
peaks <- khisto %>% 
  mutate(m = kcounts * kfrequency) %>% 
  filter(kcounts > 5 & kcounts < 100) %>% 
  group_by(tech, klen) %>% 
  filter(m == max(m)) %>% 
  filter(kcounts > 10)

# Plotting
khisto %>% 
  mutate(kcounts = case_when(kcounts > 200 ~ 200, TRUE ~ kcounts)) %>% 
  group_by(tech, klen, kcounts) %>% 
  summarise(kfrequency = sum(kfrequency), .groups = 'drop') %>%
  filter((kcounts > 1 & kcounts <= 100 & kfrequency < 5e7 & tech == "Illumina") | 
           (kcounts > 2 & kcounts <= 100 & tech == "PB") | 
           (kcounts > 3 & kcounts <= 100 & kfrequency < 5e7 & tech == "ONT")) %>%
  ggplot(aes(x = kcounts, y = kfrequency)) + 
  geom_area(stat = "identity", fill = "lightblue") + 
  geom_vline(data = peaks, aes(xintercept = kcounts)) +
  geom_label(data = peaks, aes(x = kcounts, y = kfrequency * 1.5, label = kcounts)) + 
  theme_bw() +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  facet_grid(tech ~ klen, scales = "free_y") +
  xlab("Kmer counts") + 
  ylab("Frequency") +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1))

# Additional plots
khisto %>% 
  filter(kcounts > 2 & kcounts <= 200 & kfrequency < 1e8) %>%
  ggplot(aes(x = kcounts, y = kfrequency)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  facet_wrap(tech ~ klen, scales = "free_y")

khisto %>% 
  mutate(kcounts = case_when(kcounts > 200 ~ 200, TRUE ~ kcounts)) %>% 
  group_by(tech, klen, kcounts) %>% 
  summarise(kfrequency = sum(kfrequency), .groups = 'drop') %>%
  ggplot(aes(x = kcounts, y = kfrequency)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  scale_y_log10() + 
  facet_wrap(tech ~ klen, scales = "free_y")
