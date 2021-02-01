# plot germline vars

library(tidyverse)

snps <- read.delim('combined_snps.txt', header = T)

library(stringr)

snps$status <- as.factor(ifelse(stringr::str_detect(snps$samples, ','), 'recurrent', 'private'))

total <- nrow(snps)

types_count <- snps %>% 
  dplyr::filter(!samples %in% c('D050R07-1')) %>%
  dplyr::group_by(status) %>% 
  dplyr::summarise(count = n(), perc = count/total * 100)


types_count %>% ggplot(.) +
  geom_bar(aes(status, perc), stat = 'identity') +
  scale_y_continuous("Percent contribution") +
  theme_bw()

snps %>% 
  dplyr::filter(status == 'private') %>% 
  dplyr::group_by(samples) %>% 
  dplyr::summarise(count = n(), perc = count/total * 100) %>% 
  ggplot(.) +
  geom_bar(aes(fct_reorder(samples, -count), count), stat = 'identity') +
  scale_y_continuous("Private SNP count") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y =element_text(size=15)
  )
