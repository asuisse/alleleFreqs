# plot germline vars

library(tidyverse)
library(cowplot)

infile <- 'snps_dirty_combined.txt'
infile <- 'snps_filtered_combined.txt'

snps <- read.delim(infile, header = T)

total <- nrow(snps)

dark_blue <- '#21608A'
mint_green <- '#26AD7E'
red <- '#D9564A'
purple <- '#B57CE6'

colours <- c(dark_blue, mint_green, red, purple)

types_count <- snps %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::summarise(count = n(), perc = count/total * 100)

# Percentage contribution of genotypes
# dirty_genotype_contributions

types_count %>% ggplot(.) +
  geom_bar(aes(fct_reorder(genotype, -perc), perc), fill=mint_green, alpha=0.6, stat = 'identity') +
  scale_y_continuous("Percent contribution") +
  theme_bw()


dfsObj <- paste0( "~/Desktop/script_test/Riddiford_et_al_2020/muts.RData")
load(dfsObj)


# Count SNVs per sample
snv_count <- male_snvs %>% 
  filter(af > 0) %>% 
  droplevels() %>% 
  dplyr::mutate(sample = sample_paper) %>% 
  group_by(sample) %>% 
  dplyr::summarise(somatic_tumour = n())


# Count SNPs per sample
snp_count <- snps %>% 
  dplyr::filter(genotype == 'germline_private') %>% 
  droplevels() %>% 
  dplyr::mutate(sample = samples) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(germline_private = n())

normal_count <- snps %>% 
  dplyr::filter(genotype == 'somatic_normal') %>% 
  droplevels() %>% 
  dplyr::mutate(sample = samples) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(normal_only = n())

tumour_count <- snps %>% 
  dplyr::filter(genotype == 'somatic_tumour') %>% 
  droplevels() %>% 
  dplyr::mutate(sample = samples) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(tumour_only = n())


gr_depth <- snps %>% 
  dplyr::filter(genotype == 'germline_private') %>%
  droplevels() %>% 
  dplyr::mutate(depth = as.numeric(as.character(depth)),
                sample = samples) %>% 
  dplyr::select(sample, depth)

  # dplyr::group_by(sample) %>%
  # summarise(count = n(), av_depth = mean(depth))


order <- levels(fct_reorder(snp_count$sample, -snp_count$germline_private))
gr_depth$sample <- factor(gr_depth$sample, levels = order, labels=order)

depth_p <- gr_depth %>% 
  dplyr::filter(depth < 100) %>% 
  ggplot(.) +
  geom_boxplot(aes(sample, depth), fill='grey',  alpha=0.6) +
  scale_y_continuous("Depth") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    legend.position = 'none'
  )

# combined <- plyr::join_all(list(snp_count, normal_count, tumour_count, snv_count), by='sample')
combined <- plyr::join_all(list(snp_count, snv_count), by='sample')

df_long <- tidyr::gather(combined, type, val, germline_private:somatic_tumour, factor_key=TRUE) %>% 
  dplyr::group_by(sample) %>% 
  tidyr::complete(type, fill = list(n=0, perc=0))

df_long$val[is.na(df_long$val)] <- 0

df_long$sample <- factor(df_long$sample, levels = order, labels=order)

# Fig 1 - private SNPs vs somatic SNVs


threshold = 200

count_p <- df_long %>% 
  ggplot(.) +
  geom_bar(aes(sample, val, fill=type),  alpha=0.6, stat = 'identity') +
  geom_hline(aes(yintercept=threshold), color = red, linetype = "solid") +
  geom_text(aes(20, threshold,label = threshold, vjust = -1)) +
  scale_y_continuous("Mutation count") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    legend.position = 'none'
  ) +
  scale_fill_manual(values = colours) +
  facet_wrap(~type, ncol = 1, scales = 'free_y')

count_p


plot_grid(count_p, depth_p, ncol = 1, rel_heights = c(2,1))

# # plot SNVs per sample
# snv_p <- combined %>% 
#   ggplot(.) +
#   geom_bar(aes(fct_reorder(sample, -snp), snv), fill=mint_green, alpha=0.6, stat = 'identity') +
#   scale_y_continuous("Somatic tumour SNV count") +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
#     axis.text.y = element_text(size=12),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )
# 
# # plot gerSNPs per sample
# snp_p <- combined %>% 
#   ggplot(.) +
#   geom_bar(aes(fct_reorder(sample, -snp), snp), fill=dark_blue, alpha=0.6, stat = 'identity') +
#   scale_y_continuous("Private germline SNP count") +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
#     axis.text.y = element_text(size=12),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )
# 
# normal_p <- combined %>% 
#   ggplot(.) +
#   geom_bar(aes(fct_reorder(sample, -snp), normal), fill=red, alpha=0.6, stat = 'identity') +
#   scale_y_continuous("Somatic normal SNV count") +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
#     axis.text.y = element_text(size=12),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )


# library(cowplot)
# 
# labels = c('Tumour-only', 'Germline-private', 'Normal-only')
# plot_grid(snv_p, snp_p, normal_p, labels=labels,  ncol = 1)

# Germline events shared by multiple samples
snps %>% 
  dplyr::filter(genotype %in% c('germline_private', 'germline_recurrent')) %>% 
  dplyr::group_by(sharedby) %>% 
  dplyr::tally() %>% 
  dplyr::mutate(sharedby = as.factor(sharedby)) %>% 
  ggplot(.) +
  geom_bar(aes(sharedby, n), fill=dark_blue, alpha=0.6, stat = 'identity') +
  scale_y_continuous("Number of SNPS") +
  scale_x_discrete("Shared between samples") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title =element_text(size=15)
  )


# Average depth VS count for germline_private events
snps %>% 
  dplyr::filter(genotype == 'germline_private') %>%
  droplevels() %>% 
  dplyr::mutate(depth = as.numeric(as.character(depth))) %>% 
  dplyr::group_by(samples) %>%
  summarise(count = n(), av_depth = mean(depth)) %>%
  ggplot(.,aes(av_depth, count)) +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  scale_x_continuous("Average depth") +
  scale_y_continuous("Germline private SNP count") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=15)
  )


## VAF of SNPs/SNVs
hum_1 <- read.delim('filtered/B241R59_snps.txt', header=T)

gl_p <- hum_1 %>% 
  dplyr::filter(genotype %in% c('germline')) %>% 
  ggplot(.) +
  geom_density(aes(tumour_vaf, fill=purple), alpha=.6) +
  scale_x_continuous("Allele Frequency") +
  scale_y_continuous("Density") +
  scale_fill_manual(values = purple) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=15),
    legend.position = 'none'
  ) +
  ggtitle("Germline")


snv_p <- male_snvs %>% 
  ggplot(.) +
  geom_density(aes(af, fill=mint_green), alpha=.6) +
  scale_x_continuous("Allele Frequency") +
  scale_y_continuous("Density") +
  scale_fill_manual(values = mint_green) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=15),
    legend.position = 'none'
  ) +
  ggtitle("Somatic")

plot_grid(gl_p, snv_p, ncol = 1)


# Write somatic calls for indels + snvs
male_snvs %>% 
  dplyr::mutate(genotype = 'somatic_tumour') %>% 
  dplyr::select(sample_paper,  chrom, pos, ref, alt, af, caller, gene, feature, id, fpkm) %>% 
  dplyr::mutate(sample = sample_paper, allele_frequency = af, calledby = caller) %>% 
  dplyr::select(sample,  chrom, pos, ref, alt, allele_frequency, calledby, gene, feature, id, fpkm) %>% 
  write.table(., file = 'somatic_snvs.txt', quote=FALSE, sep='\t', row.names = FALSE)


male_indels %>% 
  dplyr::mutate(genotype = 'somatic_tumour') %>% 
  dplyr::select(sample_paper, chrom, pos, ref, alt, af, type, caller, gene, feature, id, fpkm) %>% 
  dplyr::mutate(sample = sample_paper, allele_frequency = af, calledby = caller) %>% 
  dplyr::select(sample, chrom, pos, ref, alt, type, allele_frequency, calledby, gene, feature, id, fpkm) %>% 
  write.table(., file = 'somatic_indels.txt', quote=FALSE, sep='\t', row.names = FALSE)



