library(ggplot2)
library(dplyr)
library(plyr)
require(data.table)

parseVarscan <- function(file){
  cat("Parsing Varscan file", file, '\n')
  snvs <- read.table(file, header = T, stringsAsFactors = F)
  
  snvs <- snvs %>% 
    dplyr::mutate(pos = as.numeric(as.character(position)),
                  nFreq = suppressWarnings(round(as.numeric(gsub("\\%", "", normal_var_freq)))),
                  tFreq = suppressWarnings(round(as.numeric(gsub("\\%", "", tumor_var_freq)))),
                  ndepth = as.integer(normal_reads1) + as.integer(normal_reads2),
                  tdepth = as.integer(tumor_reads1) + as.integer(tumor_reads2)
                  ) %>% 
    dplyr::mutate(n_dist = abs(50-nFreq),
                  t_dist = abs(50-tFreq)
                  ) %>% 
    dplyr::filter(!is.na(pos),
                  ndepth >= 20, tdepth >= 20,
                  somatic_status %in% c("Germline", "LOH")) %>% 
    dplyr::select(chrom, pos, ndepth, nFreq, n_dist, tdepth, tFreq, t_dist) %>% 
    droplevels()
  
  l_snvs <- melt(setDT(snvs), measure.vars = patterns("_dist$", "Freq$"), variable.name = "condition", value.name = c("dist", "freq"))
  l_snvs$condition <- factor(ifelse(l_snvs$condition == 1, 'normal', 'tumour'))

  return(as.data.frame(l_snvs))
}


alleleFractionDepth <- function(df=NULL){
  if(is.null(df)) df <- getAllFreqs()
  
}


plotFreq <- function(df=NULL, inFile='data/D050R01.snp.LOH.hc', sample=NULL, tissue="tumour"){
  if(is.null(df)) df <- parseVarscan(file=inFile) 
  if(is.null(sample)) sample <- strsplit(basename(inFile), "[.]")[[1]][1]
  cat("Plotting genomic snv frequencies\n")
  
  snvs <- df %>%
    dplyr::filter(chrom %in% c("2L", "2R", "3L", "3R", "X")) %>%
    dplyr::filter(condition == tissue) %>% 
    dplyr::mutate(pos = pos/1e6,
                  allele_freq = freq/100) %>% 
    droplevels()
 
  p <- ggplot(snvs)
  p <- p + geom_point(aes(pos, allele_freq), size=0.1, alpha=.3)
  p <- p + scale_y_continuous("Variant Allele Frequency")
  p <- p + scale_x_continuous("Genomic position")
  p <- p + facet_wrap(~chrom, ncol = 5, scales = 'free_x')
  p <- p + ggtitle(paste(sample, " - ", tissue))
  p <- p + theme_bw()
  print(p)
}

compare_freqs <- function(df=NULL, inFile='data/D050R01.snp.LOH.hc', sample=NULL){
  if(is.null(df)) df <- parseVarscan(file=inFile) 
  if(is.null(sample)) sample <- strsplit(basename(inFile), "[.]")[[1]][1]
  
  df <- df %>% 
    dplyr::filter(!chrom %in% c("X", "Y")) %>% 
    dplyr::group_by(condition) %>% 
    dplyr::summarise(av_freq = mean(freq, na.rm=TRUE)) %>% 
    dplyr::mutate(sample = sample)
  return (df)
}


getAllFreqs <- function(varScanDir = '/Volumes/perso/Analysis/Analysis/Varscan'){
  
  groups = c("HUM", "A572", "A370", "A512", "B241", "A785")
  
  allFreqs <- list()
  
  for(i in 1:length(groups)){
    group <- groups[[i]]
    cat(group, "\n")
    files <- dir(paste(varScanDir, group, 'high_conf/', sep='/'), pattern = ".snp.hc$")
    
    group_wise_freqs <- list()
    for (j in 1:length(files)){
      file <- files[[j]]
      freqs <- compare_freqs(inFile=paste(varScanDir, group, 'high_conf', file, sep='/'))
      group_wise_freqs[[j]] <- freqs
    }
    
    gw_df <- as.data.frame(do.call(rbind, group_wise_freqs))
    allFreqs[[i]] <- gw_df
  }
  return(as.data.frame(do.call(rbind, allFreqs)))
}


plotDen <- function(df=NULL, inFile='data/D050R01.snps.txt'){
  if(is.null(df)) df <- parseVarscan(file=inFile) 
  cat("Plotting snv density\n")
  
  snvs <- df %>%
    dplyr::filter(chrom %in% c("2L", "2R", "3L", "3R", "X")) %>%
    dplyr::mutate(pos = pos/1e6) %>% 
    droplevels()
  
  p <- ggplot(snvs)
  p <- p + geom_density(aes(freq, fill = condition), alpha=.5)
  p <- p + facet_wrap(~chrom, ncol = 2, scales = 'free_x')
  print(p)
}

plotDist <- function(df=NULL, inFile='data/D050R01.snp.LOH.hc', sample=NULL){
  if(is.null(df)) df <- parseVarscan(file=inFile) 
  if(is.null(sample)) sample <- strsplit(basename(inFile), "[.]")[[1]][1]
  
  cat("Plotting genomic snv distances from hetrozygosity\n")
  
  snvs <- df %>%
    dplyr::filter(chrom %in% c("2L", "2R", "3L", "3R")) %>%
    dplyr::mutate(pos = pos/1e6) %>% 
    droplevels()
  
  palette = c("#E7B800", "#00AFBB")
  p <- ggplot(snvs)
  p <- p + geom_jitter(aes(pos, dist, colour = condition), size=0.2, alpha=0.1)
  p <- p + geom_smooth(aes(pos, dist, colour = condition), span = 0.1)
  p <- p + scale_y_continuous("Distance from Hetrozygosity")
  p <- p + scale_x_continuous("Genomic position")
  p <- p + facet_wrap(~chrom, ncol = 5, scales = 'free_x')
  p <- p + ggtitle(sample)
  p <- p + scale_colour_manual(values = palette)
  p <- p + theme_bw()
  print(p)
}


plotAllFreqs <- function(group='HUM', varScanDir = '/Volumes/perso/Analysis/Analysis/Varscan'){

  files <- dir(paste(varScanDir, group, 'high_conf/', sep='/'), pattern = ".snp.hc$")
  
  for (f in files){
    plotFreq(inFile=paste(varScanDir, group, 'high_conf', f, sep='/'))
  }

}



