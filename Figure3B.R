library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(parallel)
library(data.table)
library(rtracklayer)
library(ggprism)
#------------------------------------------------------------------------------#
TF_files <- list.files('Figure3B/', pattern = "\\.bed$", full.names = TRUE)
TFSymbol <- sub(".*/(.*)\\.DMR\\.TFBS\\.bed$", "\\1", TF_files)
#------------------------------------------------------------------------------#
for (tf in TFSymbol){
  Data <- read.table(paste0('Figure3B/AML.Normal.', tf, '.TF.mean.meth.tsv'),
                     sep='\t',
                     header=TRUE)
  #------------------------------------------------------------------------------#
  Plot <- ggplot(Data, aes(Bin, Methylation)) +
    geom_smooth(aes(color=Group), span=0.25, se = FALSE, size=.65) +  #Hyper는 span 0.1로 조정할 것
    scale_color_manual(values=c("CR"="forestgreen", "NR"="firebrick1", "Normal"="orange"),
                       breaks = c("CR", "NR", "Normal")) +
    scale_x_continuous(labels=c("-5000","0", "+5000"), breaks=c(1,51,101)) +
    ggtitle(tf) +
    xlab('Distance to TF border (bp)') +
    ylab('DNA methylation') +
    # ylim(0.5, 1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=15, hjust=0.5, face='italic'),
          axis.title.x = element_text(size = 13, color="black", margin=margin(t=10)),
          axis.title.y = element_text(size = 13, color="black", margin=margin(r=10)),
          axis.ticks=element_line(color="black"),
          axis.text.x=element_text(size=11, color="black"),
          axis.text.y=element_text(size=11, color="black"),
          legend.title = element_blank(),
          legend.position = 'top',
          legend.direction = 'horizontal',
          legend.text = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'))
  #------------------------------------------------------------------------------#
}