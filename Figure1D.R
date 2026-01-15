library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(parallel)
library(data.table)
library(rtracklayer)
library(ggprism)
#------------------------------------------------------------------------------#
Data <- read.table('Figure1D.txt',
                   sep='\t',
                   header=TRUE)
#------------------------------------------------------------------------------#
Plot <- ggplot(Data, aes(Bin, Methylation)) +
  geom_line(aes(color=Group), lwd=0.9) +
  scale_color_manual(values=c("CR"="forestgreen", "NR"="firebrick1", "Normal"="orange"),
                     breaks=c('CR', 'NR', 'Normal')) +
  scale_x_continuous(labels=c("TSS"),breaks=c(41)) +
  geom_vline(xintercept=41,linetype="longdash",color="black", linewidth=0.3) +
  ggtitle('Promoter methylation') +
  xlab('') +
  ylab('DNA methylation') +
  # scale_y_continuous(guide = "prism_offset",
  #                    limits = c(0, 1)) +
  # scale_x_continuous(guide = "prism_offset") + 
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5),
        axis.title.y=element_text(size=13, margin = margin(r=10), color='black'),
        axis.text=element_text(color="black", size=11), 
        axis.ticks=element_line(color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position = c(.9, .2))
  # annotate("rect", xmin = 22, xmax = 62, ymin = -Inf, ymax = 0.04,
  #          alpha = .7, fill = "firebrick2") +
  # annotate("text", x=41, y=0, label = "Promoter",
  #          color="black",size=2) +
  # annotate("rect", xmin = 1, xmax = 22, ymin = -Inf, ymax = 0.04,
  #          alpha = .7, fill = "darkgreen") +
  # annotate("text", x=11.5, y=0, label = "Upstream",
  #          color="black",size=2) +
  # annotate("rect", xmin = 62, xmax = 200, ymin = -Inf, ymax = 0.04,
  #          alpha = .7, fill = "darkgreen") +
  # annotate("text", x=131, y=0, label = "Gene body",
  #          color="black",size=2)
#------------------------------------------------------------------------------#