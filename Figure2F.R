library(ggpubr)
library(ggplot2)
library(stats)
library(multcomp)
library(ggpubr)
library(rstatix)
library(ggsignif)
#------------------------------------------------------------------------------#
GOGene_list <- read.table('Figure2F/GO.GeneList.tsv', 
                          sep='\t', 
                          header=T)
#------------------------------------------------------------------------------#
for (gene in sort(GOGene_list$GeneSymbol)){
  Meth <- read.table(paste0('Figure2F/', gene, '.boxplot.input.tsv'), 
                     sep='\t',
                     header=T)[4:30]
  Meth <- t(Meth)
  Meth <- as.data.frame(Meth)
  nums <- c(1:ncol(Meth))
  Prefix <- paste(rep("Case", length(nums)), nums, sep = "")
  colnames(Meth) <- Prefix
  Meth[,1:ncol(Meth)] <- Meth[, 1:ncol(Meth)] / 100
  Meth$Group <- rep(c('CR','NR','Normal'), c(12, 12, 3))
  #----------------------------------------------------------------------------#
  for (order in Prefix){
    Meth_sub <- Meth[, c(order, "Group")]
    colnames(Meth_sub) <- c('order', 'Group')
    Plot <- ggplot(Meth_sub, aes(x=Group, y=order, fill=Group)) + 
      geom_boxplot(lwd = 0.9, width = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = c("CR"="forestgreen", "NR"="firebrick1", "Normal"="orange"),
                        breaks = c('Normal', 'CR', 'NR')) +
      # scale_x_discrete(labels = c('Normal', 'CR', 'NR')) + 
      # scale_color_manual(values=c("forestgreen", "firebrick1")) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 15, hjust = 0.5),
            axis.title.y = element_text(size=13, margin = margin(r=10), color='black'),
            axis.text=element_text(size=11, color='black'),
            legend.text = element_text(size = 13),
            legend.title = element_blank(),
            legend.position = "none") +
      xlab('') +
      ylab(bquote('DNA methylation ( ' * italic(.(gene)) * ' )')) +
      geom_signif(comparisons = list(c("CR", "Normal") ,c("Normal", "NR"), c("CR", "NR")),
                  color="black",
                  vjust = -0.1,
                  step_increase = 0.08,
                  textsize = 3,
                  test = "t.test",
                  map_signif_level=TRUE,
                  margin_top = 0.1)
    #------------------------------------------------------------------------------#
  }
}