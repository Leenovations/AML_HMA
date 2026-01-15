library(ggplot2)
library(smplot2)
#------------------------------------------------------------------------------#
TF_list <- read.table('Figure3D.txt', sep='\t', header=T)
for (TF in TF_list$GeneSymbol){
  Pearson_data <- read.table(paste0('Figure3C/', TF, '.pearson.tsv'),
                             sep='\t',
                             header=TRUE)  # Pearson_data$Expression <- log(Pearson_data$Expression, 10)
  Pearson_data$Expression <- log(Pearson_data$Expression * 1000000, 10)
  #------------------------------------------------------------------------------#
  Plot <- ggplot(Pearson_data, aes(x=Methylation, y=Expression)) +
    geom_point(pch=1, size =2, stroke = 2, aes(color=Group)) +
    scale_color_manual(values=c('CR'='forestgreen', 'NR'='firebrick1')) +
    theme_bw() +
    ggtitle(TF) +
    xlab('DNA methylation (TF binding site)') + 
    ylab('Gene expression (TPM)') +
    sm_statCorr(color = 'black', 
                corr_method = 'pearson',
                linetype = 'dashed',
                separate_by='\n',
                text_size = 4.5,
                label_x=max(Pearson_data$Methylation)- 0.09,
                label_y=max(Pearson_data$Expression) - 0.1) +
    # annotate("text", x=max(Pearson_data$Methylation)-0.05, y=max(Pearson_data$Expression), label= TF, size=4.5) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face='italic'),
          axis.title.x = element_text(margin = margin(t = 10,b=5), size=13),
          axis.title.y = element_text(margin = margin(r = 10,l=5), size=13), legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.7, 'cm'),
          legend.direction = 'horizontal',
          legend.position = "top")
}