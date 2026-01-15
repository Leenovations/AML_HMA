library(ggplot2)
library(ggpubr)
#------------------------------------------------------------------------------#
Data <- read.table('Figure1E.txt',
                   sep='\t',
                   header=TRUE)
#------------------------------------------------------------------------------#
Data <- subset(Data, Data$Group == 'CR' | Data$Group == 'NR')
Data$Promoter <- Data$Promoter / 100
#------------------------------------------------------------------------------#
my_comparisons <- list(c("CR", "NR"))
#------------------------------------------------------------------------------#
Plot <- ggplot(Data, aes(x=Group, y=Promoter, fill=Group)) + 
  geom_boxplot(lwd = 0.9, width = 0.5) +
  scale_color_manual(values=c("forestgreen", "firebrick1")) +
  scale_fill_manual(values=c("forestgreen", "firebrick1")) +
  geom_jitter(width=0.1) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test") +
  ggtitle("Promoter") +
  xlab('') +
  ylab('DNA methylation') +
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=15, hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=13, margin = margin(r=10), color='black'),
        axis.text.x=element_text(size=11, color="black"),
        axis.text.y=element_text(size=11, color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position = "top")
#------------------------------------------------------------------------------#