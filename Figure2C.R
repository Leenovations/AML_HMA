library(readxl)
library(ggplot2)
library(tidyverse)
#--------------------------------------------------------------------#
Data2 <- read.table('Figure2C.txt',
                    sep='\t',
                    col.names=c('Region', 'Type', 'Counts'))
#--------------------------------------------------------------------#
Data2$Region <- c('1 to 5kb','1 to 5kb','3UTR', '3UTR', '5UTR', '5UTR', 'Exons', 'Exons', 'Intergenic' ,'Intergenic', 'Intronexonboundaries', 'Intronexonboundaries', 'Introns', 'Introns', 'Promoters', 'Promoters')
#--------------------------------------------------------------------#
plot <- ggplot(Data2, aes(x=Region, y=Counts, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(), color='black', width = .8) +
  scale_fill_manual(values=c("forestgreen", "firebrick1"), labels = c("CR Hypermethylation","CR Hypomethylation")) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('') +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5),
        axis.title.y=element_text(size=13, margin = margin(r=10), color='black'),
        axis.text.x = element_text(face ="bold", size = 10, angle = 45, hjust=1, color='black'),
        axis.text.y = element_text(color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'top')
#--------------------------------------------------------------------#