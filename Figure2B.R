library(readxl)
library(ggplot2)
library(tidyverse)
library(genomation)
library(annotatr)
#--------------------------------------------------------------------#
Data1 <- read.table('Figure2B.txt',
                    sep='\t',
                    col.names=c('Region','Counts'))

Data1$Perc <- round((Data1$Counts / sum(Data1$Counts)) * 100, 2)
Data1 <- Data1[order(Data1$Perc), ]

Data1 <- Data1 %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

Data1 <- Data1[order(Data1$Perc), ]
Data1$pos[8] <- 25
Data1$Region <- c('3UTR', '5UTR', 'Promoters', '1 to 5kb', 'Intronexonboundaries', 'Exons', 'Introns', 'Intergenic')
#-----------------------------------------------------------------------#
Plot <- ggplot(Data1, aes(x = "" , y = Perc, fill = fct_inorder(Region))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = rev(brewer.pal(n = length(unique(Data1$Region)), name = "Spectral"))) +
  geom_label_repel(data = Data1,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void() +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank())