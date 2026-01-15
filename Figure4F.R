library(readxl)
library(ggplot2)
#------------------------------------------------------------------------------#
Metascape_NR <- read_excel('Figure4F.xlsx',
                           sheet = 'Enrichment')
Metascape_NR <- as.data.frame(Metascape_NR)
Metascape_NR <- Metascape_NR[grep('Summary', Metascape_NR$GroupID), ]
Metascape_NR <- Metascape_NR[,c(3, 4, 5)]
Metascape_NR$LogP <- as.numeric(-Metascape_NR$LogP)
Metascape_NR$Description <- toupper(Metascape_NR$Description)
Metascape_NR$Annotation <- paste0(Metascape_NR$Description, ' \n(', Metascape_NR$Term, ')')
Metascape_NR <- Metascape_NR[1:11,]
Metascape_NR$Color <- rep(c('black', 'red'), c(8, 3))
#------------------------------------------------------------------------------#
Plot <- ggplot(Metascape_NR, aes(x=LogP, y=reorder(Annotation, LogP), fill='red')) +
  geom_bar(stat="identity", width = 0.5, color='black') +
  scale_fill_manual(values='firebrick1') +
  labs(x=(expression(-log[10] ~ italic(P)))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5),
        axis.title.x=element_text(size=13, color='black', margin=margin(t=10)),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size = 11, color='black'),
        axis.text.y = element_text(size = 8, colour = Metascape_NR$Color, face = 'bold'),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'none')
#------------------------------------------------------------------------------#
