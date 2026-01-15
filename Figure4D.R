library(readxl)
library(ggplot2)
#------------------------------------------------------------------------------#
Metascape_CR <- read_excel('Figure4D.xlsx',
                           sheet = 'Enrichment')
Metascape_CR <- as.data.frame(Metascape_CR)
Metascape_CR <- Metascape_CR[grep('Summary', Metascape_CR$GroupID), ]
Metascape_CR <- Metascape_CR[,c(3, 4, 5)]
Metascape_CR$LogP <- as.numeric(-Metascape_CR$LogP)
Metascape_CR$Description <- toupper(Metascape_CR$Description)
Metascape_CR$Annotation <- paste0(Metascape_CR$Description, ' \n(', Metascape_CR$Term, ')')
Metascape_CR <- Metascape_CR[c(1,2,3,4,7,9,10,13,15,16,17),]
Metascape_CR$Color <- c('red','black','red','red','black','black','black','black','black','black','black')
#------------------------------------------------------------------------------#
Plot <- ggplot(Metascape_CR, aes(x=LogP, y=reorder(Annotation, LogP), fill='red')) +
  geom_bar(stat="identity", width = 0.5, color='black') +
  scale_fill_manual(values='forestgreen') +
  labs(x=(expression(-log[10] ~ italic(P)))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5),
        axis.title.x=element_text(size=13, color='black', margin=margin(t=10)),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size = 11, color='black'),
        axis.text.y = element_text(size = 8, color=Metascape_CR$Color, face = 'bold'),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'none')
#------------------------------------------------------------------------------#