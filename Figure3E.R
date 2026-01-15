library(readxl)
library(ggplot2)
#------------------------------------------------------------------------------#
Metascape_TF <- read_excel('Figure3E.xlsx',
                           sheet = 'Sheet1')
Metascape_TF <- as.data.frame(Metascape_TF)
# Metascape_TF <- Metascape_TF[grep('Summary', Metascape_TF$GroupID), ]
Metascape_TF <- Metascape_TF[c(1, 8, 11, 22, 55, 59, 64, 71, 76, 206, 279),]
Metascape_TF <- Metascape_TF[,c(3, 4, 5)]
Metascape_TF$LogP <- as.numeric(-Metascape_TF$LogP)
Metascape_TF$Description <- toupper(Metascape_TF$Description)
Metascape_TF$Annotation <- paste0(Metascape_TF$Description, ' \n(', Metascape_TF$Term, ')')
# Metascape_TF <- Metascape_TF[1:11,]
#------------------------------------------------------------------------------#
Plot <- ggplot(Metascape_TF, aes(x=LogP, y=reorder(Annotation, LogP), fill='ivory2')) +
  geom_bar(stat="identity", width = 0.6, color='black') +
  scale_fill_manual(values='ivory2') +
  labs(x=(expression(-log[10] ~ italic(P)))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5),
        axis.title.x=element_text(size=13, color='black', margin=margin(t=10)),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size = 11, color='black'),
        axis.text.y = element_text(size = 8, color='black', face = 'bold'),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'none')
#------------------------------------------------------------------------------#