library(readxl)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
#------------------------------------------------------------------------------#
Data <- read.table('/labmed/02.AML/01.Input/Figure4B/deseq2_DEG.Over.50.xls', header=T)
#------------------------------------------------------------------------------#
Data$"DEG" <- ifelse(Data$p.value<0.01 & Data$m.value > 1.5, "NR",
                     ifelse(Data$p.value<0.01 & Data$m.value < -1.5, "CR","None"))
Data_sub <- subset(Data, Data$DEG != "None")
Gene <- Data_sub$gene_id
#------------------------------------------------------------------------------#
Expression <- read.table('/labmed/02.AML/01.Input/Figure4B/count.deseq2.Over.50.mtx', header=T, row.names = 1)
Exp <- Expression[Gene,]
#------------------------------------------------------------------------------#
pca_res <- prcomp(t(Exp))
#------------------------------------------------------------------------------#
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
samp_data <- data.frame(
  "Group" = rep(c('CR','NR'),c(8,5)),
  row.names = colnames(Exp)
)
#------------------------------------------------------------------------------#
Plot <- ggbiplot(pca_res, ellipse = F, 
                 var.axes = F, 
                 groups = samp_data$Group,
                 var.scale = 1) +
  geom_point(aes(colour = samp_data$Group, fill = samp_data$Group), size = 5, stroke=1.5, shape=21) +
  theme_bw() +
  scale_fill_manual(values=c("CR"="forestgreen", "NR"="firebrick1")) +
  scale_colour_manual(values=c("CR"="black", "NR"="black")) +
  guides(colour = FALSE) +
  coord_fixed(ratio = 1.1) +
  ggtitle("PCA of DEG") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.x = element_text(size = 13, margin = margin(t=10)),
        axis.title.y = element_text(size = 13, margin = margin(r=10)),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "top")
#------------------------------------------------------------------------------#
ggsave(Plot, file="/labmed/02.AML/00.Code/Figure4C.pdf", height=5, width=5)
#------------------------------------------------------------------------------#