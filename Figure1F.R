#------------------------------------------------------------------------------#
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(ggrepel) 
#------------------------------------------------------------------------------#
Data <- read.table('Figure1F.txt',
                   sep='\t',
                   header=TRUE)
Data <- Data[,c(1:24)]
Data <- Data/100
Data <- na.omit(Data)
#------------------------------------------------------------------------------#
pca_res <- prcomp(t(Data))
#------------------------------------------------------------------------------#
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
samp_data <- data.frame(
  "Group" = rep(c('CR','NR'),c(12,12)),
  row.names = colnames(Data)
)
#------------------------------------------------------------------------------#
Plot <- ggbiplot(pca_res, ellipse = F, 
         var.axes = F, 
         groups = samp_data$Group,
         var.scale = 1) +
  geom_point(aes(colour = samp_data$Group, fill = samp_data$Group), size = 5, stroke=1.5, shape=21) +
  # geom_label_repel(aes(label = rownames(samp_data)),
  #                  box.padding   = 0.5,
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = Inf) +
  theme_bw() +
  coord_fixed(ratio = 0.8) +
  scale_fill_manual(values=c("CR"="forestgreen", "NR"="firebrick1")) +
  scale_colour_manual(values=c("CR"="black", "NR"="black")) +
  guides(colour = FALSE) +
  labs(color = "Cancer type") +
  ggtitle("PCA of Promoter") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(size=13, margin = margin(r=10), color='black'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "top")
#------------------------------------------------------------------------------#
Data <- read.table('Figure1F.txt',
                   sep='\t',
                   header=TRUE)
Data <- Data/100
Data <- na.omit(Data)
#------------------------------------------------------------------------------#
pca_res <- prcomp(t(Data))
#------------------------------------------------------------------------------#
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1", "Normal" = "orange")
samp_data <- data.frame(
  "Group" = rep(c('CR','NR','Normal'),c(12,12,3)),
  row.names = colnames(Data)
)
#------------------------------------------------------------------------------#
Plot <- ggbiplot(pca_res, ellipse = F, 
                 var.axes = F, 
                 groups = samp_data$Group,
                 var.scale = 1) +
  geom_point(aes(colour = samp_data$Group, fill = samp_data$Group), size = 5, stroke=1.5, shape=21) +
  # geom_label_repel(aes(label = rownames(samp_data)),
  #                  box.padding   = 0.5,
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = Inf) +
  theme_bw() +
  scale_fill_manual(values=c("CR"="forestgreen", "NR"="firebrick1", "Normal"="orange"),
                    breaks = c("CR", "NR", "Normal")) +
  scale_colour_manual(values=c("CR"="black", "NR"="black", "Normal"="black")) +
  guides(colour = FALSE) +
  coord_fixed(ratio = 0.8) +
  labs(color = "Cancer type") +
  ggtitle("PCA of Promoter") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(size=13, margin = margin(r=10), color='black'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "top")
#------------------------------------------------------------------------------#