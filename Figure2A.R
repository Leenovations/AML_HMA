library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
library(gplots)
library(grid)
library(gridExtra)
#------------------------------------------------------------------------------#
ht_opt$message = FALSE
#------------------------------------------------------------------------------#
Data <- read.table('Figure2A.txt',
                   sep='\t',
                   header=TRUE)
Data <- Data[1:2000,c(2:25)]
Data <- na.omit(Data) / 100
Data <- as.matrix(Data)
#------------------------------------------------------------------------------#
colnames(Data) <- 
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
samp_data <- data.frame(
  "GROUP" = rep(c('CR','NR'),c(12,12)), 
  row.names = colnames(Data)
)
#------------------------------------------------------------------------------#
Bottom_anno = HeatmapAnnotation(Group = samp_data$GROUP,
                                col = list(Group = c("CR" = "forestgreen",
                                                     "NR" = "firebrick1")),
                                show_annotation_name = F,
                                annotation_label=c('HMA response'),
                                annotation_name_gp= gpar(fontsize = 12),
                                annotation_legend_param = list(Group=list(title='HMA response', direction = "horizontal")))
#------------------------------------------------------------------------------#
Plot <- pheatmap(mat = Data,
         name = 'Methylation',
         cellwidth=10,
         color = colorRampPalette(c("blue1","blue1","blue1", "yellow1","yellow1","yellow1"))(7),
         annotation_col = samp_data,
         annotation_colors = list("GROUP" = samp_colors),
         annotation_names_col = F,
         clustering_method = 'complete',
         annotation_legend = F,
         legend = F,
         # legend_breaks = c(-3,3),
         # legend_labels = c("Low", "High"),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = F,
         scale = "row",
         border_color=NA)
#------------------------------------------------------------------------------#
lgd1 = Legend(labels = c('CR', 'NR'), 
              title = 'HMA response',
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_gp = gpar(fill = c('forestgreen', 'firebrick1')),
              title_position = "topleft")

col_fun <- colorRamp2(c(-1, 0, 1), c("blue1", "gray", "yellow1"))

lgd2 = Legend(col_fun = col_fun,
              title = "Methylation",
              direction = "vertical",
              at = c(-2, 2),
              labels = c("Low", "High"),
              # border = "black",
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(2.5, "cm"),
              grid_width = unit(0.45, "cm"),
              title_position = "topleft")

lgd = packLegend(lgd1, lgd2, direction = "vertical")
#------------------------------------------------------------------------------#
Plot
draw(lgd, x = unit(15.7, "cm"), y = unit(7, "cm"))
dev.off()
#------------------------------------------------------------------------------#
DMR <- read.table('Figure2A.DMR.txt', sep='\t', header=T)
DMR <- DMR[order(DMR$meth.diff),]
DMR$Position <- paste0('chr',DMR$chr, '_', DMR$start, '_', DMR$end)
DMR_CR <- subset(DMR, DMR$meth.diff < 0)
DMR_NR <- subset(DMR, DMR$meth.diff > 0)
dim(DMR_CR)
Rank <- c(DMR_CR$Position[1:250], DMR_NR$Position[1:250])

Data <- read.table('Figure2A.perc.txt',
                   sep='\t',
                   header=TRUE)
Data$Info <- paste0('chr',Data$Chr, '_', Data$Start, '_', Data$End)
col_order <- c('Info', colnames(Data)[4:27])

Data <- Data[, col_order]
Data <- Data[Data$Info %in% Rank, ]
Data <- Data[, c(2:25)]
Data <- na.omit(Data) / 100
Data <- as.matrix(Data)
#------------------------------------------------------------------------------#
colnames(Data) <- 
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
samp_data <- data.frame(
  "GROUP" = rep(c('CR','NR'),c(12,12)), 
  row.names = colnames(Data)
)
#------------------------------------------------------------------------------#
Bottom_anno = HeatmapAnnotation(Group = samp_data$GROUP,
                                col = list(Group = c("CR" = "forestgreen",
                                                     "NR" = "firebrick1")),
                                show_annotation_name = F,
                                annotation_label=c('HMA response'),
                                annotation_name_gp= gpar(fontsize = 12),
                                annotation_legend_param = list(Group=list(title='HMA response', direction = "horizontal")))
#------------------------------------------------------------------------------#
Plot <- pheatmap(mat = Data,
                 name = 'Methylation',
                 cellwidth=10,
                 color = colorRampPalette(c("blue1","blue1","blue1", "yellow1","yellow1","yellow1"))(7),
                 annotation_col = samp_data,
                 annotation_colors = list("GROUP" = samp_colors),
                 annotation_names_col = F,
                 clustering_method = 'complete',
                 annotation_legend = F,
                 legend = F,
                 # legend_breaks = c(-3,3),
                 # legend_labels = c("Low", "High"),
                 show_rownames = F,
                 show_colnames = F,
                 cluster_rows = T,
                 cluster_cols = F,
                 scale = "row",
                 border_color=NA)
#------------------------------------------------------------------------------#
lgd1 = Legend(labels = c('CR', 'NR'), 
              title = 'HMA response',
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_gp = gpar(fill = c('forestgreen', 'firebrick1')),
              title_position = "topleft")

col_fun <- colorRamp2(c(-1, 0, 1), c("blue1", "gray", "yellow1"))

lgd2 = Legend(col_fun = col_fun,
              title = "Methylation",
              direction = "vertical",
              at = c(-2, 2),
              labels = c("Low", "High"),
              # border = "black",
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(2.5, "cm"),
              grid_width = unit(0.45, "cm"),
              title_position = "topleft")

lgd = packLegend(lgd1, lgd2, direction = "vertical")
#------------------------------------------------------------------------------#
draw(lgd, x = unit(15.7, "cm"), y = unit(11, "cm"))