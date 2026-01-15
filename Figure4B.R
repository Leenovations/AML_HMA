library(readxl)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
#------------------------------------------------------------------------------#
Data <- read.table('Figure4B.xls', header=T)
#------------------------------------------------------------------------------#
Data$"DEG" <- ifelse(Data$p.value<0.01 & Data$m.value > 2, "NR",
                      ifelse(Data$p.value<0.01 & Data$m.value < -2, "CR","None"))
Data_sub <- subset(Data, Data$DEG != "None")
Gene <- Data_sub$gene_id
#------------------------------------------------------------------------------#
Expression <- read.table('Figure4B.mtx', header=T, row.names = 1)
Exp <- Expression[Gene,]
Sample_name <- colnames(Expression)
#------------------------------------------------------------------------------#
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
samp_data <- data.frame(
  "GROUP" = rep(c('CR','NR'),c(8,5)), 
  row.names = Sample_name
)
#------------------------------------------------------------------------------#
lgd1 = Legend(labels = c('CR', 'NR'), 
              title = 'HMA response',
              labels_gp = gpar(fontsize = 10),
              legend_gp = gpar(fill = c('forestgreen', 'firebrick1')),
              title_position = "topleft")

lgd = packLegend(lgd1, direction = "vertical")
#------------------------------------------------------------------------------#
Exp <- as.matrix(Exp)
Exp <- log(Exp+1,2)
Exp <- na.omit(Exp)
Exp <- as.matrix(Exp)
#------------------------------------------------------------------------------#
Plot <- pheatmap(mat = Exp,
         name = 'Expression',
         cellwidth=25,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_col = samp_data,
         annotation_colors = list("GROUP" = samp_colors),
         clustering_method = 'complete',
         legend_breaks = c(-2.5,2.5),
         legend_labels = c("Low", "High"),
         annotation_legend = F,
         legend = T,
         show_rownames = F, 
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         scale = "row",
         annotation_names_col=F,
         border_color = NA)
#------------------------------------------------------------------------------#
pdf("/labmed/02.AML/00.Code/Figure4B.pdf")
draw(Plot)
draw(lgd, x = unit(15.85, "cm"), y = unit(10, "cm"))
dev.off()
#------------------------------------------------------------------------------#