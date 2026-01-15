library(ComplexHeatmap)
library(circlize)
library(gplots)
library(ggplot2)
#------------------------------------------------------------------------------#
Total_Data <- read.table('Figure1A.txt', header=T)
Total_Data[1:nrow(Total_Data), 4:ncol(Total_Data)] <- Total_Data[1:nrow(Total_Data), 4:ncol(Total_Data)] / 100
#------------------------------------------------------------------------------#
# TEST <- subset(Total_Data, Total_Data$Chromosome == 'chr22')
#------------------------------------------------------------------------------#
CR <- Total_Data[,c(1, 2, 3, 4:15)]
NR <- Total_Data[,c(1, 2, 3, 16:27)]
#------------------------------------------------------------------------------#
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22)) 
col_fun = colorRamp2(c(0, 1), c("mediumblue","yellow"), space='sRGB', transparency=0)
circos.genomicHeatmap(CR, col = col_fun, side = "inside",
                      line_col = NA, heatmap_height=0.2, connection_height = mm_h(.5))
circos.genomicHeatmap(NR, col = col_fun, side = "inside",
                      line_col = NA, heatmap_height=0.2, connection_height = mm_h(.5)) +
  lines(c(0.3, 0.45), c(0.67, 0.9), col = "black", lwd = 2) +
  lines(c(0.45, 0.9), c(0.9, 0.9), col = "black", lwd = 2) +
  points(c(0.3), c(0.67), col = "black", pch = 19, cex = 1) +
  text(x = 0.9, y = 0.9, labels = " CR", col = "black", cex = 1, pos = 4) + 
  #------------------------------------------------------------------------------#
  lines(c(0.3, 0.6), c(0.37, 0.8), col = "black", lwd = 2) +
  lines(c(0.6, 0.9), c(0.8, 0.8), col = "black", lwd = 2) +
  points(c(0.3), c(0.37), col = "black", pch = 19, cex = 1) +
  text(x = 0.9, y = 0.8, labels = " NR", col = "black", cex = 1, pos = 4)
#------------------------------------------------------------------------------#
lgd = Legend(at = c(0,1), 
             col_fun = col_fun, 
             labels = c("Low", "High"),
             legend_width = unit(2, "cm"),
             grid_height = unit(0.5, "cm"),
             title_position = "topcenter", 
             title = "DNA methylation", 
             direction = "horizontal",
             title_gap=unit(3,'mm'))
draw(lgd)
#------------------------------------------------------------------------------#
circos.clear()
#------------------------------------------------------------------------------#