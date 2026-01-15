library(data.table)
library(pheatmap)
library(ComplexHeatmap)
#------------------------------------------------------------------------------#
ssGSEA_result<-fread("Figure4H.gct")
ssGSEA_result<-as.data.frame(ssGSEA_result)
rownames(ssGSEA_result)<-ssGSEA_result$Name
rownames(ssGSEA_result) <- toupper(gsub("_", " ", rownames(ssGSEA_result)))
ssGSEA_result<-ssGSEA_result[,-c(1:2)]
ssGSEA_result<-na.omit(ssGSEA_result)
col_annotation <- data.frame(type=rep(c('CR','NR'),c(8,5)))
#------------------------------------------------------------------------------#
rownames(col_annotation) <- colnames(ssGSEA_result)
samp_colors <- c("CR" = "forestgreen", "NR" = "firebrick1")
Sample_name <- colnames(ssGSEA_result)
samp_data <- data.frame(
  "GROUP" = rep(c('CR','NR'),c(8,5)), 
  row.names = Sample_name
)
#------------------------------------------------------------------------------#
P_value <- c()
for (i in 1:nrow(ssGSEA_result)){
  Data <- ssGSEA_result[i,]
  Data <- t(Data)
  Data <- as.data.frame(Data)
  Data$Group <- rep(c('CR','NR'),c(8,5))
  Value <- t.test(Data[,1]~ Data[,2])
  P_value <- c(P_value, Value$p.value)
}
#------------------------------------------------------------------------------#
P_ssgsea <- ssGSEA_result
P_ssgsea$p_value <- P_value

Sig_ssgsea <- subset(P_ssgsea, P_ssgsea$p_value < 0.05)
Sig_ssgsea <- Sig_ssgsea[,-14]
#------------------------------------------------------------------------------#
lgd1 = Legend(labels = c('CR', 'NR'), 
              title = 'HMA response',
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_gp = gpar(fill = c('forestgreen', 'firebrick1')),
              title_position = "topleft")

col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "red"))

lgd2 = Legend(col_fun = col_fun,
              title = "Enrichment",
              direction = "vertical",
              at = c(-2, 2),
              labels = c("Low", "High"),
              # border = "black",
              title_gp = gpar(fontsize = 11, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(3, "cm"),
              grid_width = unit(0.45, "cm"),
              title_position = "topleft")

lgd = packLegend(lgd1, lgd2, direction = "vertical")
#------------------------------------------------------------------------------#
Plot <- pheatmap(Sig_ssgsea, 
         annotation_col = samp_data,
         cellwidth=20,
         cellheight=20,
         color=colorRampPalette(c("navy", "white", "red"))(100),
         annotation_colors = list("GROUP" = samp_colors),
         scale= 'row',
         annotation_legend = F,
         legend = F,
         # legend_breaks = c(-2,2),
         # legend_labels = c("Low", "High"),
         show_colnames = F,
         treeheight_col= 0,
         show_rownames = T,
         cluster_cols = T,
         cluster_rows = T,
         border_color = NA,
         annotation_names_col=F)
