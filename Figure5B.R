library(gplots)
library(RColorBrewer)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

ht_opt$message = FALSE
Data <- read.table('', header=T)
head(Data)
Data <- Data[complete.cases(Data),]
sum(is.na(Data))

#데이터프레임 형변환#
cname <- colnames(Data[4:27])
Data[ , cname] <- lapply (Data[ , cname], as.numeric)
Data <- Data[,c(4:27)]

#데이터 검토#
Zero_row <- apply(Data,1,mean)
Data <- Data[Zero_row!=0, ]

Zero_col <- apply(Data,2,mean)
Data <- Data[Zero_col!=0,]

Zero_sd <- apply(Data,1,sd)
Data <- Data[Zero_sd!=0, ]

#Dataframe to Maxtrix#
Data <- as.matrix(Data)
Data_mean <- t(as.matrix(apply(Data, 2, mean)))

#Annotation정보#
colnames(Data) <- 

#Heatmap Annotation
samp_data <- data.frame(
  "DNMT3A" = rep(c('WT','mut','WT','mut','WT','mut','WT'),c(3,1,9,2,4,1,4)),
  "DNMT3B" = rep(c('WT','mut','WT'),c(5,1,18)),
  "IDH2" = rep(c('WT','mut','WT','mut','WT'),c(1,1,13,1,8)),
  "TET2" = rep(c('WT','mut','WT','mut','WT','mut','WT','mut','WT'),c(4,1,3,1,1,1,7,2,4)),
  "KMT2C" = rep(c('WT','mut','WT'),c(20,1,3)),
  "KMT2D" = rep(c('WT','mut','WT'),c(14,1,9)),
  "NSD1" = rep(c('WT','mut','WT'),c(18,1,5)),
  "NSD3" = rep(c('WT','mut','WT','mut','WT'),c(5,1,10,1,7)),
  "SETBP1" = rep(c('WT','mut','WT'),c(12,1,11)),
  "SUZ12" = rep(c('WT','mut','WT'),c(21,1,2)),
  "Group" = rep(c('CR','NR'),c(12,12)),
  row.names = colnames(Data)
)

Bottom_anno = HeatmapAnnotation(Group = samp_data$Group,
                                col = list(Group = c("CR" = "forestgreen",
                                                     "NR" = "firebrick1")),
                                show_annotation_name = F,
                                annotation_label=c('HMA response'),
                                annotation_name_gp= gpar(fontsize = 12),
                                annotation_legend_param = list(Group=list(title='HMA response', direction = "horizontal")))

Top_anno = HeatmapAnnotation(DNMT3A = samp_data$DNMT3A,
                             DNMT3B = samp_data$DNMT3B,
                             IDH2 = samp_data$IDH2,
                             TET2 = samp_data$TET2,
                             KMT2C = samp_data$KMT2C,
                             KMT2D = samp_data$KMT2D,
                             NSD1 = samp_data$NSD1,
                             NSD3 = samp_data$NSD3,
                             SETBP1 = samp_data$SETBP1,
                             SETD2 = samp_data$SET2D,
                             SUZ12 = samp_data$SUZ12,
                             col = list(DNMT3A = c('mut' = "black",
                                                   'WT' = "#CCCCCC"),
                                        DNMT3B = c('mut' = "black",
                                                   'WT' = "#CCCCCC"),
                                        IDH2 = c('mut' = "black",
                                                 'WT' = "#CCCCCC"),
                                        TET2 = c('mut' = "black",
                                                 'WT' = "#CCCCCC"),
                                        KMT2C = c('mut' = "black",
                                                  'WT' = "#CCCCCC"),
                                        KMT2D = c('mut' = "black",
                                                  'WT' = "#CCCCCC"),
                                        NSD1 = c('mut' = "black",
                                                 'WT' = "#CCCCCC"),
                                        NSD3 = c('mut' = "black",
                                                 'WT' = "#CCCCCC"),
                                        SETBP1 = c('mut' = "black",
                                                   'WT' = "#CCCCCC"),
                                        SETD2 = c('mut' = "black",
                                                  'WT' = "#CCCCCC"),
                                        SUZ12 = c('mut' = "black",
                                                  'WT' = "#CCCCCC")),
                             show_annotation_name = T,
                             annotation_label=c('DNMT3A', 'DNMT3B','IDH2','TET2','KMT2C','KMT2D','NSD1','NSD3','SETD2','SUZ12'),
                             annotation_name_gp= gpar(fontsize = 10, fontface = "italic"),
                             show_legend = rep(c(TRUE,FALSE),c(1,9)),
                             annotation_legend_param = list(DNMT3A=list(title='Mutation')),
                             height = unit(5, "cm"),
                             gap = unit(0.1, "mm"),
                             border = F,
                             gp = gpar(col = "white"),
                             simple_anno_size_adjust = TRUE)
#------------------------------------------------------------------------------#
#Draw Heatmap
# Heatmap_col = colorRamp2(c(0, 1), c("blue", "yellow"), space = "sRGB")
Heatmap_col = colorRampPalette(c("blue1","yellow1"))(100)
Plot <- Heatmap(Data_mean,
                col = Heatmap_col,
                cluster_rows = F,
                cluster_columns = T,
                show_row_names = F,
                show_column_names = F,
                width = ncol(Data_mean)*unit(5, "mm"), 
                height = nrow(Data_mean)*unit(5, "mm"),
                show_heatmap_legend = F,
                top_annotation = c(Top_anno, Bottom_anno),
                heatmap_legend_param = list(title="Methylation", 
                                            at = c(0,1), 
                                            labels = c('Low', 'High'),
                                            direction = "vertical",
                                            legend_width = unit(2, "cm"),
                                            title_position = "topleft")
)

col_fun <- colorRamp2(c(70, 90), c("blue", "yellow"), space = "sRGB")
lgd2 = Legend(col_fun = Heatmap_col,
              title = "Global\nmethylation",
              direction = "horizontal",
              at = c(70, 90),
              labels = c("Low", "High"),
              # border = "black",
              title_gp = gpar(fontsize = 10, fontface = "bold"),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(2.5, "cm"),
              grid_width = unit(0.45, "cm"),
              title_position = "topleft")

draw(Plot, heatmap_legend_side = "left", 
     annotation_legend_side = "left", 
     merge_legend = TRUE)
draw(lgd2, x = unit(3, "cm"), y = unit(1.2, "cm"))
#------------------------------------------------------------------------------#