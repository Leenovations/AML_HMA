library(ComplexHeatmap)
library(readxl)
library(circlize)
library(ensembldb)
library(AnnotationFilter)
#------------------------------------------------------------------------------#
Onco.Data <- read.table('Figure5A.txt',
                        sep='\t',
                        row.names = 1,
                        header = T)
Onco.Data[is.na(Onco.Data)] <- ''
#------------------------------------------------------------------------------#
Function <- read_excel('Figure5A.xlsx',
                       sheet='Enrichment')

Function <- as.data.frame(Function)
Function <- Function[grepl('Summary', Function[,1]),]

hemopoiesis <- Function[4, 8]
hemopoiesis <- unlist(strsplit(hemopoiesis, ','))

chromatin_organization <- Function[2,8]
chromatin_organization <- unlist(strsplit(chromatin_organization, ','))
chromatin_organization <- setdiff(chromatin_organization, hemopoiesis)

stem_cell_proliferation <- Function[13, 8]
stem_cell_proliferation <- unlist(strsplit(stem_cell_proliferation, ','))
stem_cell_proliferation <- Reduce(setdiff, list(chromatin_organization, hemopoiesis, stem_cell_proliferation))

immune_system_development <- Function[18, 8]
immune_system_development <- unlist(strsplit(immune_system_development, ','))
immune_system_development <- Reduce(setdiff, list(immune_system_development, chromatin_organization, hemopoiesis, stem_cell_proliferation))

Gene <- rownames(Onco.Data)

hemopoiesis <- Onco.Data %>%
  dplyr::filter(Gene %in% hemopoiesis)

chromatin_organization <- Onco.Data %>%
  dplyr::filter(Gene %in% chromatin_organization)

stem_cell_proliferation <- Onco.Data %>%
  dplyr::filter(Gene %in% stem_cell_proliferation)

immune_system_development <- Onco.Data %>%
  dplyr::filter(Gene %in% immune_system_development)

Onco.func.Data <- rbind(hemopoiesis, chromatin_organization, stem_cell_proliferation, immune_system_development)
#---------------------------------------------------------------------------------------------------#
#Annotation
Annotation <- read_excel('')
Annotation <- as.data.frame(Annotation)
Group = Annotation$Diagnosis
Sex = Annotation$Sex
col_fun = colorRamp2(c(66,84), c("seashell", "seashell4"))
ha = HeatmapAnnotation(Diagnosis = Annotation$Diagnosis,
                       Sex = Annotation$Sex,
                       HMA = Annotation$HMA,
                       Age = Annotation$Age,
                       col = list(Diagnosis = c("CR" = "forestgreen",
                                                "NR" = "firebrick1"),
                                  Sex = c('F'='lightcoral','M'='steelblue1'), 
                                  HMA = c('Decitabine' = 'goldenrod1'),
                                  Age = col_fun),
                       show_annotation_name = T,
                       annotation_label=c('HMA response','Sex','Treatment','Age'),
                       annotation_name_gp=gpar(fontsize = 10),
                       annotation_legend_param=list(Diagnosis=list(title='HMA response'),
                                                    Sex=list(title='Sex'),
                                                    HMA = list(title='HMA'),
                                                    Age = list(title='Age',direction = "horizontal")))

ha1 = HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = F,
                                                                height = unit(3, "cm")))


#Sample order
sample_order = colnames(Onco.Data)

#color지정
col.Onco = c("Missense" = "forestgreen",
             "Nonsense" = "red",
             "Frameshift" = "orange",
             "Insertion" = "dodgerblue4",
             "Deletion" = "black",
             "Duplication" = "mediumpurple3",
             "Copy_number_deletion" = "royalblue1",
             "Copy_number_amplification" = "deeppink")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  # big green
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Missense"], col = NA))
  },
  # big red
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w*0.93, h*0.33, 
              gp = gpar(fill = col.Onco["Nonsense"], col = NA))
  },
  # big orange
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9,  
              gp = gpar(fill = col.Onco["Frameshift"], col = NA))
  },
  # big black
  Insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Insertion"], col = NA))
  },
  # big yellow
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Deletion"], col = NA))
  },
  # big purple
  Duplication = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Duplication"], col = NA))
  },
  # big purple
  Copy_number_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Copy_number_deletion"], col = NA))
  },
  # big purple
  Copy_number_amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.948, h*0.9, 
              gp = gpar(fill = col.Onco["Copy_number_amplification"], col = NA))
  }
)

column_title = "Altered in 28 of 28 AML samples"
heatmap_legend_param = list(title = "Alternations", at = c("Missense",
                                                           "Nonsense",
                                                           "Frameshift",
                                                           'Insertion',
                                                           "Deletion",
                                                           "Duplication",
                                                           "Copy_number_deletion",
                                                           "Copy_number_amplification"), 
                            labels = c("Missense",
                                       "Nonsense",
                                       'Frameshift',
                                       "Insertion",
                                       "Deletion",
                                       "Duplication",
                                       "Copy number deletion", 
                                       "Copy number amplification"),
                            direction = "horizontal")
#------------------------------------------------------------------------------#
a <- oncoPrint(Onco.func.Data,
               alter_fun = alter_fun, col = col.Onco,
               column_order = sample_order,
               column_title = column_title, heatmap_legend_param = heatmap_legend_param,
               pct_side = "right", row_names_side = "left",
               show_column_names = F,
               remove_empty_columns = FALSE, remove_empty_rows = TRUE,
               alter_fun_is_vectorized = FALSE,
               bottom_annotation = ha,
               top_annotation = ha1,
               row_names_gp = gpar(fontsize = 7, fontface = "italic"),
               pct_gp = gpar(fontsize = 8))
draw(a, heatmap_legend_side = "right",
     annotation_legend_side = "right", merge_legend = TRUE,
     row_split = rep(c("Hematopoiesis", "Chromatin organization", "Stem cell\nproliferation", "Immune\nsystem\ndevelopment"),
                     c(dim(hemopoiesis)[1], dim(chromatin_organization)[1], dim(stem_cell_proliferation)[1], dim(immune_system_development)[1])),
     row_gap = unit(3, "mm"))
#------------------------------------------------------------------------------#