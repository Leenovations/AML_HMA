library(ggbiplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(RColorBrewer)
#------------------------------------------------------------------------------#
GOGene_list <- read.table('Figure2E.txt', 
                          sep='\t', 
                          header=T)
HMM_list <- c('K562', 'GM12878', 'H1hESC', 'HMEC', 'HSMM', 'HUVEC', 'HepG2', 'NHEK', 'NHLF')
#------------------------------------------------------------------------------#
library(ggplot2) 
library(grid)
library(gridExtra) 
#------------------------------------------------------------------------------#
Manual_legend_value <- c(1,1,1,1,1,1,1,1,1,1,1,1)
HMM_Region <- c('Active Promoter',
                'Weak Promoter',
                'Poised Promoter',
                'Strong Enhancer',
                'Weak Enhancer',
                'Repetitive/CNV',
                'Heterochrom/lo',
                'Insulator',
                'Txn Transition',
                'Txn Elongation',
                'Weak Txn',
                'Repressed')

CHROMHMM <- data.frame('Region' = HMM_Region,
                       'value' = Manual_legend_value)

my_hist <- ggplot(CHROMHMM, aes(value, fill = Region)) + 
  geom_bar() +
  scale_fill_manual(values = c('Active Promoter'='firebrick1',
                               'Weak Promoter'='firebrick3',
                               'Poised Promoter'='purple',
                               'Strong Enhancer'='darkorange',
                               'Weak Enhancer'='#FFFF99',
                               'Repetitive/CNV'='gray88',
                               'Heterochrom/lo'='gray78',
                               'Insulator'='#1F78B4',
                               'Txn Transition'='#33A02C',
                               'Txn Elongation'='#33A06C',
                               'Weak Txn'='#B2DF8A',
                               'Repressed'='gray45'),
                    breaks = c('Active Promoter',
                               'Weak Promoter',
                               'Poised Promoter',
                               'Strong Enhancer',
                               'Weak Enhancer',
                               'Repetitive/CNV',
                               'Heterochrom/lo',
                               'Insulator',
                               'Txn Transition',
                               'Txn Elongation',
                               'Weak Txn',
                               'Repressed')) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'))

legend <- get_legend(my_hist)
#------------------------------------------------------------------------------#
for (gene in sort(GOGene_list$GeneSymbol)){
  Range <- subset(GOGene_list, GOGene_list$GeneSymbol==gene)
  Position_Start = Range$Start
  Position_End = Range$End
  #------------------------------------------------------------------------------#
  Input <- read.table(paste0('Figure2H/', gene, '.bin.meth.tsv'),
                      sep='\t',
                      header=T)
  Input <- cbind(Input[,1:3], Input[,4:6]/100)
  Input$Start <- ifelse(Input$Start < Position_Start, Position_Start, Input$Start)
  Input$End <- ifelse(Input$End > Position_End, Position_End, Input$End)
  #------------------------------------------------------------------------------#
  spline_CR <- abs(as.data.frame(spline(Input$End, Input$CR)))
  colnames(spline_CR) <- c('End','CR')
  spline_NR <- abs(as.data.frame(spline(Input$End, Input$NR)))
  colnames(spline_NR) <- c('End','NR')
  spline_Normal <- abs(as.data.frame(spline(Input$End, Input$Normal)))
  colnames(spline_Normal) <- c('End','Normal')
  #------------------------------------------------------------------------------#
  p1 <- ggplot() +
    geom_line(data = spline_NR,  aes(x = End, y = NR, colour='NR'), size=1) +
    geom_line(data = spline_Normal,  aes(x = End, y = Normal, colour='Normal'), size=1) +
    geom_line(data = spline_CR,  aes(x = End, y = CR, colour='CR'), size=1) +
    scale_colour_manual(values = c("CR"="forestgreen", "NR"="firebrick1", "Normal"="orange"),
                        breaks = c("CR", "NR", "Normal")) +
    xlim(Position_Start, Position_End) + 
    ylab('DNA\nmethylation\n (WGBS)') +
    # scale_y_continuous(breaks = c(0,1), labels = c(0,1)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=10, margin=margin(r=10)),
          axis.text=element_text(color="black"), 
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title = element_blank(),
          legend.position = "right")
  #------------------------------------------------------------------------------#
  ExonIntron <- read.table(paste0('Figure2H/', gene,'.ExonIntron.tsv'),
                           sep='\t', 
                           header=T)
  
  ExonIntron$Region <- gsub("Intron_\\d+", "Intron", ExonIntron$Region)
  ExonIntron$Region <- gsub("Exon_\\d+", "Exon", ExonIntron$Region)
  ExonIntron['Color'] = ifelse(ExonIntron$Region=='Exon', "black", "gray")
  ExonIntron['Size'] = ifelse(ExonIntron$Region=='Exon', 1, 0.5)
  NM_number <- unique(ExonIntron$NM)
  #------------------------------------------------------------------------------#
  p2 <- ggplot(ExonIntron) +
    geom_segment(aes(x = Start, xend = End,
                     y = 1, yend = 1,
                     colour = Color, size= Size)) +
    scale_color_manual(values = c("navy", 'navy'), 
                       labels = c("Exon", "Intron")) + 
    xlim(Position_Start, Position_End) + 
    ylim(1,1) +
    xlab('') + 
    ylab(NM_number) +
    ggtitle(paste0(gene, ' (', unique(ExonIntron$Chromosome), ':', format(min(ExonIntron$Start), big.mark = ","), '-', format(max(ExonIntron$End), big.mark = ","), ' ,', unique(ExonIntron$Strand), ' strand', ')')) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=10, hjust=0.5, face='italic'),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=10, margin=margin(r=10)),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.background = element_blank(),
          legend.position = "none",
          legend.key = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank())
  #------------------------------------------------------------------------------#
  DMR <- read.table(paste0('Figure2H/', gene,'.DMR.tsv'),
                    sep='\t', 
                    header=T)
  
  DMR$Start_methylkit <- ifelse(DMR$Start_methylkit < Position_Start, Position_Start, DMR$Start_methylkit)
  DMR$End_methylkit <- ifelse(DMR$End_methylkit > Position_End, Position_End, DMR$End_methylkit)
  
  DMR['Color'] = 'darkorchid2'
  DMR['Size'] = 1
  
  p4 <- ggplot(DMR) +
    geom_segment(aes(x = Start_methylkit, xend = End_methylkit,
                     y = 1, yend = 1,
                     colour = Color, size= Size)) +
    scale_color_manual(values = c("darkorchid2"), 
                       labels = c("TF")) + 
    xlim(Position_Start, Position_End) + 
    ylim(1,1) +
    xlab('') + 
    ylab('DMR') +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=10, hjust=0.5, face='italic'),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=10, margin=margin(r=10)),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.background = element_blank(),
          legend.position = "none",
          legend.key = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank())
  #------------------------------------------------------------------------------#
  if (file.exists(paste0('Figure2H/', gene,'.bin.CpG.tsv'))) {
    CpG <- read.table(paste0('Figure2H/', gene,'.bin.CpG.tsv'),
                      sep='\t', 
                      header=T)
    
    CpG$Start <- ifelse(CpG$Start < Position_Start, Position_Start, CpG$Start)
    CpG$End <- ifelse(CpG$End > Position_End, Position_End, CpG$End)
    
    CpG['Color'] = "goldenrod1"
    CpG['Size'] = 1
    
    p5 <- ggplot(CpG) +
      geom_segment(aes(x = Start, xend = End,
                       y = 1, yend = 1,
                       colour = Color, size= Size)) +
      scale_color_manual(values = c("goldenrod1"), 
                         labels = c("CpG_island")) + 
      xlim(Position_Start, Position_End) + 
      ylim(1,1) +
      xlab('') + 
      ylab('CpG island') +
      theme_bw() +
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=10, hjust=0.5, face='italic'),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=10, margin=margin(r=10)),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.background = element_blank(),
            legend.position = "none",
            legend.key = element_blank(),
            legend.text = element_blank(),
            legend.title = element_blank())
  }
  #------------------------------------------------------------------------------#
  order <- 6
  for (hmm in HMM_list){
    plot_names <- paste0("p", order)
    
    HMM <- read.table(paste0('Figure2H/', gene,'.bin.', hmm ,'.HMM.tsv'),
                      sep='\t', 
                      header=T)
    
    HMM$Start <- ifelse(HMM$Start < Position_Start, Position_Start, HMM$Start)
    HMM$End <- ifelse(HMM$End > Position_End, Position_End, HMM$End)
    
    HMM$Region <- gsub("_", " ", HMM$Region)
    HMM['Size'] = 1
    
    p <- ggplot() +
      geom_segment(data=HMM, aes(x = Start, xend = End,
                                 y = 1, yend = 1,
                                 colour = Region, 
                                 size= Size)) +
      scale_color_manual(values = c('Active Promoter'='firebrick1',
                                    'Weak Promoter'='firebrick3',
                                    'Poised Promoter'='purple',
                                    'Strong Enhancer'='darkorange',
                                    'Weak Enhancer'='#FFFF99',
                                    'Repetitive/CNV'='gray88',
                                    'Heterochrom/lo'='gray78',
                                    'Insulator'='#1F78B4',
                                    'Txn Transition'='#33A02C',
                                    'Txn Elongation'='#33A06C',
                                    'Weak Txn'='#B2DF8A',
                                    'Repressed'='gray45')) +
      # scale_color_manual(values = rev(brewer.pal(n = length(unique(HMM$Region)), name = "Paired"))) +
      guides(size = "none", color = guide_legend(override.aes = list(linewidth = 6))) +
      xlim(Position_Start, Position_End) + 
      ylim(1,1) +
      xlab('') + 
      ylab(paste0(hmm,'\nChromHMM')) +
             theme_bw() +
             guides(size = F) +
             theme(panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=10, margin=margin(r=10)),
                   axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.background = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "none")
    #--------------------------------------------------------------------#
      assign(plot_names, p)
      order <- order + 1
    }
}