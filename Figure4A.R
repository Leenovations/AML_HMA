library(ggplot2)
library(ggrepel)
#------------------------------------------------------------------------------#
Table <- read.table("Figure4A.txt", header=T)
Table <- subset(Table, log2(p.value) != 0)
# Table$"DEG" <- ifelse(Table$p.value<0.05 & Table$m.value > 2, "NR Up-reg",
#                       ifelse(Table$p.value<0.05 & Table$m.value < -2, "CR Up-reg","None"))
TLR_subset <- Table[grep("^TLR", Table$GeneSymbol), ]
#------------------------------------------------------------------------------#
volcano <- ggplot(Table, aes(x = m.value, y= -log(p.value,10), col = CustomedDEG)) +
  ggtitle("") +
  geom_point(alpha=0.8, size=3) + 
  scale_color_manual(values=c("CR"= "forestgreen", "NR"="firebrick1", "None"="#999999")) +
  labs(x=(log[2]~fc), y=(-log[10]~p~value)) +
  xlim(-11, 11) + 
  ylim(0, 13) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5),
        axis.title.y=element_text(size=13, margin = margin(r=10), color='black'),
        axis.title.x=element_text(size=13, margin = margin(t=10), color='black'),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom") +
  geom_text_repel(data = subset(TLR_subset, CustomedDEG == 'NR'),
                  aes(label = GeneSymbol),
                  size = 5,
                  box.padding = unit(1, "lines"),
                  fontface = "italic",
                  nudge_x = 2.3,                        
                  nudge_y = 7,               
                  show.legend = FALSE,
                  min.segment.length = unit(0.2, "lines"),
                  point.padding = unit(1, "lines"))
#------------------------------------------------------------------------------#