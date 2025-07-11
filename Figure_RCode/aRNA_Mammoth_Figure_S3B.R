#######################################
#                                     #
#      Figure S3B. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggridges)


setwd("") #Insert path to Figure auxiliary files



Table_Transcripts = read.delim("Figure_S3B_Table.txt", header=T)
Table_Transcripts$ID = factor(Table_Transcripts$ID, levels=c("Mammoth 1", "Mammoth 4", "Mammoth 10"))
Table_Transcripts$Type = factor(Table_Transcripts$Type, levels=c("Protein-coding mRNAs", "rRNAs", "tRNAs", "snRNAs",
                                                                 "snoRNAs", "lncRNAs"))





setwd("") #Insert path to output Figure PDFs

#Figure S3B
pdf("aRNA_Mammoth_Figure_S3B.pdf", height = 150/25.4, width = 82.5/25.4)
ggplot(Table_Transcripts, aes(x = Length, y = ID, fill = ID, color= ID)) +
  geom_density_ridges(scale=4, bandwidth = 0.5, alpha=0.7) +
  facet_wrap( ~ Type, ncol=1, strip.position="top") +
  
  theme_classic() +
  
  scale_fill_manual(values=c("Mammoth 1" = "#DB6443", "Mammoth 4" = "plum", "Mammoth 10" = "#F0C57F"), 
                    labels=c("Mammoth 1", "Mammoth 4", "Mammoth 10")) +
  scale_color_manual(values=c("#DB6443", "plum", "#F0C57F")) +
  
  xlab("Sequence length (nt)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line.x = element_line(linewidth = 0.3),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0, l = 0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5, b=10),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2)))

dev.off()
