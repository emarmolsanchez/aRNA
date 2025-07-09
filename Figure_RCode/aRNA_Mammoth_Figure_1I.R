#######################################
#                                     #
#       Figure 1I. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Counts = read.delim("Figure1I_Table_Exonic-Intronic-Intergenic_Counts_23nt.txt", header=T)
Counts$Type = factor(Counts$Type, levels=c("Exonic", "Intronic", "Intergenic"))
Counts$ID = factor(Counts$ID, levels = c("Mammoth1", "Mammoth4", "Mammoth10"))
Counts$Color = factor(Counts$Color, levels = c("#5D987B", "mediumpurple", "orange2"))





setwd("") #Insert path to output Figure PDFs

#Figure 1I
pdf("aRNA_Mammoth_Figure_1I.pdf", height = 40/25.4, width = 45/25.4)
ggplot(Counts, aes(x = ID, y = Count, fill = Type, Color=Type)) +
  geom_bar(stat = "identity", position="fill", width=0.6, alpha=0.7, color=rep(c("#5D987B", "mediumpurple", "orange2"),3)) +
  
  theme_classic() +
  
  xlab("") + ylab("Mapped aRNAs (%)") +
  scale_fill_manual(values=c("Exonic" = "#5D987B", "Intronic" = "mediumpurple", "Intergenic" = "orange2"), labels=c("Exonic", "Intronic", "Intergenic")) +
  scale_color_manual(values=c("#5D987B", "mediumpurple", "orange2")) +
  
  scale_y_continuous(expand=c(0,0), labels = scales::label_number(scale = 100)) +
  scale_x_discrete(labels = c("Mammoth 1", "Mammoth 4", "Mammoth 10")) +
  
  guides(fill = guide_legend(override.aes = list(color =c("#5D987B", "mediumpurple", "orange2")))) +
  
  theme(legend.key.size = unit(0.65, 'lines'),
        legend.margin = margin(l=-5),
        legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=-0.25, l=0.25, r=0.25, "cm")) +
  
  ggtitle("Genomic alignment") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()