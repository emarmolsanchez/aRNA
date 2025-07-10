#######################################
#                                     #
#      Figure S3A. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggridges)


setwd("") #Insert path to Figure auxiliary files



Lengths = read.table("Figure_S3A_Table.txt", header=T)
Lengths$IDs = factor(Lengths$IDs, levels=c("Mammoth 1", "Mammoth 4", "Mammoth 10"))
Lengths$Type = factor(Lengths$Type, levels=c("Exonic", "Intronic", "Intergenic"))





setwd("") #Insert path to output Figure PDFs

#Figure S3A
pdf("aRNA_Mammoth_Figure_S3A.pdf", height = 100/25.4, width = 82.5/25.4)
ggplot(Lengths, aes(x = V1, y = Type, fill = Type, color= Type)) +
  geom_density_ridges(scale=4, bandwidth = 0.5, alpha=0.7) +
  facet_wrap( ~ IDs, ncol=1, strip.position="top") +
  
  theme_classic() +
  
  scale_fill_manual(values=c("Exonic" = "#5D987B", "Intronic" = "mediumpurple", "Intergenic" = "orange2"), labels=c("Exonic", "Intronic", "Intergenic")) +
  scale_color_manual(values=c("#5D987B", "mediumpurple", "orange2")) +
  
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
