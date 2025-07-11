#######################################
#                                     #
#       Figure 1H. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Length_Count = read.table("Figure_1H_Table.txt", header=T)
Length_Count$ID = factor(Length_Count$ID, levels=c("Mammoth1", "Mammoth4", "Mammoth10"))
Length_Count$Type = factor(Length_Count$Type, levels=c("Exonic", "Intronic", "Intergenic"))





setwd("") #Insert path to output Figure PDFs

#Figure 1H
pdf("aRNA_Mammoth_Figure_1H.pdf", height = 40/25.4, width = 75/25.4)
ggplot(Length_Count, aes(x=Length, y=log10(Count), group=Type, fill=Type, color=Type)) +
  geom_line(linewidth=1, alpha=0.7) +
  facet_wrap( ~ ID, ncol=3, strip.position="top") +
  
  theme_classic() +
  
  scale_fill_manual(values=c("Exonic" = "#5D987B", "Intronic" = "mediumpurple", "Intergenic" = "orange2"), labels=c("Exonic", "Intronic", "Intergenic")) +
  scale_color_manual(values=c("#5D987B", "mediumpurple", "orange2")) +
  
  xlab("Sequence length (nt)") +
  ylab(expression(paste("log" ["10"], "(Sequence count)"))) +
  
  scale_x_continuous(limits = c(18, 25), breaks = seq(18, 25, 2),
                     labels = c("18", "20", "22", "24")) +
  scale_y_continuous(limits = c(2, 6), breaks = seq(2, 6, 1),
                     labels = c("2", "3", "4", "5", "6")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = margin(t = 0.75, r = 0.25, b = -0.2, l = 0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5, b=10),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  ggtitle("Mammoth 1        Mammoth 4        Mammoth 10") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
