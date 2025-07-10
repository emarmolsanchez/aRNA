#######################################
#                                     #
#      Figure S2D. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggpubr)

setwd("") #Insert path to Figure auxiliary files



PC = read.table("Figure_S2D_Table1.txt", header=T)
PC$ID = factor(PC$ID, levels=c("Mammoth10","Mammoth9","Mammoth8","Mammoth7","Mammoth6","Mammoth5",
                               "Mammoth4","Mammoth3","Mammoth2","Mammoth1"))
PC$Class = factor(PC$Class, levels=c("Strong", "Weak"))

NC = read.table("Figure_S2D_Table2.txt", header=T)
NC$ID = factor(NC$ID, levels=c("Mammoth10","Mammoth9","Mammoth8","Mammoth7","Mammoth6","Mammoth5",
                               "Mammoth4","Mammoth3","Mammoth2","Mammoth1"))
NC$Class = factor(NC$Class, levels=c("Strong", "Weak"))





setwd("") #Insert path to output Figure PDFs

#Figure S2D
PC_plot = ggplot(PC, aes(x = Count, y= ID, fill = Class, color=Class)) + 
  geom_bar(stat="identity", alpha=0.7, width = 0.5) +
  
  theme_classic() + 
  
  xlab("") + ylab("") +
  scale_fill_manual(values=c("Strong" = "#5D987B", "Weak" = "#C0D8CC"), labels=c("Coverage ≥5%", "Coverage <5%")) +
  scale_color_manual(values=c("#5D987B","#C0D8CC"), labels=c("Coverage ≥5%", "Coverage <5%")) +
  
  scale_x_continuous(expand=c(0,0), limits = c(0, 20000), breaks = seq(0, 20000, 5000),
                     labels = c("0", "5K", "10K", "15K", "20K")) +
  scale_y_discrete(labels = c("Mammoth 10","Mammoth 9","Mammoth 8","Mammoth 7","Mammoth 6","Mammoth 5",
                              "Mammoth 4","Mammoth 3","Mammoth 2","Mammoth 1")) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5, b=10),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.x =element_text(size=7, margin = margin(t = 5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=0, l=-0.3, r=0.5, "cm")) +
  ggtitle("Protein-coding") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))


NC_plot = ggplot(NC, aes(x = Count, y= ID, fill = Class, color=Class)) + 
  geom_bar(stat="identity", alpha=0.7, width = 0.5) +
  
  theme_classic() + 
  
  xlab("") + ylab("") +
  scale_fill_manual(values=c("Strong" = "#5D987B", "Weak" = "#C0D8CC"), labels=c("Coverage ≥5%", "Coverage <5%")) +
  scale_color_manual(values=c("#5D987B","#C0D8CC"), labels=c("Coverage ≥5%", "Coverage <5%")) +
  
  scale_x_continuous(expand=c(0,0), limits = c(0, 120000), breaks = seq(0, 120000, 30000),
                     labels = c("0", "30K", "60K", "90K", "120K")) +
  scale_y_discrete(labels = c("Mammoth 10","Mammoth 9","Mammoth 8","Mammoth 7","Mammoth 6","Mammoth 5",
                              "Mammoth 4","Mammoth 3","Mammoth 2","Mammoth 1")) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5, b=10),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.x =element_text(size=7, margin = margin(t = 5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0, b=0.25, l=-0.3, r=0.5, "cm")) +
  ggtitle("Noncoding") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))



pdf("aRNA_Mammoth_Figure_S2D.pdf", height = 75/25.4, width = 82.5/25.4)
ggarrange(PC_plot, NC_plot, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

dev.off()
