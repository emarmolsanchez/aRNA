#######################################
#                                     #
#       Figure 1J. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Enrich = read.table("Figure_1J_Table.txt", header=T)
Enrich$Type = factor(Enrich$Type, levels=c("Exonic", "Intronic"))
Enrich$ID = factor(Enrich$ID, levels=c("Mammoth1", "Mammoth4", "Mammoth10"))





setwd("") #Insert path to output Figure PDFs

#Figure 1J
pdf("aRNA_Mammoth_Figure_1J.pdf", height = 40/25.4, width = 45/25.4)
ggplot(Enrich, aes(x=ID, y=log2FC, fill=Type, color=Type)) + 
  geom_bar(stat="identity", width=0.7, position=position_dodge(0.8), alpha=0.7) +
  
  theme_classic() + 
  
  xlab("") + ylab(expression(paste("log" ["2"], "(Fold change)"))) +
  scale_y_continuous(limits = c(-3, 6), breaks = seq(-3, 6, 3), labels=c("-3", "0", "3", "6")) +
  scale_x_discrete(labels = c("Mammoth 1", "Mammoth 4", "Mammoth 10")) +
  scale_fill_manual(values=c("Exonic" = "#5D987B", "Intronic" = "mediumpurple"), labels=c("Exonic", "Intronic")) +
  scale_color_manual(values=c("#5D987B", "mediumpurple"), labels=c("Exonic", "Intronic")) +
  
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size=7, margin = margin(l = 2)),
        legend.key.size = unit(0.65, 'lines'),
        legend.margin = margin(l=-5)) +
  
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=-0.25, l=0.25, r=0.2, "cm")) +
  
  geom_hline(yintercept=0, linewidth=0.5, linetype='solid', color='grey60') +
  
  ggtitle("Coding enrichment") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
