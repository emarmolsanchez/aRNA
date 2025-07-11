#######################################
#                                     #
#      Figure 2C. Mammoth aRNA        #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(showtext)
showtext_auto()
female = intToUtf8(9792)
male = intToUtf8(9794)

setwd("") #Insert path to Figure auxiliary files



Table_sex = read.table("Figure_2C_Table.txt", header=T)
Table_sex$ID = factor(Table_sex$ID, levels=c("Mammoth1", "Mammoth2", "Mammoth3", "Mammoth4", "Mammoth5", 
                                             "Mammoth6", "Mammoth7", "Mammoth8", "Mammoth9", "Mammoth10"))





setwd("") #Insert path to output Figure PDFs

#Figure 2C
pdf("aRNA_Mammoth_Figure_2C.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(Table_sex, aes(x = ID, y = Ratio)) + 
  geom_rect(data=NULL,aes(xmin=0.4,xmax=10.6,ymin=0.8,ymax=1), fill="grey95") +
  geom_rect(data=NULL,aes(xmin=0.4,xmax=10.6,ymin=0.4,ymax=0.6), fill="grey95") +
  geom_point(size=2.5, color= "#4D4DA1", alpha=0.7) +
  
  
  theme_classic() + 
  
  xlab("") + ylab("ChrX vs. Autosomal depth Ratio") +
  scale_y_continuous(expand=c(0,0), limits = c(0, 1.1), breaks = seq(0, 1.1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_x_discrete(labels = c("Mammoth 1", "Mammoth 2", "Mammoth 3", "Mammoth 4", "Mammoth 5", 
                              "Mammoth 6", "Mammoth 7", "Mammoth 8", "Mammoth 9", "Mammoth 10")) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(size=6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=-0.3, l=0.5, r=1, "cm")) +
  
  ggtitle("Genomic sex estimation") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold")) +
  
  annotate("text", x = 11, y = 0.9, label = female, size = 8) +
  geom_hline(yintercept = 1, linetype='twodash', color='grey60', linewidth=0.4) +
  geom_hline(yintercept = 0.8, linetype='twodash', color='grey60', linewidth=0.4) +
  annotate("text", x = 11, y = 0.5, label = male, size = 8, angle=315) +
  geom_hline(yintercept = 0.6, linetype='twodash', color='grey60', linewidth=0.4) +
  geom_hline(yintercept = 0.4, linetype='twodash', color='grey60', linewidth=0.4) +
  coord_cartesian(xlim=c(1,10), clip = "off")

dev.off()
