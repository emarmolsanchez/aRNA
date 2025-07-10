#######################################
#                                     #
#     Figures S3C-F. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



Exon_exon = read.table("Figure_S3C_Table.txt")
Exon_exon = subset(Exon_exon, Exon_exon$V4>=0.05)

GC1 = read.table("Figure_S3DE_Table.txt")
GC2 = read.table("Figure_S3F_Table.txt")





setwd("") #Insert path to output Figure PDFs

#Figure S3C
pdf("aRNA_Mammoth_Figure_S3C.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(Exon_exon, aes(x=log10(V2), y=log10(V3))) +
  geom_point(size=2, color="#5D987B", alpha=0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", fill="grey80") +
  
  theme_classic() +
  
  xlab(expression(paste("log" ["10"], "(Exon-exon spanning sequence count)"))) + 
  ylab(expression(paste("log" ["10"], "(Sequence count)"))) +
  scale_x_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5),
                     labels = c("0","0.5","1","1.5","2","2.5")) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5),
                     labels = c("0","0.5","1","1.5","2","2.5","3","3.5")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_text(label=expression(paste(italic("r "), "= 0.771")), x=0.2, y=3, color = "black", size=2.5) +
  
  ggtitle("Exon-exon spanning aRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure S3D
pdf("aRNA_Mammoth_Figure_S3D.pdf", height = 50/25.4, width = 55/25.4)
ggplot(GC1, aes(x = V4, y=log10(V2+1))) +
  geom_point(size=1.5, color="#5D987B", alpha=0.7) +
  
  theme_classic() +
  
  xlab("GC Content (%)") + 
  ylab(expression(paste("log" ["10"], "(Sequence count)"))) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = c("0","25","50","75","100")) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1),
                     labels = c("0", "1", "2", "3", "4", "5")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Sequence count bias") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure S3E
pdf("aRNA_Mammoth_Figure_S3E.pdf", height = 50/25.4, width = 55/25.4)
ggplot(GC1, aes(x = V3, y= V4)) +
  geom_point(size=1.5, color="#5D987B", alpha=0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", fill="grey80") +
  
  theme_classic() +
  
  ylab("GC Content (%)") + 
  xlab("Coverage (%)") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = c("0","25","50","75","100")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = c("0", "25", "50", "75", "100")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_text(label=expression(paste(italic("r "), "= 0.147")), x=0.2, y=90, color = "black", size=2.5) +
  
  ggtitle("Coverage bias") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()




#Figure S3F
pdf("aRNA_Mammoth_Figure_S3F.pdf", height = 50/25.4, width = 55/25.4)
ggplot(GC2, aes(x = V2, y= V3)) +
  geom_point(size=1.5, color="#5D987B", alpha=0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", fill="grey80") +
  
  theme_classic() +
  
  xlab("Covered regions (%)") + 
  ylab("Uncovered regions (%)") +
  scale_x_continuous(limits = c(25, 75), breaks = seq(25, 75, 10),
                     labels = c("25","35", "45", "55", "65", "75")) +
  scale_y_continuous(limits = c(25, 75), breaks = seq(25, 75, 10),
                     labels = c("25","35", "45", "55", "65", "75")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_text(label=expression(paste(italic("r "), "= 0.568")), x=35, y=70, color = "black", size=2.5) +
  
  ggtitle("GC content bias") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
