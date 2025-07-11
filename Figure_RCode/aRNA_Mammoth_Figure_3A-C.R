#######################################
#                                     #
#      Figures 3A-C. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



Sample1 = read.table("Figure_3A_Table.txt", header=T)
Sample4 = read.table("Figure_3B_Table.txt", header=T)
Sample10 = read.table("Figure_3C_Table.txt", header=T)





setwd("") #Insert path to output Figure PDFs

#Figure 3A1
pdf("aRNA_Mammoth_Figure_3A1.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample1) +
  geom_line(aes(x=X1, y=CT5), color="red3", linewidth=0.3) +
  geom_line(aes(x=X1, y=CT_CpG_5), color="red3", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5), labels=c("1","5","10","15","20", "25", "30")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 5' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 1 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()


#Figure 3A2
pdf("aRNA_Mammoth_Figure_3A2.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample1) +
  geom_line(aes(x=X1, y=GA3), color="cornflowerblue", linewidth=0.3) +
  geom_line(aes(x=X1, y=GA_CpG_3), color="cornflowerblue", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_reverse(limits = c(30, 0), breaks = seq(30, 0, -5), labels=c("30","25","20","15","10", "5", "0")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 3' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 1 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure 3B1
pdf("aRNA_Mammoth_Figure_3B1.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample4) +
  geom_line(aes(x=X1, y=CT5), color="red3", linewidth=0.3) +
  geom_line(aes(x=X1, y=CT_CpG_5), color="red3", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5), labels=c("1","5","10","15","20", "25", "30")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 5' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 4 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()


#Figure 3B2
pdf("aRNA_Mammoth_Figure_3B2.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample4) +
  geom_line(aes(x=X1, y=GA3), color="cornflowerblue", linewidth=0.3) +
  geom_line(aes(x=X1, y=GA_CpG_3), color="cornflowerblue", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_reverse(limits = c(30, 0), breaks = seq(30, 0, -5), labels=c("30","25","20","15","10", "5", "0")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 3' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 4 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure 3C1
pdf("aRNA_Mammoth_Figure_3C1.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample10) +
  geom_line(aes(x=X1, y=CT5), color="red3", linewidth=0.3) +
  geom_line(aes(x=X1, y=CT_CpG_5), color="red3", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5), labels=c("1","5","10","15","20", "25", "30")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 5' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 10 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()


#Figure 3C2
pdf("aRNA_Mammoth_Figure_3C2.pdf", height = 40/25.4, width = 41.25/25.4)
ggplot(data=Sample10) +
  geom_line(aes(x=X1, y=GA3), color="cornflowerblue", linewidth=0.3) +
  geom_line(aes(x=X1, y=GA_CpG_3), color="cornflowerblue", linewidth=0.3, linetype="twodash") +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,0.1) +
  
  scale_x_reverse(limits = c(30, 0), breaks = seq(30, 0, -5), labels=c("30","25","20","15","10", "5", "0")) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02), labels=c("0", "2", "4", "6", "8", "10")) +
  xlab("Distance from 3' end") +
  ylab("Mismatch frequency (%)") +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 10 (aDNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
