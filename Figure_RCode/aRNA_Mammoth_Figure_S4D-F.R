#######################################
#                                     #
#     Figures S4D-F. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



Damage1 = read.table("Figure_S4D_Table.txt", header=T)
Damage4 = read.table("Figure_S4E_Table.txt", header=T)
Damage10 = read.table("Figure_S4F_Table.txt", header=T)





setwd("") #Insert path to output Figure PDFs

#Figure S4D
pdf("aRNA_Mammoth_Figure_S4D.pdf", height = 40/25.4, width = 82.5/25.4)
ggplot(data=Damage1) +
  geom_line(aes(x=z, y=CT5), color="red3", linewidth=0.5) +
  
  theme_classic() +
  
  ylab("% Mismatch frequency (%)") +
  xlab("Distance from 5' read end (nt)") +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05),
                     labels = c("0", "5", "10", "15", "20", "25", "30")) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, 10),
                     labels = c("0", "10", "20", "30", "40", "50", "60", "70")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_segment(aes(x = 0, y = 0.25, xend = 6, yend = 0.25), colour="red3", linewidth=0.75) +
  geom_text(label="C>U", x=10, y=0.25, color = "black", size=2) +
  geom_text(label="Avg. C>U deaminated sites = 1.97%", x=18.5, y=0.3, color = "black", size=2) +
  
  ggtitle("Mammoth 1 (aRNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure S4E
pdf("aRNA_Mammoth_Figure_S4E.pdf", height = 40/25.4, width = 82.5/25.4)
ggplot(data=Damage4) +
  geom_line(aes(x=z, y=CT5), color="red3", linewidth=0.5) +
  
  theme_classic() +
  
  ylab("% Mismatch frequency (%)") +
  xlab("Distance from 5' read end (nt)") +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05),
                     labels = c("0", "5", "10", "15", "20", "25", "30"))  +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, 10),
                     labels = c("0", "10", "20", "30", "40", "50", "60", "70")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_segment(aes(x = 0, y = 0.25, xend = 6, yend = 0.25), colour="red3", linewidth=0.75) +
  geom_text(label="C>U", x=10, y=0.25, color = "black", size=2) +
  geom_text(label="Avg. C>U deaminated sites = 1.04%", x=18.5, y=0.3, color = "black", size=2) +
  
  ggtitle("Mammoth 4 (aRNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure S4F
pdf("aRNA_Mammoth_Figure_S4F.pdf", height = 40/25.4, width = 82.5/25.4)
ggplot(data=Damage10) +
  geom_line(aes(x=z, y=CT5), color="red3", linewidth=0.5) +
  
  theme_classic() +
  
  ylab("% Mismatch frequency (%)") +
  xlab("Distance from 5' read end (nt)") +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05),
                     labels = c("0", "5", "10", "15", "20", "25", "30"))  +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, 10),
                     labels = c("0", "10", "20", "30", "40", "50", "60", "70")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_segment(aes(x = 0, y = 0.25, xend = 6, yend = 0.25), colour="red3", linewidth=0.75) +
  geom_text(label="C>U", x=10, y=0.25, color = "black", size=2) +
  geom_text(label="Avg. C>U deaminated sites = 2.02%", x=18.5, y=0.3, color = "black", size=2) +
  
  
  ggtitle("Mammoth 10 (aRNA)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
