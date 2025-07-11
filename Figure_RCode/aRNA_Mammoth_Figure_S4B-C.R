#######################################
#                                     #
#     Figures S4B-C. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



S1_mirna = read.table("Figure_S4B_Table.txt")
colnames(S1_mirna) = c("Length")

Sample = c(rep("Mammoth 1", nrow(S1_mirna)))
Table_hist = as.data.frame(cbind(S1_mirna, Sample))
Table_hist$Sample = factor(Table_hist$Sample, levels=c("Mammoth 1"))

Damage1 = read.table("Figure_S4C_Table.txt", header=T)





setwd("") #Insert path to output Figure PDFs

#Figure S4B
pdf("aRNA_Mammoth_Figure_S4B.pdf", height = 40/25.4, width = 75/25.4)
ggplot(Table_hist, aes(x = Length, y = after_stat(density), fill = Sample, color=Sample)) +
  geom_density(alpha = 0.6, adjust = 2, linewidth=0.2) +
  
  xlab("Sequence length (nt)") +
  ylab("Density") +
  
  scale_x_continuous(limits = c(15, 30), breaks = seq(15, 30, 5), labels=c("15","20","25","30")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  #coord_cartesian(xlim = c(min(Table_hist$Length) + 0, max(Table_hist$Length) + 0)) +
  
  scale_fill_manual(values = c("#DB6443")) +
  
  theme_classic() +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x =element_text(size=7, margin = margin(t = 5)),
        axis.line.x = element_line(linewidth = 0.3),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  ggtitle("Mammoth 1 microRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure S4C
pdf("aRNA_Mammoth_Figure_S4C.pdf", height = 40/25.4, width = 75/25.4)
ggplot(data=Damage1) +
  geom_line(aes(x=z, y=CT5), color="red3", linewidth=0.5) +
  
  theme_classic() +
  
  ylab("Mismatch frequency (%)") +
  xlab("Distance from 5' read end (nt)") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     labels = c("0", "10", "20", "30", "40", "50")) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     labels = c("0", "5", "10", "15", "20")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  geom_segment(aes(x = 0, y = 0.4, xend = 2, yend = 0.4), colour="red3", linewidth=0.75) +
  geom_text(label="C>U", x=3, y=0.4, color = "black", size=2) +
  geom_text(label="Avg. C>U deaminated sites = 6.42%", x=5.8, y=0.5, color = "black", size=2) +
  
  ggtitle("Mammoth 1 (microRNAs)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
