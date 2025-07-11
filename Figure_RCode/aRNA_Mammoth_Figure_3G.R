#######################################
#                                     #
#       Figures 3G. Mammoth aRNA      #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



prob_S1 = read.table("Figure_3G_Table1.txt")
prob_S4 = read.table("Figure_3G_Table2.txt")





setwd("") #Insert path to output Figure PDFs

#Figure 3G1
pdf("aRNA_Mammoth_Figure_3G1.pdf", height = 40/25.4, width = 82.5/25.4)
ggplot(prob_S1, aes(x = V2)) +
  stat_bin(geom = "bar", bins=200, fill="red3") +
  
  theme_classic() +
  
  xlab("Endogenous probability (%)") + ylab("Read count") +
  scale_y_continuous(limits = c(0,40000), breaks = seq(0, 40000, 10000),
                     labels = c("0", "10K", "20K", "30K", "40K")) +
  scale_x_continuous(limits = c(0.9985,1), breaks = seq(0.9985, 1, 0.0005),
                     labels = c("99.85", "99.90", "99.95", "100")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  theme(legend.position = "none", 
        legend.title = element_blank()) +
  
  ggtitle("Mammoth 1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure 3G2
pdf("aRNA_Mammoth_Figure_3G2.pdf", height = 40/25.4, width = 82.5/25.4)
ggplot(prob_S4, aes(x = V2)) +
  stat_bin(geom = "bar", bins=100, fill="red3") +
  
  theme_classic() +
  
  xlab("Endogenous probability (%)") + ylab("Read count") +
  scale_y_continuous(limits = c(0,15000), breaks = seq(0, 15000, 5000),
                     labels = c("0", "5K", "10K", "15K")) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2),
                     labels = c("0","20", "40", "60", "80","100")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  theme(legend.position = "none", 
        legend.title = element_blank()) +
  
  ggtitle("Mammoth 4") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
