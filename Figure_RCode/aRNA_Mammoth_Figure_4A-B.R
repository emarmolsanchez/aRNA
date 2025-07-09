#######################################
#                                     #
#      Figures 4A-B. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Count_Cov_PC = read.table("Figure4A_Table.txt",
                          header=F)
colnames(Count_Cov_PC) = c("ID","Gene","Count", "Coverage")

Count_Cov_NC = read.table("Figure4B_Table.txt",
                          header=F)
colnames(Count_Cov_NC) = c("ID","Gene", "Type", "Count", "Coverage")

Count_Cov_PC = subset(Count_Cov_PC, Count_Cov_PC$Coverage>=0.05)
Count_Cov_NC = subset(Count_Cov_NC, Count_Cov_NC$Coverage>=0.05)





setwd("") #Insert path to output Figure PDFs

#Figure 4A-B
pdf("aRNA_Mammoth_Figure_4A.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(Count_Cov_PC, aes(x=Coverage, y=log10(Count+1))) +
  geom_point(size=1.5, color="#DB6443", alpha=0.5) +
  
  theme_classic() +
  
  xlab("Coverage (%)") + ylab(expression(paste("log" ["10"] ,"(Sequence count)"))) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1),
                     labels = c("0", "1", "2", "3", "4")) +
  scale_x_continuous(expand=c(0.015, 0.015), breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "25", "50", "75", "100"), limits = c(0,1)) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 1 protein-coding mRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



pdf("aRNA_Mammotn_Figure_4B.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(Count_Cov_NC, aes(x=Coverage, y=log10(Count+1))) +
  geom_point(size=1.5, color= "#DB6443", alpha=0.5) +
  
  theme_classic() +
  
  xlab("Coverage (%)") + ylab(expression(paste("log" ["10"] ,"(Sequence count)"))) +
  scale_y_continuous(limits = c(0, 4.2), breaks = seq(0, 4.2, 1),
                     labels = c("0", "1", "2", "3", "4")) +
  scale_x_continuous(expand=c(0.015, 0.015), breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "25", "50", "75", "100"), limits = c(0,1)) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mammoth 1 noncoding RNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
