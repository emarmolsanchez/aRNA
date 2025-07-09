#######################################
#                                     #
#     Figures 1B-D. Mammoth aRNA      #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Kraken = read.delim("Figure1BD_Table_KrakenUniq_Counts_Mammoth_Afrotheria_level.txt", header=T)
Kraken$ID = factor(Kraken$ID, levels=c("Mammoth1", "Mammoth2", "Mammoth3", "Mammoth4", "Mammoth5", 
                                       "Mammoth6", "Mammoth7", "Mammoth8", "Mammoth9", "Mammoth10"))
Kraken$Type = factor(Kraken$Type, levels=c("aDNA", "aRNA"))
Kraken_aDNA = subset(Kraken, Kraken$Type=="aDNA")
Kraken_aDNA$alpha = c(0.7, 0.2, 0.2, 0.7, 0.2, 0.2, 0.2, 0.2, 0.2, 0.7)
Kraken_aRNA = subset(Kraken, Kraken$Type=="aRNA")
Kraken_aRNA$alpha = c(0.7, 0.2, 0.2, 0.7, 0.2, 0.2, 0.2, 0.2, 0.2, 0.7)




setwd("") #Insert path to output Figure PDFs

#Figure 1B
pdf("aRNA_Mammoth_Figure_1B.pdf", height = 40/25.4, width = 80/25.4)
ggplot(Kraken, aes(x = ID, y = log10(Kmers+1), fill=Type, color=Type)) + 
  geom_bar(stat="identity", width=0.7, position=position_dodge(0.84), alpha=0.7) +
  
  theme_classic() + 
  
  xlab("") + ylab(expression(paste("log" ["10"], "(Afrotherian ", italic("k-"), "mers)"))) +
  scale_fill_manual(values=c("aRNA" = "#5D987B", "aDNA" = "#4D4DA1")) +
  scale_color_manual(values=c("#4D4DA1", "#5D987B")) +
  
  scale_y_continuous(expand=c(0,0), limits = c(0, 6.5), breaks = seq(0, 6.5, 2),
                     labels = c("0", "2", "4", "6")) +
  scale_x_discrete(labels = c("Mammoth 1", "Mammoth 2", "Mammoth 3", "Mammoth 4", "Mammoth 5", 
                              "Mammoth 6", "Mammoth 7", "Mammoth 8", "Mammoth 9", "Mammoth 10")) +
  
  theme(legend.key.size = unit(0.65, 'lines'),
        legend.spacing.x = unit(0.2,"lines"),
        legend.margin = margin(l=-5),
        legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=-0.2, l=0.25, r=0.25, "cm")) +
  
  labs(title= "KrakenUniq sequence assignment") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



#Figure 1C
pdf("aRNA_Mammoth_Figure_1C.pdf", height = 40/25.4, width = 40/25.4)
ggplot(Kraken_aDNA, aes(x=log10(TaxReads+1), y=log10(Kmers+1))) +
  geom_point(size=2, color="#4D4DA1", alpha=Kraken_aDNA$alpha) +
  
  theme_classic() +
  
  xlab(expression(paste("log" ["10"], "(Taxonomic sequences)"))) + ylab(expression(paste("log" ["10"], "(", italic("k-"), "mers)"))) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1),
                     labels = c("0","1", "2", "3", "4", "5")) +
  scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2),
                     labels = c("0", "2", "4", "6")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("aDNA coverage") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold")) +
  
  geom_text(label="1", x=3.7, y=5.75, color = "black", size=2) +
  geom_text(label="4", x=3.2, y=4.35, color = "black", size=2) +
  geom_text(label="10", x=4.7, y=6.5, color = "black", size=2)

dev.off()



#Figure 1D
pdf("aRNA_Mammoth_Figure_1D.pdf", height = 40/25.4, width = 40/25.4)
ggplot(Kraken_aRNA, aes(x=log10(TaxReads+1), y=log10(Kmers+1))) +
  geom_point(size=2, alpha=Kraken_aRNA$alpha, color="#5D987B") +
  
  theme_classic() +
  
  xlab(expression(paste("log" ["10"], "(Taxonomic sequences)"))) + ylab(expression(paste("log" ["10"], "(", italic("k-"), "mers)"))) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1),
                     labels = c("0","1", "2", "3", "4", "5")) +
  scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2),
                     labels = c("0", "2", "4", "6")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("aRNA coverage") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold")) +
  
  geom_text(label="1", x=4.2, y=4.5, color = "black", size=2) +
  geom_text(label="4", x=4.22, y=3.3, color = "black", size=2) +
  geom_text(label="10", x=1.5, y=2.7, color = "black", size=2)

dev.off()



