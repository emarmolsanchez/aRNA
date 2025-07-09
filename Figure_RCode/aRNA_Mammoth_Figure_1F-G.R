#######################################
#                                     #
#      Figures 1F-G. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggpubr)

setwd("") #Insert path to Figure auxiliary files



Align = read.table("Figure1F_Table_Alignment_Mismatch_Comparison.txt", header=T)
Align2 = read.table("Figure1G_Table_Alignment_Read_Length_Distribution.txt", header=T)





setwd("") #Insert path to output Figure PDFs

#Figure 1F-G
Align_plot = ggplot(Align, aes(x = Length, y = Mismatch, color=Algorithm)) +
  geom_point(size = 0.2, alpha=0.7) + geom_line(linewidth=0.2) +
  
  theme_classic() +
  
  ylab("Mismatch frequency (%)") +
  xlab("Sequence length (nt)") +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.02),
                     labels = c("0", "2", "4", "6", "8", "10", "12")) +
  scale_x_continuous(limits = c(15, 65), breaks = seq(15, 65, 10),
                     labels = c("15", "25", "35", "45", "55", "65")) +
  
  scale_fill_manual(values=c("Bowtie" = "dodgerblue", "Bowtie2" = "red3", "BWA" = "grey60"), labels=c("Bowtie", "Bowtie2", "BWA aln")) +
  scale_color_manual(values=c("dodgerblue", "red3", "grey60")) +
  guides(color = guide_legend(override.aes = list(size = 3, linewidth = 0)),
         fill = guide_legend(override.aes = list(alpha = 0.7))) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.key.height = unit(1.5, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Algorithm comparison") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

Align_plot2 = ggplot(Align2, aes(x = Length, y = Count, color=Algorithm, fill=Algorithm)) +
  geom_ribbon(aes(ymin = 0, ymax = Count), alpha = 0.7) +
  geom_line(linewidth=0.5) + 
  
  theme_classic() +
  
  ylab("Sequence count") +
  xlab("Sequence length (nt)") +
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 20000),
                     labels = c("0", "20K", "40K", "60K", "80K", "100K")) +
  scale_x_continuous(limits = c(15, 65), breaks = seq(15, 65, 10),
                     labels = c("15", "25", "35", "45", "55", "65")) +
  
  scale_fill_manual(values=c("Bowtie" = "dodgerblue", "Bowtie2" = "red3", "BWA" = "grey60"), labels=c("Bowtie", "Bowtie2", "BWA aln")) +
  scale_color_manual(values=c("dodgerblue", "red3", "grey60")) +
  guides(color = guide_legend(override.aes = list(shape = NA, linetype=1, fill=NA)), fill = "none") +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.key.height = unit(1.5, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-5),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mapped aRNAs length distribution") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))


pdf("aRNA_Mammoth_Figure_1F-G.pdf", height = 40/25.4, width = 125/25.4)
ggarrange(Align_plot, Align_plot2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

dev.off()
