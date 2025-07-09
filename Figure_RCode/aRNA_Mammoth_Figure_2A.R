#######################################
#                                     #
#      Figure 2A. Mammoth aRNA        #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



USP9Y = read.table("Figure2A_Sample1_Trimmed_Collapsed_UMIheader_Bowtie2_AsianElephant_Transcriptome_PC-NC-miRNAs_UMIDEDUP_CLEAN_proteincoding_PerBase_Depth_Coverage_LOC126069875_USP9Y.txt")
USP9Y$ID = rep("Mammoth1", nrow(USP9Y))





setwd("") #Insert path to output Figure PDFs

#Figure 2A
pdf("aRNA_Mammoth_Figure_2A.pdf", height = 25/25.4, width = 82.5/25.4)
ggplot(USP9Y, aes(x=V2, y=V3, color=ID)) +
  geom_line(linewidth=0.4) +
  
  theme_classic() +
  
  xlab("") + ylab("Sequence count") +
  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, 1),
                     labels = c("0", "1", "2", "3", "4")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(USP9Y$V2), max(USP9Y$V2)), labels = c("-5'", "-3'")) +
  scale_color_manual(labels=c("Mammoth 1"),
                     values=c("#5D987B")) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7)) +
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(0.2,"lines"),
        legend.margin = margin(l=-5)) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = 0.25, "cm")) +
  
  ggtitle("LOC126069875 (USP9Y)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()