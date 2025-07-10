#######################################
#                                     #
#      Figure S2C. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(gtools)

setwd("") #Insert path to Figure auxiliary files




Table_final = read.delim("Figure_S2C_Table.txt", header=T)

Table_final$FC = foldchange(Table_final$Mammoth.1.deep, Table_final$Mammoth.1)
Table_final$log2FC = log2(foldchange(Table_final$Mammoth.1.deep, Table_final$Mammoth.1))
Table_final$Type = factor(Table_final$Type, levels=c("Raw sequences", "Trimmed sequences", "PCR-deduplicated", "Mapped & UMI-deduplicated",
                                                     "Mapped to Mitochondrial genome", "Mapped to Nuclear genome",
                                                     "Mapped to gene loci", "Mapped to exonic loci",
                                                     "Mapped to intronic loci", "Mapped to intergenic regions",
                                                     "Mapped to Transcriptome & UMI-deduplicated", "Mapped to PC transcripts",
                                                     "Mapped to NC transcripts", "Number of PC loci with Coverage ≥5%",
                                                     "Number of NC loci with Coverage ≥5%"))





setwd("") #Insert path to output Figure PDFs

#Figure S2C
pdf("aRNA_Mammoth_Figure_S2C.pdf", height = 75/25.4, width = 82.5/25.4)
ggplot(Table_final, aes(y = FC, x= Type, color=Type)) + 
  geom_bar(stat="identity", color= "#5D987B", fill="#5D987B", alpha=0.7, width = 0.5) +
  
  theme_classic() + 
  
  ylab("Fold change") + xlab("") +
  scale_y_continuous(expand=c(0,0), limits = c(0, 10), breaks = seq(0, 10, 2),
                     labels = c("0", "2", "4", "6", "8", "10")) +
  
  theme(axis.text.x = element_text(size=6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, l=0.25, b=-0.3, r=0.25, "cm")) +
  
  ggtitle("Deep sequencing data gain") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
