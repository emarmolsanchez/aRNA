#######################################
#                                     #
#     Figures S2A-B. Mammoth aRNA     #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggpubr)

setwd("") #Insert path to Figure auxiliary files



Align = read.delim("Figure_S2A_Table.txt", header=T)
Align$Reference = factor(Align$Reference, levels=c("Mammoth 10 3D-guided", "mEleMax1"))

Length = read.delim("Figure_S2B_Table.txt", header=T)
Length$Reference = factor(Length$Reference, levels=c("Mammoth 10 3D-guided", "mEleMax1"))




setwd("") #Insert path to output Figure PDFs

#Figure S2A-B
Align1 = ggplot(Align, aes(x = Length, y = Mismatch, color=Reference, fill=Reference)) +
  geom_point(size = 0.2, alpha=0.7) + geom_line(linewidth=0.2) +
  
  theme_classic() +
  
  ylab("Mismatch frequency (%)") +
  xlab("Read length (nt)") +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02),
                     labels = c("0", "2", "4", "6", "8")) +
  scale_x_continuous(limits = c(15, 65), breaks = seq(15, 65, 10),
                     labels = c("15", "25", "35", "45", "55", "65")) +
  
  scale_fill_manual(values=c("mEleMax1" = "grey60", "Mammoth 10 3D-guided" = "darkorange2"), 
                    labels=c("mEleMax1" = "mEleMax1", "Mammoth 10 3D-guided" = "Mammoth 10 3D-guided")) +
  scale_color_manual(values=c("mEleMax1" = "grey60", "Mammoth 10 3D-guided" = "darkorange2"), 
                     labels=c("mEleMax1" = "mEleMax1", "Mammoth 10 3D-guided" = "Mammoth 10 3D-guided")) +
  guides(color = guide_legend(override.aes = list(size = 3, linewidth = 0)),
         fill = guide_legend(override.aes = list(alpha = 0.7))) +
  
  theme(legend.key.size = unit(0.65, 'lines'),
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(1, "lines")) + 
  theme(legend.spacing.x = unit(1,"lines"),
        legend.margin = margin(t=-5)) +
  
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(r = 5))) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Genome Reference comparison") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))


Align2 = ggplot(Length, aes(x = Length, y = Count, color=Reference, fill=Reference)) +
  geom_ribbon(aes(ymin = 0, ymax = Count), alpha = 0.4) +
  geom_line(linewidth=0.5) + 
  
  theme_classic() +
  
  ylab("Read count") +
  xlab("Read length (nt)") +
  scale_y_continuous(limits = c(0, 150000), breaks = seq(0, 150000, 25000),
                     labels = c("0", "25K", "50K", "75K", "100K", "125K", "150K")) +
  scale_x_continuous(limits = c(15, 65), breaks = seq(15, 65, 10),
                     labels = c("15", "25", "35", "45", "55", "65")) +
  
  scale_fill_manual(values=c("mEleMax1" = "grey60", "Mammoth 10 3D-guided" = "darkorange2"), 
                    labels=c("mEleMax1" = "mEleMax1", "Mammoth 10 3D-guided" = "Mammoth 10 3D-guided")) +
  scale_color_manual(values=c("mEleMax1" = "grey60", "Mammoth 10 3D-guided" = "darkorange2"), 
                     labels=c("mEleMax1" = "mEleMax1", "Mammoth 10 3D-guided" = "Mammoth 10 3D-guided")) +
  guides(color = guide_legend(override.aes = list(size = 3, linewidth = 0)),
         fill = guide_legend(override.aes = list(alpha = 0.7))) +
  
  theme(legend.key.size = unit(0.65, 'lines'),
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(1, "lines")) + 
  theme(legend.spacing.x = unit(1,"lines"),
        legend.margin = margin(t=-5)) +
  
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(r = 5))) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, "cm")) +
  
  ggtitle("Mapped read length distribution") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))



pdf("aRNA_Mammoth_Figure_S2A-B.pdf", height = 45/25.4, width = 165/25.4)
ggarrange(Align1, Align2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

dev.off()
