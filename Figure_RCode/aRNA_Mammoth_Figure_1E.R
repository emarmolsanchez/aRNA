#######################################
#                                     #
#      Figure 1E. Mammoth aRNA        #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)

setwd("") #Insert path to Figure auxiliary files



Age_Length = read.table("Figure1E_Table_Age_vs_Fragment_Length.txt", header=T)





setwd("") #Insert path to output Figure PDFs

#Figure 1E
pdf("aRNA_Mammoth_Figure_1E.pdf", height = 40/25.4, width = 40/25.4)
ggplot(Age_Length, aes(x=Age, y=Length, group=Type, fill=Type, color=Type)) +
  geom_point(size=2, alpha=0.7) + geom_line(linewidth=0.5) +
  
  theme_classic() +
  
  ylab("Average length (nt)") + xlab("Estimated age") +
  scale_y_continuous(limits = c(15, 55), breaks = seq(15, 55, 10),
                     labels = c("15","25", "35", "45", "55")) +
  scale_x_continuous(limits = c(35000, 50000), breaks = seq(35000, 50000, 5000),
                     labels = c("35K", "40K", "45K", "50K")) +
  
  scale_fill_manual(values=c("aRNA" = "#5D987B", "aDNA" = "#4D4DA1")) +
  scale_color_manual(values=c("#4D4DA1", "#5D987B")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = 0, l = 0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, "lines"),
        legend.key.height = unit(1.5, "lines"),
        legend.spacing.x = unit(0.4,"lines"),
        legend.margin = margin(t=-10),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7, margin = margin(l = 2))) +
  
  ggtitle("Fragment length") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
