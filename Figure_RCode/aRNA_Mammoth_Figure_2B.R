#######################################
#                                     #
#      Figure 2B. Mammoth aRNA        #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(scales)

setwd("") #Insert path to Figure auxiliary files



Ema_4729 = read.table("Figure2B_Table1.txt")
L386 = read.table("Figure2B_Table2.txt")
L414 = read.table("Figure2B_Table3.txt")
Yuka = read.table("Figure2B_Table4.txt")

Ema_4729$V4 = Ema_4729$V3/119
L386$V4 = L386$V3/146
L414$V4 = L414$V3/75
Yuka$V4 = Yuka$V3/274

Mean_Ema_4729 = mean(Ema_4729$V4)
SD_Ema_4729 = sd(Ema_4729$V4)
Mean_L386 = mean(L386$V4)
SD_L386 = sd(L386$V4)
Mean_L414 = mean(L414$V4)
SD_L414 = sd(L414$V4)
Mean_Yuka = mean(Yuka$V4)
SD_Yuka = sd(Yuka$V4)

Ema_4729$V5 = rescale((Ema_4729$V4 - Mean_Ema_4729) / SD_Ema_4729)
L386$V5 = rescale((L386$V4 - Mean_L386) / SD_L386)
L414$V5 = rescale((L414$V4 - Mean_L414) / SD_L414)
Yuka$V5 = rescale((Yuka$V4 - Mean_Yuka) / SD_Yuka)

Table_cov = as.data.frame(c(Yuka$V5, L386$V5, L414$V5, Ema_4729$V5))
Table_cov$ID = c(rep("Yuka", 695), rep("L386", 695), rep("L414", 695), rep("Asian elephant", 695))
Table_cov$Length = rep(1:695, 4)
Table_cov$Normalized_counts = c(Yuka$V4,  L386$V4, L414$V4, Ema_4729$V4)
Table_cov$Counts = c(Yuka$V3, L386$V3, L414$V3, Ema_4729$V3)
colnames(Table_cov) = c("Value", "ID", "Length", "Value2", "Counts")
Table_cov$ID = factor(Table_cov$ID, levels=c("Yuka", "L386", "L414", "Asian elephant"))





setwd("") #Insert path to output Figure PDFs

#Figure 2B
pdf("aRNA_Mammoth_Figure_2B.pdf", height = 25/25.4, width = 82.5/25.4)
ggplot(Table_cov, aes(x = Length , y = Value, color=ID)) +
  geom_line(linewidth=ifelse(Table_cov$ID == "Yuka", 0.7, 0.3),
            alpha = ifelse (Table_cov$ID == "Yuka", 0.8, 0.4)) +
  
  theme_classic() +
  
  scale_color_manual(labels=c("Mammoth 1", "Mammoth L386", "Mammoth L414", "Asian elephant"),
                     values=c("#4D4DA1", "goldenrod2", "tomato", "tan")) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7)) +
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(0.2,"lines"),
        legend.margin = margin(l=-5)) +
  
  xlab("") + ylab("Depth score") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(1), max(695)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7, margin = margin(r=2.5)),
        axis.title.x = element_text(size = 7, margin = margin(t=5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = 0.25, "cm")) +
  
  ggtitle("LOC126069867 (SRY)") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()
