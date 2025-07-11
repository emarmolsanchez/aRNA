#######################################
#                                     #
#       Figure 1K. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(tidyverse)

setwd("") #Insert path to Figure auxiliary files



#Load the RNA data
data = read.table("Figure_1K_Table1.txt", header=F)
colnames(data) = c("name", "chr", "start", "end", "counts")
#Annotate the chromosomes
chr_alt = paste0("Chr",seq(1:27))
chr_alt = c(chr_alt, "ChrX", "ChrY")
a = as.vector(table(data$chr))
data$chr_alt = rep(chr_alt, times=a)

data$pos = seq(1,dim(data)[1])*500000 #make a column called position, for the bins in the chromosomes this is a running addition of 250.000 
data$cumsum = cumsum(data$counts)
data_somatic = subset(data, data$chr_alt != "ChrY")
data_somatic = subset(data_somatic, data_somatic$chr_alt != "ChrX")
data_X =  subset(data, data$chr_alt == "ChrX")
data_Y =  subset(data, data$chr_alt == "ChrY")
data_sex = rbind(data_X, data_Y)
data_somatic$alpha = ifelse(data_somatic$counts > 92, "#5D987B", "#D9E8E1")
data_sex$alpha = ifelse(data_sex$counts > 34, "#5D987B", "#D9E8E1")
data = rbind(data_somatic, data_sex)
data$Type = rep("aRNA", nrow(data))

#Load the DNA data
datar = read.table("Figure_1K_Table2.txt", header=F)
colnames(datar) = c("name", "chr", "start", "end", "counts")
#Annotate the chromosomes
chr_alt = paste0("Chr",seq(1:27))
chr_alt = c(chr_alt, "ChrX", "ChrY")
a = as.vector(table(datar$chr))
datar$chr_alt = rep(chr_alt, times=a)

#Arrange the df
datar$pos = seq(1,dim(datar)[1])*500000 #make a column called position, for the bins in the chromosomes this is a running addition of 250.000 
datar$cumsum = cumsum(datar$counts)
datar_somatic = subset(datar, datar$chr_alt != "ChrY")
datar_somatic = subset(datar_somatic, datar_somatic$chr_alt != "ChrX")
datar_X =  subset(datar, datar$chr_alt == "ChrX")
datar_Y =  subset(datar, datar$chr_alt == "ChrY")
datar_sex = rbind(datar_X, datar_Y)
datar_somatic$alpha = ifelse(datar_somatic$counts > 92, "#4D4DA1", "grey85")
datar_sex$alpha = ifelse(datar_sex$counts > 34, "#4D4DA1", "grey85")
datar = rbind(datar_somatic, datar_sex)
datar$Type = rep("aDNA", nrow(datar))

namef = factor(chr_alt, levels = chr_alt, ordered = TRUE)

###########All data
#Find the mean postition of each chromosome for plotting the label
mean_pos = data %>% group_by(chr_alt) %>% summarize(mean_pos = min(pos) + ((max(pos)-min(pos))/2))
mean_pos2 = c(122250000, 362750000, 592250000, 802000000, 984000000, 1139750000, 1281500000, 1418250000, 1549500000, 1673500000,
              1792500000, 1906750000, 2015500000, 2119000000, 2217500000, 2313500000, 2407250000, 2495750000, 2581250000, 2666500000, 
              2750750000, 2833750000, 2914750000, 2988750000, 3057250000, 3116250000, 3170000000, 3270000000, 3368750000)





setwd("") #Insert path to output Figure PDFs

#Figure 1K
pdf("aRNA_Mammoth_Figure_1K.pdf", height = 60/25.4, width = 165/25.4)
ggplot() + 
  geom_point(data = data, aes(x = pos, y = log10(as.numeric(counts)+1), fill=Type), color = data$alpha, size=2, alpha=0.6) +
  geom_point(data = datar, aes(x = pos, y = log10(as.numeric(counts)+1), fill=Type), color = datar$alpha, size=2, alpha=0.6) +
  
  geom_segment(aes(x = 1700000000, y = 4, xend = 2640000000, yend = 4), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_segment(aes(x = 1700000000, y = 5.3, xend = 2640000000, yend = 5.3), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_segment(aes(x = 1700000000, y = 4, xend = 1700000000, yend = 5.3), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_segment(aes(x = 2640000000, y = 4, xend = 2640000000, yend = 5.3), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_text(aes(label="rRNA gene\n clusters", x=2150000000, y=4.6, label="Center"), color = "black", size=2.5, family="Helvetica") +
  
  theme_classic() + 
  
  scale_x_continuous(expand=c(0.015, 0.015), name = element_blank(), labels = namef, breaks = mean_pos2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  ylab(expression(paste("log" ["10"], "(Sequence count)"))) +
  
  theme(legend.key.size = unit(0.65, 'lines'),
        legend.spacing.y = unit(0.4,"lines"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.9, 0.95), 
        legend.title = element_blank(), 
        legend.text = element_text(size=6),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.box.margin = margin(t=0.05, b=0.05, l=0.05, r=0.1, unit = "cm")) +
  
  guides(fill = guide_legend(override.aes=list(shape=19, alpha= 0.6, colour=c("#4D4DA1", "#5D987B"), size=4), nrow=2)) +
  
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.title.y =element_text(size=7, margin = margin(r = 2.5)),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=0, l=0.25, r=0.25, "cm")) +
  
  ggtitle("Mammoth 1 Genome-wide alignment") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))


dev.off()
