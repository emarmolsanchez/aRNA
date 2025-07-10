#######################################
#                                     #
#       Figure S6B. Mammoth aRNA      #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



Atlas_melted = read.table("Figure_S6B_Table.txt", header=T)
Atlas_melted$variable = factor(Atlas_melted$variable, levels=c("Human_Stomach",
                                                               "Human_Illeum",
                                                               "Human_Jejunum",
                                                               "Human_Duodenum",
                                                               "Human_Liver",
                                                               "Human_Pancreas",
                                                               "Human_Ovary",
                                                               "Human_Placenta",
                                                               "Human_Plasma",
                                                               "Human_Skin",
                                                               "Human_Lung",
                                                               "Human_Cartilage",
                                                               "Human_Brain", 
                                                               "Human_Thyroid",
                                                               "Human_Heart", 
                                                               "Human_Muscle", 
                                                               "Yuka_muscle"))
genes = Atlas_melted$genes[1:30]
Atlas_melted$genes = factor(Atlas_melted$genes, levels=genes)




setwd("") #Insert path to output Figure PDFs

#Figure S6B
pdf("aRNA_Mammoth_Figure_S6B.pdf", height = 60/25.4, width = 165/25.4)
ggplot(Atlas_melted, aes(x = genes, y = variable, fill = value)) +
  geom_tile(color = "grey60", lwd = 0.2, linetype = 1) +
  scale_fill_gradient(low = "white", high = "#DB6443") +
  
  xlab("") + ylab("") +
  scale_y_discrete(labels = c("Stomach",
                              "Illeum",
                              "Jejunum",
                              "Duodenum",
                              "Liver",
                              "Pancreas",
                              "Ovary", 
                              "Placenta",
                              "Plasma", 
                              "Skin",
                              "Lung",      
                              "Cartilage",
                              "Brain",
                              "Thyroid gland", 
                              "Heart",         
                              "Skeletal muscle",
                              "Mammoth 1 (muscle)")) +
  
  theme_classic() +
  
  guides(fill = guide_colourbar(ticks=TRUE, ticks.colour = "grey60", draw.ulim = FALSE, draw.llim = FALSE,
                                barwidth = 0.5, barheight = 5, frame.colour = "grey60",
                                frame.linewidth = 0.2)) +
  
  theme(legend.position = "right", 
        legend.text = element_text(size=6),
        legend.title = element_blank()) +
  
  ggtitle("Mammoth 1 microRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold")) +
  
  theme(axis.text.x = element_text(size=6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=6),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t=0.25, b=0, l=0, r=0.25, "cm")) +
  
  coord_cartesian(xlim=c(1,30), clip = "off") +
  annotate("text", x = 31.3, y = 8.6, label = expression(paste("log" ["2"] ,"(Normalized Read count)")),
           size = 2, angle=90)

dev.off()
