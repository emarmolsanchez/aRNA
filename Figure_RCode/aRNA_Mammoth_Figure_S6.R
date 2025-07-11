#######################################
#                                     #
#      Figures S6. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(ggh4x)
library(ggstar)

setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
TTN_Expr = read.table("Figure_S6_Table3.txt", row.names=1, header=T)


#Keep Humans and Mammoth 1
Expr_1 = TTN_Expr[,1:104]
Expr_2 = TTN_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                       "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$TTN_Expr = as.numeric(Expr)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_TTN.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=TTN_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$TTN_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,20),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("TTN") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()



setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
XIRP2_Expr = read.table("Figure_S6_Table4.txt", row.names=1, header=T)

#Keep Humans and Mammoth 1
Expr_1 = XIRP2_Expr[,1:104]
Expr_2 = XIRP2_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                             "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$XIRP2_Expr = as.numeric(Expr)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_XIRP2.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=XIRP2_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$XIRP2_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,20),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("XIRP2") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()






setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
ACTA1_Expr = read.table("Figure_S6_Table5.txt", row.names=1, header=T)

#Keep Humans and Mammoth 1
Expr_1 = ACTA1_Expr[,1:104]
Expr_2 = ACTA1_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                             "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$ACTA1_Expr = as.numeric(Expr)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_ACTA1.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=ACTA1_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$ACTA1_Expr),
                       breaks = c(0, 5, 10, 15, 20, 25),
                       labels = c("0", "5", "10", "15", "20", "25"),
                       limits = c(0,25),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("ACTA1") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()






setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
NEB_Expr = read.table("Figure_S6_Table6.txt", row.names=1, header=T)

#Keep Humans and Mammoth 1
Expr_1 = NEB_Expr[,1:104]
Expr_2 = NEB_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                             "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$NEB_Expr = as.numeric(Expr)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_NEB.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=NEB_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$NEB_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,25),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("NEB") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()





setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
TNNC1_Expr = read.table("Figure_S6_Table7.txt", row.names=1, header=T)

#Keep Humans and Mammoth 1
Expr_1 = TNNC1_Expr[,1:104]
Expr_2 = TNNC1_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                             "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$TNNC1_Expr = as.numeric(Expr)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_TNNC1.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=TNNC1_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$TNNC1_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,25),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("TNNC1") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()





setwd("") #Insert path to Figure auxiliary files



UMAP_table_PC = read.delim("Figure_S6_Table1.txt", row.names=1)
TNNT1_Expr = read.table("Figure_S6_Table8.txt", row.names=1, header=T)

#Keep Humans and Mammoth 1
Expr_1 = TNNT1_Expr[,1:104]
Expr_2 = TNNT1_Expr[,120]
Expr = as.vector(cbind(Expr_1, Expr_2))

UMAP_table1 = UMAP_table_PC[1:104,]
UMAP_table2 = UMAP_table_PC[113,]
UMAP_table_PC = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_PC$Tissue = factor(UMAP_table_PC$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                             "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table_PC$Shape = factor(UMAP_table_PC$Shape, levels=c(15, 1))
UMAP_table_PC$Size = c(rep(2, 104), rep(4, 1))
UMAP_table_PC$Alpha = c(rep(1, 105))
UMAP_table_PC$TNNT1_Expr = as.numeric(Expr)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))




setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_TNNT1.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_PC, aes(x=UMAP1, y=UMAP2, fill=TNNT1_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_PC$Size, alpha=UMAP_table_PC$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_PC$TNNT1_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,25),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-35,25)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("TNNT1") * " mRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()





setwd("") #Insert path to Figure auxiliary files



UMAP_table_small = read.delim("Figure_S6_Table2.txt", row.names=1)
MIR1_Expr = read.table("Figure_S6_Table9.txt", row.names=1, header=T)


#Keep Humans and Mammoth 1
MIR1_1 = MIR1_Expr[,1:71]
MIR1_2 = MIR1_Expr[,73]
MIR1_Expr = as.vector(cbind(MIR1_1, MIR1_2))

UMAP_table1 = UMAP_table_small[1:71,]
UMAP_table2 = UMAP_table_small[73,]
UMAP_table_small = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_small$Tissue = factor(UMAP_table_small$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Brain", "Ovary",
                                                       "Placenta", "Cartilage", "Pancreas", "Plasma", "Stomach", "Jejunum", "Ileum",
                                                       "Duodenum", "Thyroid gland"))


UMAP_table_small$Shape = factor(UMAP_table_small$Shape, levels=c(15, 1))
UMAP_table_small$Size = c(rep(2, 71), rep(4, 1))
UMAP_table_small$Alpha = c(rep(1, 72))
UMAP_table_small$MIR1_Expr = as.numeric(MIR1_Expr)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))





setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_MIR1.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_small, aes(x=UMAP1, y=-UMAP2, fill=MIR1_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_small$Size, alpha=UMAP_table_small$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_small$MIR1_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,20),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-25,35)) + xlim(c(-25,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("Mir-1-P3") * " microRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()






setwd("") #Insert path to Figure auxiliary files



UMAP_table_small = read.delim("Figure_S6_Table2.txt", row.names=1)
MIR122_Expr = read.table("Figure_S6_Table10.txt", row.names=1, header=T)


#Keep Humans and Mammoth 1
MIR122_1 = MIR122_Expr[,1:71]
MIR122_2 = MIR122_Expr[,73]
MIR122_Expr = as.vector(cbind(MIR122_1, MIR122_2))

UMAP_table1 = UMAP_table_small[1:71,]
UMAP_table2 = UMAP_table_small[73,]
UMAP_table_small = as.data.frame(rbind(UMAP_table1, UMAP_table2))


UMAP_table_small$Tissue = factor(UMAP_table_small$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Brain", "Ovary",
                                                                   "Placenta", "Cartilage", "Pancreas", "Plasma", "Stomach", "Jejunum", "Ileum",
                                                                   "Duodenum", "Thyroid gland"))


UMAP_table_small$Shape = factor(UMAP_table_small$Shape, levels=c(15, 1))
UMAP_table_small$Size = c(rep(2, 71), rep(4, 1))
UMAP_table_small$Alpha = c(rep(1, 72))
UMAP_table_small$MIR122_Expr = as.numeric(MIR122_Expr)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))





setwd("") #Insert path to output Figure PDFs

#Figures S6
pdf("aRNA_Mammoth_Figure_S6_MIR122.pdf", height = 50/25.4, width = 82.5/25.4)
ggplot(UMAP_table_small, aes(x=UMAP1, y=-UMAP2, fill=MIR122_Expr, starshape=Shape)) + 
  geom_star(size=UMAP_table_small$Size, alpha=UMAP_table_small$Alpha) + 
  
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "red3", midpoint = median(UMAP_table_small$MIR122_Expr),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "5", "10", "15", "20"),
                       limits = c(0,20),
                       
                       guide = guide_colorbar(ticks=TRUE, ticks.colour = "grey30", draw.ulim = FALSE, draw.llim = FALSE,
                                              barwidth = 0.5, barheight = 4, frame.colour = "grey30",
                                              frame.linewidth = 0.2, order=1)) +
  
  scale_starshape_manual(labels = c("Human", "Mammoth"),
                         values = c(15, 1),
                         guide = guide_legend(override.aes = list(shape=15, fill = "black", size=2.5))) +
  
  theme_minimal() +
  
  ylim(c(-25,35)) + xlim(c(-25,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(2,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="center") +
  
  ggtitle(bquote(italic("Mir-122") * " microRNA abundance")) + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
