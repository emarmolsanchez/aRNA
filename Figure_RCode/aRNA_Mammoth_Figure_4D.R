#######################################
#                                     #
#       Figure 4D. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(umap)
library(ggh4x)
library(ggstar)

setwd("") #Insert path to Figure auxiliary files



miRNAs_Human_Elephant_Mammoth_Dog_Thy_tmm_log2 = read.table("Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_log2TMM_Normalized_Expression_microRNAs.txt")

#Define tissue and color identities
tissue = c(rep("Cartilage", 5), rep("Skeletal muscle", 4), rep("Heart", 4), rep("Skin", 8), rep("Plasma", 6), rep("Lung", 4), rep("Brain", 4), rep("Liver", 4), rep("Ovary", 4),
           rep("Thyroid gland", 4), rep("Placenta", 4), rep("Pancreas", 4), rep("Stomach", 4), rep("Jejunum", 4), rep("Ileum", 4), rep("Duodenum", 4))

color = c(rep("lightskyblue3", 5), rep("#DB6443", 4), rep("plum", 4), rep("#F0C57F", 8), rep("lemonchiffon2", 6), rep("goldenrod4", 4),
          rep("lavender", 4), rep("goldenrod3", 4), rep("lavenderblush3", 4), rep("darkolivegreen3", 4), rep("rosybrown2", 4),
          rep("ivory4", 4), rep("palegreen1", 4), rep("palegreen3", 4), rep("olivedrab2", 4), rep("olivedrab4", 4))


#Remove data for running UMAP on Human expression tissue atlas only
Human_tmm = miRNAs_Human_Elephant_Mammoth_Dog_Thy_tmm_log2[,1:71]
Elephant_Mammoth_tmm = miRNAs_Human_Elephant_Mammoth_Dog_Thy_tmm_log2[,72:74]
Dog_tmm = miRNAs_Human_Elephant_Mammoth_Dog_Thy_tmm_log2[,75:79]
Thy_tmm =  miRNAs_Human_Elephant_Mammoth_Dog_Thy_tmm_log2[,80:81]

#Run UMAP algorithm and predict projection of Elephant, Mammoth, Thylacine and Canid tissues
umap = umap(t(Human_tmm), n_neighbors=10, metric="pearson", spread=10, random_state=40)
UMAP_table = as.data.frame(cbind(umap$layout, tissue, color))
Elephant_Mammoth_predict = predict(umap, t(Elephant_Mammoth_tmm))
Dog_predict = predict(umap, t(Dog_tmm))
Thy_predict = predict(umap, t(Thy_tmm))

#Arrange tables for plotting
colnames(UMAP_table) = c("UMAP1", "UMAP2", "Tissue", "Color")
UMAP_table$Shape = factor(rep(15, nrow(UMAP_table)))
UMAP_table$UMAP1 = as.numeric(UMAP_table$UMAP1)
UMAP_table$UMAP2 = as.numeric(UMAP_table$UMAP2)

Elephant_Mammoth_predict = as.data.frame(Elephant_Mammoth_predict)
colnames(Elephant_Mammoth_predict) = c("UMAP1", "UMAP2")
Elephant_Mammoth_predict$Tissue = c("Plasma", rep("Skeletal muscle",2))
Elephant_Mammoth_predict$Color = c("lemonchiffon2", rep("#DB6443",2))
Elephant_Mammoth_predict$Shape = factor(c(28, rep(1,2)))
Elephant_Mammoth_predict$UMAP1 = as.numeric(Elephant_Mammoth_predict$UMAP1)
Elephant_Mammoth_predict$UMAP2 = as.numeric(Elephant_Mammoth_predict$UMAP2)

Thy_predict = as.data.frame(Thy_predict)
colnames(Thy_predict) = c("UMAP1", "UMAP2")
Thy_predict$Tissue = c("Skeletal muscle", "Skin")
Thy_predict$Color = c("#DB6443", "#F0C57F")
Thy_predict$Shape = factor(rep(28, nrow(Thy_predict)))
Thy_predict$UMAP1 = as.numeric(Thy_predict$UMAP1)
Thy_predict$UMAP2 = as.numeric(Thy_predict$UMAP2)

Dog_predict = as.data.frame(Dog_predict)
colnames(Dog_predict) = c("UMAP1", "UMAP2")
Dog_predict$Tissue = c("Cartilage", "Liver", "Skeletal muscle", "Skin", "Skin")
Dog_predict$Color = c("lightskyblue3","goldenrod3", "#DB6443", "#F0C57F", "#F0C57F")
Dog_predict$Shape = factor(c(rep(2,3), rep(28, 2)))
Dog_predict$UMAP1 = as.numeric(Dog_predict$UMAP1)
Dog_predict$UMAP2 = as.numeric(Dog_predict$UMAP2)

#Reorder
#Reorder UMAP table
UMAP_table_SKM = subset(UMAP_table, grepl("Skeletal muscle", Tissue))
UMAP_table_Heart = subset(UMAP_table, grepl("Heart", Tissue))
UMAP_table_Skin = subset(UMAP_table, grepl("Skin", Tissue))
UMAP_table_Liver = subset(UMAP_table, grepl("Liver", Tissue))
UMAP_table_Lung = subset(UMAP_table, grepl("Lung", Tissue))
UMAP_table_Brain = subset(UMAP_table, grepl("Brain", Tissue))
UMAP_table_Ovary = subset(UMAP_table, grepl("Ovary", Tissue))
UMAP_table_Placenta = subset(UMAP_table, grepl("Placenta", Tissue))
UMAP_table_Cartilage = subset(UMAP_table, grepl("Cartilage", Tissue))
UMAP_table_Pancreas = subset(UMAP_table, grepl("Pancreas", Tissue))
UMAP_table_Plasma = subset(UMAP_table, grepl("Plasma", Tissue))
UMAP_table_Stomach = subset(UMAP_table, grepl("Stomach", Tissue))
UMAP_table_Jejunum = subset(UMAP_table, grepl("Jejunum", Tissue))
UMAP_table_Ileum = subset(UMAP_table, grepl("Ileum", Tissue))
UMAP_table_Duodenum = subset(UMAP_table, grepl("Duodenum", Tissue))
UMAP_table_Thyroid = subset(UMAP_table, grepl("Thyroid gland", Tissue))

UMAP_table_small = rbind(UMAP_table_SKM, UMAP_table_Heart, UMAP_table_Skin, UMAP_table_Liver, UMAP_table_Lung, UMAP_table_Brain, UMAP_table_Ovary,
                         UMAP_table_Placenta, UMAP_table_Cartilage, UMAP_table_Pancreas, UMAP_table_Plasma, UMAP_table_Stomach, UMAP_table_Jejunum,
                         UMAP_table_Ileum, UMAP_table_Duodenum, UMAP_table_Thyroid)

UMAP_table_small_final = rbind(UMAP_table_small, Elephant_Mammoth_predict, Thy_predict, Dog_predict)

UMAP_table_small = rbind(UMAP_table_small, Elephant_Mammoth_predict)
UMAP_table_small$Size = c(rep(2, 71), rep(2.5, 1), rep(4, 2))
UMAP_table_small$Alpha = c(rep(0.5, 72), rep(1, 2))

#Only select Mammoth 1
UMAP_table_small = UMAP_table_small[1:73,]


UMAP_table_small$Tissue = factor(UMAP_table_small$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Brain", "Ovary",
                                                                   "Placenta", "Cartilage", "Pancreas", "Plasma", "Stomach", "Jejunum", "Ileum",
                                                                   "Duodenum", "Thyroid gland"))
UMAP_table_small$Color = factor(UMAP_table_small$Color, levels=c("#DB6443", "plum", "#F0C57F", "goldenrod3", "goldenrod4", "lavender", 
                                                                 "lavenderblush3", "rosybrown2", "lightskyblue3", "ivory4", "lemonchiffon2",
                                                                 "palegreen1","palegreen3","olivedrab2","olivedrab4","darkolivegreen3"))

UMAP_table_small$Shape = factor(UMAP_table_small$Shape, levels=c(15, 28, 1))

#Define colors for plot
color_final_small = c(rep("#DB6443", 4), rep("plum", 4), rep("#F0C57F",8), rep("goldenrod3",4), rep("goldenrod4",4), rep("lavender",4), 
                      rep("lavenderblush3",4), rep("rosybrown2",4),  rep("lightskyblue3",5), rep("ivory4",4), rep("lemonchiffon2",6),
                      rep("palegreen1",4), rep("palegreen3",4), rep("olivedrab2",4), rep("olivedrab4",4), rep("darkolivegreen3",4), rep("black", 2))



legend_border = c("#DB6443", "plum", "#F0C57F", "goldenrod3", "goldenrod4", "lavender", 
                  "lavenderblush3", "rosybrown2", "lightskyblue3", "ivory4", "lemonchiffon2",
                  "palegreen1","palegreen3","olivedrab2","olivedrab4","darkolivegreen3")
legend_border = factor(legend_border, levels=c("#DB6443", "plum", "#F0C57F", "goldenrod3", "goldenrod4", "lavender", 
                                               "lavenderblush3", "rosybrown2", "lightskyblue3", "ivory4", "lemonchiffon2",
                                               "palegreen1","palegreen3","olivedrab2","olivedrab4","darkolivegreen3"))

#Define axis grid for UMAP arrows
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))



setwd("") #Insert path to output Figure PDFs

#Figure 4D
pdf("aRNA_Mammoth_Figure_4D.pdf", height = 70/25.4, width = 82.5/25.4)
ggplot(UMAP_table_small, aes(x=UMAP1, y=-UMAP2, fill=Color, starshape=Shape)) + 
  geom_star(size=UMAP_table_small$Size, color= color_final_small, alpha=UMAP_table_small$Alpha) + 
  
  scale_fill_manual(labels = c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Brain", "Ovary",
                               "Placenta", "Cartilage", "Pancreas", "Plasma", "Stomach", "Jejunum", "Ileum",
                               "Duodenum", "Thyroid gland"),
                    values = c("#DB6443", "plum", "#F0C57F", "goldenrod3", "goldenrod4", "lavender", 
                               "lavenderblush3", "rosybrown2", "lightskyblue3", "ivory4", "lemonchiffon2",
                               "palegreen1","palegreen3","olivedrab2","olivedrab4","darkolivegreen3")) +
  scale_starshape_manual(labels = c("Human", "Elephant", "Mammoth"),
                         values = c(15, 28, 1),
                         guide = guide_legend(override.aes =list(fill = "black", size=2.5))) +
  guides(fill = guide_legend(order=2, override.aes =list(starshape=15, colour=legend_border, size=2.5), nrow=16)) +
  
  theme_minimal() +
  
  ylim(c(-35,35)) + xlim(c(-35,35)) +
  xlab("UMAP1") + ylab("UMAP2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(linewidth=0.3, arrow = arrow(type = "closed", length = unit(5, 'pt')))) +
  
  theme(axis.title.x =element_text(size=6, margin = margin(t = 5), hjust=0),
        axis.title.y =element_text(size=6, margin = margin(r = 2.5), hjust=0),
        plot.margin = margin(t=0.25, b=0.25, l=0.25, r=0.25, "cm")) +
  
  theme(legend.key.size = unit(0.65, 'lines')) + 
  theme(legend.spacing.y = unit(1.5,"lines"),
        legend.margin = margin(-13,0,0,0)) +
  
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size=7),
        legend.justification ="top") +
  
  ggtitle("microRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
