#######################################
#                                     #
#      Figure S5A. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)
library(umap)
library(ggh4x)
library(ggstar)

setwd("") #Insert path to Figure auxiliary files



Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2 = read.table("Figure_6A_Table.txt")

#Define tissue and color identities
tissue = c(rep("Skeletal muscle", 6), rep("Heart", 9), rep("Skin", 6), rep("Liver", 5), rep("Lung", 8), rep("Kidney", 4),
           rep("Brain", 3), rep("Testis", 8), rep("Ovary", 4), rep("Placenta", 7), rep("Cartilage", 3), rep("Adipose tissue", 5), rep("Bone marrow", 8),
           rep("Lymph node", 13), rep("Thyroid gland", 9),  rep("Salivary gland", 6))

color = c(rep("#DB6443", 6), rep("plum", 9), rep("#F0C57F",6),rep("goldenrod3", 5), rep("goldenrod4", 8), rep("lightblue", 4), rep("lavender", 3), rep("lavenderblush2", 8), rep("lavenderblush3", 4),
          rep("rosybrown2", 7), rep("lightskyblue3", 3), rep("lightgoldenrod1", 5), rep("lightgoldenrod3", 8), rep("darkseagreen2", 13), rep("darkolivegreen3", 9), rep("darkseagreen4", 6))


#Remove data for running UMAP on Human expression tissue atlas only
Atlas_Human_filtered_tmm = Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2[,1:104]
Dog_tmm = Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2[,105:109]
Thy_tmm = Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2[,110:111]
Elephant_tmm = Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2[,112:119]
Mammoth_tmm = Atlas_Human_filtered_Dog_Thy_Elephant_Mammoth_tmm_log2[,120:123]

#Run UMAP algorithm and predict projection of Elephant, Mammoth, Thylacine and Canid tissues
umap = umap(t(Atlas_Human_filtered_tmm), n_neighbors=10, metric="pearson", spread=10, random_state=40)
UMAP_table = as.data.frame(cbind(umap$layout, tissue, color))
Elephant_predict = predict(umap, t(Elephant_tmm))
Mammoth_predict = predict(umap, t(Mammoth_tmm))
Thy_predict = predict(umap, t(Thy_tmm))
Dog_predict = predict(umap, t(Dog_tmm))

#Arrange tables for plotting
colnames(UMAP_table) = c("UMAP1", "UMAP2", "Tissue", "Color")
UMAP_table$Shape = factor(rep(15, nrow(UMAP_table)))
UMAP_table$UMAP1 = as.numeric(UMAP_table$UMAP1)
UMAP_table$UMAP2 = as.numeric(UMAP_table$UMAP2)

Elephant_predict = as.data.frame(Elephant_predict)
colnames(Elephant_predict) = c("UMAP1", "UMAP2")
Elephant_predict$Tissue = c("Heart", "Heart", "Lung", "Lymph node", "Ovary", "Salivary gland", "Thyroid gland", "Placenta")
Elephant_predict$Color = c("plum", "plum", "goldenrod4", "darkseagreen2", "lavenderblush3", "darkseagreen4", "darkolivegreen3", "rosybrown2")
Elephant_predict$Shape = factor(rep(28, nrow(Elephant_predict)))
Elephant_predict$UMAP1 = as.numeric(Elephant_predict$UMAP1)
Elephant_predict$UMAP2 = as.numeric(Elephant_predict$UMAP2)

Mammoth_predict = as.data.frame(Mammoth_predict)
colnames(Mammoth_predict) = c("UMAP1", "UMAP2")
Mammoth_predict$Tissue = c(rep("Skeletal muscle",3), "Skin")
Mammoth_predict$Color = c(rep("#DB6443",3), "#F0C57F")
Mammoth_predict$Shape = factor(rep(1, nrow(Mammoth_predict)))
Mammoth_predict$UMAP1 = as.numeric(Mammoth_predict$UMAP1)
Mammoth_predict$UMAP2 = as.numeric(Mammoth_predict$UMAP2)

Thy_predict = as.data.frame(Thy_predict)
colnames(Thy_predict) = c("UMAP1", "UMAP2")
Thy_predict$Tissue = c("Skeletal muscle", "Skin")
Thy_predict$Color = c("#DB6443", "#F0C57F")
Thy_predict$Shape = factor(rep(11, nrow(Thy_predict)))
Thy_predict$UMAP1 = as.numeric(Thy_predict$UMAP1)
Thy_predict$UMAP2 = as.numeric(Thy_predict$UMAP2)

Dog_predict = as.data.frame(Dog_predict)
colnames(Dog_predict) = c("UMAP1", "UMAP2")
Dog_predict$Tissue = c("Cartilage", "Liver", "Skeletal muscle", "Skin", "Skin")
Dog_predict$Color = c("lightskyblue3","goldenrod3", "#DB6443", "#F0C57F", "#F0C57F")
Dog_predict$Shape = factor(c(rep(23,3), rep(30, 2)))
Dog_predict$UMAP1 = as.numeric(Dog_predict$UMAP1)
Dog_predict$UMAP2 = as.numeric(Dog_predict$UMAP2)

UMAP_table = rbind(UMAP_table, Thy_predict, Dog_predict)
UMAP_table$Tissue = factor(UMAP_table$Tissue, levels=c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                                                       "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"))
UMAP_table$Color = factor(UMAP_table$Color, levels=c("#DB6443","plum","#F0C57F","goldenrod3", "goldenrod4", "lightblue","lavender","lavenderblush2","lavenderblush3",
                                                     "rosybrown2","lightskyblue3","lightgoldenrod1","lightgoldenrod3","darkseagreen2","darkolivegreen3","darkseagreen4"))
UMAP_table$Shape = factor(UMAP_table$Shape, levels=c(15, 11, 23, 30))
UMAP_table$Size = c(rep(2, 104), rep(4, 7))
UMAP_table$Alpha = c(rep(0.5, 104), rep(1, 7))



axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm"))


color_final = c(rep("#DB6443", 6), rep("plum", 9), rep("#F0C57F",6),rep("goldenrod3", 5), rep("goldenrod4", 8), rep("lightblue", 4), rep("lavender", 3), rep("lavenderblush2", 8), rep("lavenderblush3", 4),
                rep("rosybrown2", 7), rep("lightskyblue3", 3), rep("lightgoldenrod1", 5), rep("lightgoldenrod3", 8), rep("darkseagreen2", 13), rep("darkolivegreen3", 9), rep("darkseagreen4", 6),
                rep("black", 7))

legend_border = c("#DB6443","plum","#F0C57F","goldenrod3", "goldenrod4", "lightblue","lavender","lavenderblush2","lavenderblush3",
                  "rosybrown2","lightskyblue3","lightgoldenrod1","lightgoldenrod3","darkseagreen2","darkolivegreen3","darkseagreen4")





setwd("") #Insert path to output Figure PDFs

#Figure S5A
pdf("aRNA_Mammoth_Figure_S7A.pdf", height = 70/25.4, width = 82.5/25.4)
ggplot(UMAP_table, aes(x=UMAP1, y=UMAP2, fill=Color, starshape=Shape)) + 
  geom_star(size=UMAP_table$Size, color= color_final, alpha=UMAP_table$Alpha) + 
  
  scale_fill_manual(labels = c("Skeletal muscle", "Heart", "Skin", "Liver", "Lung", "Kidney", "Brain", "Testis", "Ovary",
                               "Placenta", "Cartilage", "Adipose tissue", "Bone marrow", "Lymph node", "Thyroid gland", "Salivary gland"),
                    values = c("#DB6443","plum","#F0C57F","goldenrod3", "goldenrod4", "lightblue","lavender","lavenderblush2","lavenderblush3",
                               "rosybrown2","lightskyblue3","lightgoldenrod1","lightgoldenrod3","darkseagreen2","darkolivegreen3","darkseagreen4")) +
  scale_starshape_manual(labels = c("Human", "Thylacine", "Tumat canid", "Wolf"),
                         values = c(15, 11, 23, 30),
                         guide = guide_legend(override.aes =list(fill = "black", size=2.5))) +
  guides(fill = guide_legend(override.aes =list(starshape=15, colour=legend_border, size=2.5), nrow=16)) +
  
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
  
  ggtitle("Protein-coding mRNAs") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold"))

dev.off()
