#######################################
#                                     #
#     Figures S6C. Mammoth aRNA       #
#   Dr. Emilio Marmol-Sanchez. 2025   #
#                                     #
#######################################


library(ggplot2)


setwd("") #Insert path to Figure auxiliary files



ACTA1 = read.table("Figure_S6C_Table1.txt")
ANKRD1 = read.table("Figure_S6C_Table2.txt")
FLNC = read.table("Figure_S6C_Table3.txt")
KLHL40 = read.table("Figure_S6C_Table4.txt")
MACF1 = read.table("Figure_S6C_Table5.txt")
MYH7 = read.table("Figure_S6C_Table6.txt")
NEB = read.table("Figure_S6C_Table7.txt")
NRAP = read.table("Figure_S6C_Table8.txt")
TNNC1 = read.table("Figure_S6C_Table9.txt")
TNNT1 = read.table("Figure_S6C_Table10.txt")
TTN = read.table("Figure_S6C_Table11.txt")
XIRP2 = read.table("Figure_S6C_Table12.txt")





setwd("") #Insert path to output Figure PDFs


pdf("aRNA_Mammoth_Figure_S6C_ACTA1.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(ACTA1, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(ACTA1$V2), max(ACTA1$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("ACTA1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()


pdf("aRNA_Mammoth_Figure_S6C_ANKRD1.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(ANKRD1, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(ANKRD1$V2), max(ANKRD1$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("ANKRD1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()


pdf("aRNA_Mammoth_Figure_S6C_FLNC.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(FLNC, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(FLNC$V2), max(FLNC$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("FLNC") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()



pdf("aRNA_Mammoth_Figure_S6C_KLHL40.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(KLHL40, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(KLHL40$V2), max(KLHL40$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("KLHL40") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()



pdf("aRNA_Mammoth_Figure_S6C_MACF1.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(MACF1, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(MACF1$V2), max(MACF1$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("MACF1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()



pdf("aRNA_Mammoth_Figure_S6C_MYH7.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(MYH7, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(MYH7$V2), max(MYH7$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("MYH7") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()




pdf("aRNA_Mammoth_Figure_S6C_NEB.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(NEB, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(NEB$V2), max(NEB$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("NEB") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()




pdf("aRNA_Mammoth_Figure_S6C_NRAP.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(NRAP, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(NRAP$V2), max(NRAP$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("NRAP") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()



pdf("aRNA_Mammoth_Figure_S6C_TNNC1.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(TNNC1, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(TNNC1$V2), max(TNNC1$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("TNNC1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()




pdf("aRNA_Mammoth_Figure_S6C_TNNT1.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(TNNT1, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(TNNT1$V2), max(TNNT1$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("TNNT1") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()




pdf("aRNA_Mammoth_Figure_S6C_TTN.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(TTN, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(TTN$V2), max(TTN$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("TTN") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()




pdf("aRNA_Mammoth_Figure_S6C_XIRP2.pdf", height = 25/25.4, width = 41.2/25.4)
ggplot(XIRP2, aes(x=V2, y=log10(V3+1))) +
  geom_line(linewidth=0.2, color= "#5D987B") +
  
  theme_classic() +
  
  xlab("") + ylab("") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  scale_x_continuous(expand=c(0.015, 0.015),
                     breaks = c(min(XIRP2$V2), max(XIRP2$V2)), labels = c("-5'", "-3'")) +
  
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(t = 0.25, r = 0.25, b = -0.25, l = -0.2, "cm")) +
  
  ggtitle("XIRP2") + theme(plot.title = element_text(size = 8)) +
  theme(plot.title = element_text(face="bold.italic"))

dev.off()
