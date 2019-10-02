#!/usr/bin/Rscript
library.path <- .libPaths()
library("ggpubr", lib.loc = library.path)
library("ggplot2", lib.loc = library.path)
library("reshape", lib.loc = library.path)
library("tidyverse", lib.loc = library.path)
library("argparse", lib.loc = library.path)

args <- commandArgs(trailingOnly = TRUE)

# Reading in output CSV file
FeGenie_heatmap_data_organized <- read.csv(args[1], header = FALSE)
l = length(FeGenie_heatmap_data_organized)
fegenie.t <- t(FeGenie_heatmap_data_organized[,3:l-1])
fegenie.matrix = as.matrix(fegenie.t)

colnames(fegenie.t) <- as.vector(FeGenie_heatmap_data_organized$V1)
FeGenie <- fegenie.t

#-------------data preparation
# convert to data frame
FeGenie.data <- as.data.frame(FeGenie)
#class(FeGenie.data)
#(FeGenie.data)

# ----------- melt data
FeGenie.data.melt <- melt(FeGenie.data, id.vars = 1:1)
#head(FeGenie.data.melt)

# rename columns
colnames(FeGenie.data.melt)[colnames(FeGenie.data.melt)=="variable"] <- "Iron_category"
colnames(FeGenie.data.melt)[colnames(FeGenie.data.melt)=="value"] <- "Gene_counts"

# ----------------- output files
outfile = paste(args[2], "Fegenie-dotplot.tiff", sep = "/", collapse = NULL)
tiff(outfile, units="in", width=12, height=10, res=300)

FeGenie.data.melt$Gene_counts = as.character(FeGenie.data.melt$Gene_counts)
FeGenie.data.melt$Gene_counts = as.numeric(FeGenie.data.melt$Gene_counts)

FeGenie.meta.plot <- ggplot(FeGenie.data.melt, aes(x = X, y = Iron_category, size = Gene_counts), alpha=0.7) +
  geom_point(aes(color=Iron_category)) +
  scale_size_area(max_size = 10) +
  labs(x="", y="Iron Category") +
  scale_y_discrete(labels=c("iron_aquisition-iron_uptake" = "Iron uptake", 
                            "iron_aquisition-heme_uptake" = "Heme uptake", 
                            "iron_aquisition-heme_lyase" = "Heme lyase", 
                            "iron_aquisition-siderophore_synthesis" = "Siderophore synthesis",
                            "iron_aquisition-siderophore_uptake" = "Siderophore uptake", 
                            "iron_gene_regulation" = "Iron gene regulation", 
                            "iron_oxidation" = "Iron oxidation", 
                            "iron_reduction" = "Iron reduction", 
                            "iron_storage" = "Iron storage")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
        panel.border = element_rect(colour="black", size=1, fill=NA),
        strip.background=element_rect(fill='white', colour='white'),
        strip.text = element_text(face="bold", size=10),
        panel.grid.major = element_line(size = 0),
        panel.grid.minor = element_line(size = 0),
        axis.text = element_text(size=12, colour="black"),
        axis.title = element_text(face="bold", size=14),
        axis.text.x = element_text(vjust = 1, angle = 45, color = "black", size = 12, hjust=1),
        legend.position="right")

FeGenie.meta.plot

dev.off()



