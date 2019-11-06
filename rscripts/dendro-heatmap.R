library.path <- .libPaths()
library("argparse", lib.loc = library.path)
library("ggdendro", lib.loc = library.path)
library("reshape2", lib.loc = library.path)
library("grid", lib.loc = library.path)
library("ggplot2", lib.loc = library.path)
library("pvclust", lib.loc = library.path)

args <- commandArgs(trailingOnly = TRUE)

fegenie <- read.csv(args[1])
l = length(fegenie)
fegenie.sub <- t(fegenie[,2:l])

fegenie.matrix <- as.matrix(fegenie.sub)

colnames(fegenie.matrix) <- as.vector(fegenie$X)


fegenie.scaled = fegenie.matrix
fegenie.scaled <- scale(fegenie.scaled)

fegenie.dendro <- as.dendrogram(hclust(d = dist(x = fegenie.scaled)))
dendro.plot <- ggdendrogram(data = fegenie.dendro, rotate = TRUE)

#outfile = paste(args[2], "Fegenie-dendro.tiff", sep = "/", collapse = NULL)
#tiff(outfile, units="in", width=12, height=10, res=300)
#dendro.plot
#dev.off()

fit <- pvclust(fegenie[,2:l], method.hclust="ward", method.dist="euclidean")
outfile = paste(args[2], "Fegenie-dendro.tiff", sep = "/", collapse = NULL)tiff(outfile, units="in", width=20, height=12, res=300)
tiff(outfile, units="in", width=20, height=12, res=300)
plot(fit, lwd = 4)
pvrect(fit, alpha=.95, lwd = 4)
dev.off()

#********************************************************************************
#***************************** DENDROGRAM *********************************
#********************************************************************************
fegenie.long <- melt(fegenie.scaled, id = c(X1))
fegenie.order <- order.dendrogram(fegenie.dendro)
fegenie.long$Var1 <- factor(x = fegenie.long$Var1,
                          levels = fegenie.long$Var1[fegenie.order], 
                          ordered = TRUE)

heatmap.plot <- ggplot(data = fegenie.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 10))

# REPOSITION LEGEND
heatmap.plot <- ggplot(data = fegenie.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 8),
        legend.position = "top") + theme(axis.text.x = element_text(vjust = 1, angle = 45, color = "black", size = 12, hjust=1))

outfile = paste(args[2], "Fegenie-heatmap.tiff", sep = "/", collapse = NULL)
tiff(outfile, units="in", width=12, height=10, res=300)
heatmap.plot
dev.off()




