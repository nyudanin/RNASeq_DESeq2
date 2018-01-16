## Load required libraries
library("BiocParallel")
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("ggplot2")
library("gplots")
library("pheatmap")
library("Rtsne")

## SETUP ## -------------------
register(MulticoreParam(4))
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/3-1E-1/3-1E-1 RProject")
date <- paste0(Sys.Date())
expNum <- paste0("3-1E-1")

colnames(countdata)
colnames(coldata)
coldata$dietGroup <- as.factor(paste0(coldata$Diet.1,"_",coldata$Diet.2))
colorder <- c(grep("ILC2", colnames(dds)), grep("Treg", colnames(dds)))

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~cell.type+dietGroup)
dds<-DESeq(dds)

write.csv(counts(dds, normalized=TRUE)[,colorder],file= paste0(date," ",expNum," ", "Normalized Counts.csv"))
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

### Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
rld <- rlogTransformation(dds)
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

### PCA ----------------------------------------------------
pca <- plotPCA(vsd, intgroup=c("dietGroup", "cell.type"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

dietgroupfill <- c("Con_Con" = "#1f78b4", 
                   "Con_Fat" = "#a6cee3", 
                   "Fat_Con" = "#b2df8a", 
                   "Fat_Fat" = "#33a02c") [pca$dietGroup]

cell.typeshapes <- c(Treg=21, ILC2=22) [pca$cell.type]
cell.typecolors <- c(Treg="#D62D20", ILC2="#212121") [pca$cell.type]

pdf(paste0(date," ",expNum," Subset & Diet PCA.pdf"), w=6, h=6)
ggplot(pca, aes(PC1, PC2, shape=cell.type, color=dietGroup)) +
  geom_point(size=4, stroke = 1, aes(shape=cell.type, color=cell.type, fill=dietGroup)) +
  scale_shape_manual(values = cell.typeshapes) +
  scale_color_manual(values = cell.typecolors) +
  scale_fill_manual(values = dietgroupfill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste0(expNum," Subset & Diet PCA"))
dev.off()

pca2 <- prcomp(assay(dds), scale=TRUE)
write.csv(pca2$x, file = "3-1E-1 PCA Genes.csv")
write.csv(pca2$rotation, file = "3-1E-1 PCA Samples.csv")

pdf(paste0(date," ",expNum," Subset & Diet PCA2.pdf"), w=7, h=7)
plot(pca2$rotation[,c(1,2)], bg=dietgroupfill, pch=cell.typeshapes, col=cell.typecolors, cex=2)
dev.off()


