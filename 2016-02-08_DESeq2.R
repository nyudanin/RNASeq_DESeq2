## Load required libraries
library("DESeq2")
library("gplots")
library("ggplot2")
library("pheatmap")
library("Rtsne")
library("RColorBrewer")
library("vsn")
library("PoiClaClu")

load("/Users/Naomi/Google Drive/Weill Cornell/Research/Data Analysis/RNAseq/5-1A-1/600SingleCellsMatrix.RData")

# Convert to matrix
countdata <- as.matrix(data_all)
write.csv(countdata, file="600 Splenocytes.csv")
# Assign conditions 
names <- colnames(countdata)
getgroup <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
group <- unlist(lapply(names, getgroup))

getqual <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
qual <- unlist(lapply(names, getqual))

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet.
(coldata <- data.frame(row.names=colnames(countdata), group, qual))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group+qual)
nt <- log2(counts(dds)+1)
se <- SummarizedExperiment(nt, colData=colData(dds))

## PCA and tSNE ----------------------------------------------------
sampleDists <- as.matrix(dist(t(nt)[,1:300]))
rtsnesample <- Rtsne(sampleDists, is_distance = TRUE)

qualln <- length(levels(colData(dds)$qual))

rtsnegene <- Rtsne(unique(nt))
groupcolor <- c(grp1="#C7302A", grp2="#4266F6")[colData(dds)$group]
qualcolor <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(qualln)
my_palette <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(100)

gene <- nt[topVarGenes,]
gene <- nt/max(nt)
gene <- round(gene*100)

genecolor <- my_palette[gene]

cairo_pdf("5-1A-1 Group tSNE.pdf", w=6, h=6)
plot(rtsnesample$Y, col=groupcolor, pch=19, cex=0.5, xlab="tSNE 1", ylab="tSNE 2", main="5-1A-1 Group tSNE")
dev.off()

cairo_pdf("5-1A-1 Quality tSNE.pdf", w=6, h=6)
plot(rtsnesample$Y, col=qualcolor, pch=19, cex=0.5, xlab="tSNE 1", ylab="tSNE 2", main="5-1A-1 Quality tSNE")
dev.off()

cairo_pdf("5-1A-1 Quality tSNE.pdf", w=6, h=6)
plot(rtsnesample$Y, col=genecolor, pch=19, cex=0.5, xlab="tSNE 1", ylab="tSNE 2", main="tSNE")
dev.off()

sepca <- plotPCA (DESeqTransform(se), intgroup=c("qual"))

pca <- plot(princomp(nt[,1:100]))
pcascores <- pca$scores
write.csv(pcascores, file="5-1A-1 PCA Scores.csv")

cairo_pdf("5-1A-1 Sample PCA.pdf", w=6, h=6)
plotPCA (DESeqTransform(se), intgroup=c("qual"))
dev.off()

## Sample distance heatmap ----------------------------------------------------
sampleDists <- as.matrix(dist(t(nt)))
my_palette <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(100)

cairo_pdf("5-1A-1 Sample Distance Heatmap.pdf", w=6, h=6)
pheatmap (sampleDists,
         cluster_rows=TRUE,
         scale="row",
         breaks = c(seq(-1, 1, length.out = 99)),
         border_color = NA,
         drop_levels=TRUE,
         color = my_palette,
         show_rownames=FALSE,
         show_colnames = FALSE,
         cluster_cols=TRUE,
         annotation_legend=FALSE,
         height=h,
         width=w)
dev.off()



## Top PCA Genes Heatmap ----------------------------------------------------

pcagenes <- read.table("~/Google Drive/Weill Cornell/Research/Data Analysis/RNAseq/5-1A-1/5-1A-1 PCA Genes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
pcagenes <- pcagenes$V1
topPCAgenes <- pcagenes[1:50]

cairo_pdf("5-1A-1 PCA Genes Heatmap.pdf", w=6, h=6)
pheatmap (nt[topPCAgenes,],
          cluster_rows=TRUE,
          scale="row",
          border_color = NA,
          drop_levels=TRUE,
          color = my_palette,
          show_rownames=TRUE,
          show_colnames = FALSE,
          cluster_cols=TRUE,
          annotation_legend=FALSE,
          height=h,
          width=w)
dev.off()

ann_colors = c(
  group = groupcolor,
  qual = qualcolor)



