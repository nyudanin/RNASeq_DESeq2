## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("ggplot2")
library("GMD")
library("pheatmap")

# Convert to matrix----------------------------------------------------
PPGrawcounts <- read.delim("2-2B-1 Raw Counts.txt", row.names=1)
countdata <- as.matrix(PPGrawcounts)
head(countdata)

# Assign conditions ----------------------------------------------------
names <- colnames(countdata)
gettissue <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
tissue <- unlist(lapply(names, gettissue))

getsubset <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
subset <- unlist(lapply(names, getsubset))

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet.
(coldata <- data.frame(row.names=colnames(countdata),tissue,subset))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~tissue+subset)
dds <- dds[ rowSums(counts(dds)) > 1, ]
# Run DESeq normalization ----------------------------------------------------
dds<-DESeq(dds)

write.csv(counts(dds, normalized=TRUE),file="2-2B-1 Normalized Counts.csv")

## Regularized log and variance stabilizing transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds,fast=TRUE)
vsd <- varianceStabilizingTransformation(dds)
write.csv(assay(rld),file="2-2B-1 Log Transformed Normalized Counts.csv")

## PCA ----------------------------------------------------
cairo_pdf("PPG Tissue Subset PCA.pdf", w=6, h=4)
data <- plotPCA(rld, intgroup=c("tissue","subset"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=tissue, shape=subset)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

## PCA (subset only)
cairo_pdf("2-2B-1 Subset PCA.pdf", w=6, h=4)
data <- plotPCA(rld, intgroup=c("subset"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=subset)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

## PCA (tissue only)
cairo_pdf("2-2B-1 Tissue PCA.pdf", w=6, h=4)
data <- plotPCA(rld, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=tissue)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

## Sample distance heatmap ----------------------------------------------------
sampleDists <- as.matrix(dist(t(assay(rld))))
cairo_pdf("PPG Sample Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(sampleDists), 
          scale="row",
          key=F, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-2, 2, length.out = 256)),
          ColSideColors = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(rld)$subset ],
          RowSideColors = c( Spleen="black", Lung="light gray", Jejunum="#707070" )[colData(rld)$tissue ]
          )
dev.off()

## Top Variable Genes Heatmap ----------------------------------------------------
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 50 )
my_palette <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255)
resheatmap(vsd,
           height = 10,
           width = 10,
           genes = topVarGenes,
           samples=-grep("Spleen.*ILC2.*|Lung.*ILC3.*",colnames(vsd)),
           filename="2-2B-1 Subset Signature Genes Heatmap.pdf")
dev.off()

## Differential Expression (tissues)----------------------------------------------------
res.Spleen.Jejunum <- results(dds, contrast=c("tissue","Spleen","Jejunum"))
res.Spleen.Lung <- results(dds, contrast=c("tissue","Spleen","Lung"))
res.Jejunum.Lung <- results(dds,  contrast=c("tissue","Jejunum","Lung"))

#Order by row variance
res.Spleen.Jejunum <- res.Spleen.Jejunum[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.Spleen.Lung <- res.Spleen.Lung[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.Jejunum.Lung <- res.Jejunum.Lung[order(rowVars(assay(vsd)), decreasing=TRUE),]

#Filter top 50 genes by padj<0.01 & LFC >2
sdg.Spleen.Jejunum <- subset(res.Spleen.Jejunum, padj < 0.05 & abs(log2FoldChange)>1.5 & baseMean >10)
sdg.Spleen.Lung <- subset(res.Spleen.Lung, padj < 0.05 & abs(log2FoldChange)>1.5 & baseMean >10)
sdg.Jejunum.Lung <- subset(res.Jejunum.Lung, padj < 0.05 & abs(log2FoldChange)>1.5 & baseMean >10)

## Write DE results
write.csv(sdg.Spleen.Jejunum,file="2-2B-1 Spleen.Jejunum-DE.csv")
write.csv(sdg.Spleen.Lung,file="2-2B-1 Spleen.Lung-DE.csv")
write.csv(sdg.Jejunum.Lung,file="2-2B-1 Jejunum.Lung-DE.csv")

## Tissue Signature Genes ---------------------------------------------------------------------
signatureGenes <- unique(c(rownames(sdg.Spleen.Jejunum),rownames(sdg.Spleen.Lung),rownames(sdg.Jejunum.Lung)))
signatureGenes <- signatureGenes[which(rowMeans(counts(dds)[signatureGenes,])>50, useNames=TRUE)]
signatureGenes <- signatureGenes[order(rowVars(assay(vsd)[signatureGenes,]), decreasing=TRUE)]
write.csv(counts(dds, normalized=TRUE)[signatureGenes,], file="2-2B-1 Tissue Signature Genes.csv")

resheatmap(vsd, 
           genes = rownames(sdg.Spleen.Jejunum),
           samples= -grep("ILC1",colnames(vsd)),
           height = 20,
           filename="2-2B-1 Tissue Signature Genes Heatmap.pdf")
dev.off()

## Sig. DE Genes heatmaps (tissue)----------------------------------------------------

resheatmap(vsd, 
           genes = rownames(sdg.Spleen.Jejunum),
           samples=-grep("Lung|ILC2",colnames(vsd)),
           height=9,
           width=7,
           filename="2-2B-1 Spleen.Jejunum-DEG Heatmap.pdf")
dev.off()

resheatmap(rld, 
           genes = sdg.Spleen.Lung,
           samples=grep("Spleen|Lung",colnames(rld)),
           filename="2-2B-1 Spleen.Lung-DEG Heatmap.pdf")
dev.off()

resheatmap(rld, 
           genes = sdg.Jejunum.Lung,
           samples=grep("Jejunum|Lung",colnames(rld)),
           filename="2-2B-1 Jejunum.Lung-DEG Heatmap.pdf")
dev.off()

## Differential Expression (subsets in each tissue)----------------------------------------------------
res.ILC1.ILC2 <- results(dds, tidy=FALSE, contrast=c("subset","ILC1","ILC2"))
res.ILC1.ILC3 <- results(dds, tidy=FALSE, contrast=c("subset","ILC1","ILC3"))
res.ILC1.NK <- results(dds, tidy=FALSE, contrast=c("subset","ILC1","NK"))
res.ILC2.ILC3 <- results(dds, tidy=FALSE, contrast=c("subset","ILC2","ILC3"))
res.ILC2.NK <- results(dds, tidy=FALSE, contrast=c("subset","ILC2","NK"))
res.ILC3.NK <- results(dds, tidy=FALSE, contrast=c("subset","ILC3","NK"))

#Order by row variance
res.ILC1.ILC2 <- res.ILC1.ILC2[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC1.ILC3 <- res.ILC1.ILC3[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC1.NK <- res.ILC1.NK[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC2.ILC3 <- res.ILC2.ILC3[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC2.NK <- res.ILC2.NK[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC3.NK <- res.ILC3.NK[order(rowVars(assay(vsd)), decreasing=TRUE),]

## Write filtered results
write.csv(sdg.ILC1.ILC2, file="2-2B-1 ILC1.ILC2-DE.csv")
write.csv(sdg.ILC1.ILC3, file="2-2B-1 ILC1.ILC3-DE.csv")
write.csv(sdg.ILC1.NK, file="2-2B-1 ILC1.NK-DE.csv")
write.csv(sdg.ILC2.ILC3, file="2-2B-1 ILC2.ILC3-DE.csv")
write.csv(sdg.ILC2.NK, file="2-2B-1 ILC2.NK-DE.csv")
write.csv(sdg.ILC3.NK, file="2-2B-1 ILC3.NK-DE.csv")


#Filter DE genes by mean counts, padj<0.1 & LFC >2
sdg.ILC1.ILC2 <- subset(res.ILC1.ILC2, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.ILC1.ILC3 <- subset(res.ILC1.ILC3, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10 )
sdg.ILC1.NK <- subset(res.ILC1.NK, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.ILC2.ILC3 <- subset(res.ILC2.ILC3, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10 )
sdg.ILC2.NK <- subset(res.ILC2.NK, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10 )
sdg.ILC3.NK <- subset(res.ILC3.NK, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10 )


## Subset DE Signature Genes ---------------------------------------------------------------------
signatureGenes <- unique(c(rownames(sdg.ILC1.ILC2),
                           rownames(sdg.ILC1.ILC3),
                           rownames(sdg.ILC1.NK),
                           rownames(sdg.ILC2.ILC3),
                           rownames(sdg.ILC2.NK),
                           rownames(sdg.ILC3.NK)))
signatureGenes <- signatureGenes[which(rowMeans(counts(dds)[signatureGenes,])>10, useNames=TRUE)]
signatureGenes <- signatureGenes[order(rowVars(assay(vsd)[signatureGenes,-grep("Lung", colnames(vsd))]), decreasing=TRUE)]
write.csv(counts(dds, normalized=TRUE)[signatureGenes,], file="2-2B-1 Subset Signature Genes.csv")

sigGenesCluster<- resheatmap(vsd, 
                             genes = signatureGenes[1:100], 
                             samples=-grep("Jejunum|Lung|ILC2", colnames(vsd)),
                             height = 14,
                             filename="2-2B-1 Spleen Subset Signature Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = ILC1genes,
           samples=-grep("217|Jejunum|Lung", colnames(vsd)),
           height = 6,
           filename="2-2B-1 ILC1 Signature Genes Heatmap.pdf")

resheatmap(vsd, 
           genes = ILC2genes,
           samples=-grep("217|Jejunum", colnames(vsd)),
           height = 10,
           filename="2-2B-1 ILC2 Signature Genes Heatmap.pdf")

resheatmap(vsd, 
           genes = ILC3genes,
           samples=-grep("217|Jejunum", colnames(vsd)),
           height = 10,
           filename="2-2B-1 ILC3 Signature Genes Heatmap.pdf")

sigGenesCluster<- resheatmap(vsd, 
                             genes = signatureGenes[1:100], 
                             samples=-grep("Spleen|Lung", colnames(vsd)),
                             height = 14,
                             filename="2-2B-1 Jejunum Subset Signature Genes Heatmap.pdf")
dev.off()

#Filter top 50 genes by row variance, padj<0.1 & LFC >2
sdg.ILC1.ILC2 <- rownames(subset(res.ILC1.ILC2, padj < 0.1 & abs(log2FoldChange)>2)[1:50,])
sdg.ILC1.ILC3 <- rownames(subset(res.ILC1.ILC3, padj < 0.1 & abs(log2FoldChange)>2 )[1:50,])
sdg.ILC1.NK <- rownames(subset(res.ILC1.NK, padj < 0.1 & abs(log2FoldChange)>2)[1:50,])
sdg.ILC2.ILC3 <- rownames(subset(res.ILC2.ILC3, padj < 0.1 & abs(log2FoldChange)>2 )[1:50,])
sdg.ILC2.NK <- rownames(subset(res.ILC2.NK, padj < 0.1 & abs(log2FoldChange)>2 )[1:50,])
sdg.ILC3.NK <- rownames(subset(res.ILC3.NK, padj < 0.1 & abs(log2FoldChange)>2 )[1:50,])

## Sig. DE Genes heatmaps (subset)---------------------------------------------------------------------
resheatmap <- function(rld, genes, samples, cluster_cols=TRUE, height = 6, width = 6, filename, ...){
  filtered <- assay(rld) [genes, samples]
  filtered <- filtered[rowVars(filtered)>0,]
  with(rld,
       pheatmap(filtered,
                cluster_rows=TRUE,
                scale="row",
                breaks = c(seq(-2, 2, length.out = 256)),
                border_color = NA,
                drop_levels=TRUE,
                color = my_palette,
                show_rownames=TRUE,
                cluster_cols=cluster_cols,
                annotation_col=coldata,
                annotation_colors = ann_colors,
                annotation_legend=FALSE,
                height=height,
                width=width,
                filename= filename)
       )
}


ann_colors = list(
  subset = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(vsd)$subset ],
  tissue = c( Spleen="black", Lung="light gray", Jejunum="#707070" )[colData(vsd)$tissue ])

resheatmap(vsd, 
           genes = sdg.ILC1.ILC2, 
           samples=grep("ILC1|ILC2",colnames(vsd)), 
           filename="2-2B-1 ILC1.ILC2-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC1.ILC3, 
           samples=grep("ILC1|ILC3",colnames(vsd)), 
           filename="2-2B-1 ILC1.ILC3-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC1.NK, 
           samples=grep("ILC1|NK",colnames(vsd)), 
           filename="2-2B-1 ILC1.NK-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC2.ILC3, 
           samples=grep("ILC2|ILC3",colnames(vsd)), 
           filename="2-2B-1 ILC2.ILC3-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC2.NK, 
           samples=grep("ILC2|NK",colnames(vsd)), 
           filename="2-2B-1 ILC2.NK-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC3.NK, 
           samples=grep("ILC3|NK",colnames(vsd)), 
           filename="2-2B-1 ILC3.NK-DEG Heatmap.pdf")
dev.off()

## Sig DE heatmaps ILC1 vs. NK by tissue-----------
resheatmap(vsd,
           genes = sdg.ILC1.NK,
           samples=grep("Spleen_(.*)ILC1|Spleen_(.*)NK",colnames(vsd)), 
           filename="2-2B-1 Spleen ILC1.NK-DEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = sdg.ILC1.NK, 
           samples=grep("Jejunum_(.*)ILC1|Jejunum_(.*)NK",colnames(vsd)), 
           filename="2-2B-1 Jejunum ILC1.NK-DEG Heatmap.pdf")
dev.off()


## Violin Plots --------------------- 
library(ggplot2)
library(gridExtra)
library(reshape2)

ILC1genes <- as.character(degenelist[,1], na.rm=TRUE)
ILC2genes <- as.character(degenelist[,2], na.rm=TRUE)
ILC3genes <- as.character(degenelist[,3], na.rm=TRUE)

ILC1genes <- ILC1genes[which(ILC1genes %in% rownames(assay(vsd)))]
ILC2genes <- ILC2genes[which(ILC2genes %in% rownames(assay(vsd)))]
ILC3genes <- ILC3genes[which(ILC3genes %in% rownames(assay(vsd)))]
ILC3genes <- ILC3genes[-(which(ILC3genes %in% ILC1genes))]

ILC1counts <- t(counts(dds, normalized=TRUE)[ILC1genes,])[,order(colMeans(ILC1counts), decreasing=TRUE)]
ILC2counts <- t(counts(dds, normalized=TRUE)[ILC2genes,])[,order(colMeans(ILC2counts), decreasing=TRUE)]
ILC3counts <- t(counts(dds, normalized=TRUE)[ILC3genes,])[,order(colMeans(ILC3counts), decreasing=TRUE)]



ILC1counts <- cbind2(coldata, ILC1counts)
ILC2counts <- cbind2(coldata, ILC2counts)
ILC3counts <- cbind2(coldata, ILC3counts)

meltedILC1 <- melt(ILC1counts, variable.name = "gene", value.name = "count" )
meltedILC2 <- melt(ILC2counts, variable.name = "gene", value.name = "count" )
meltedILC3 <- melt(ILC3counts, variable.name = "gene", value.name = "count" )

cairo_pdf(filename="ILC1 Gene Plots.pdf", w=15, h=15)
ggplot(data = meltedILC1,  aes(x = subset, y = log(count+1), fill = subset)) +
  geom_violin() +
  facet_wrap(~gene, nrow=4 ) +
  scale_fill_manual(values = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )) +
  ggtitle("ILC1 Genes") +
  theme(axis.title.x = element_blank(),
        legend.position='none') +
  ylab("Log Normalized Counts")
dev.off()

cairo_pdf(filename="ILC2 Gene Plots.pdf", w=20, h=15)
ggplot(data = meltedILC2[meltedILC2$gene != "CAMP" & meltedILC2$gene != "PTGDR2", ],  aes(x = subset, y = log(count+1), fill = subset)) +
  geom_violin() +
  facet_wrap(~gene, nrow=5 ) +
  scale_fill_manual(values = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )) +
  ggtitle("ILC2 Genes") +
  theme(axis.title.x = element_blank(),
        legend.position='none') +
  ylab("Log Normalized Counts")
dev.off()

cairo_pdf(filename="ILC3 Gene Plots.pdf", w=20, h=20)
ggplot(data = meltedILC3,  aes(x = subset, y = log(count+1), fill = subset)) +
  geom_violin() +
  facet_wrap(~gene, nrow=7 ) +
  scale_fill_manual(values = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )) +
  ggtitle("ILC3 Genes") +
  theme(axis.title.x = element_blank(),
        legend.position='none') +
  ylab("Log Normalized Counts")
dev.off()

ILC1vNKgenes <- as.character(ILC.NK.genes[,1], na.rm=TRUE)
ILC1vNKcounts <- t(counts(dds, normalized=TRUE)[ILC1vNKgenes,])
ILC1vNKcounts <- cbind2(coldata, ILC1vNKcounts)
meltedILC1vNK <- melt(ILC1vNKcounts, variable.name = "gene", value.name = "count" )

cairo_pdf(filename="ILC vs. NK select Plots.pdf", w=5, h=5)
ggplot(data = meltedILC1vNK[meltedILC1vNK$gene == c(NKselect,ILC1select) & meltedILC1vNK$subset != "ILC3",],  aes(x = subset, y = log(count+1), fill = subset)) +
  geom_violin(scale="width") +
  facet_wrap(~gene, nrow=3) +
  scale_fill_manual(values = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )) +
  ggtitle("ILC1 vs. NK Genes") +
  theme(axis.title.x = element_blank(),
        legend.position='none') +
  ylab("Log Normalized Counts")
dev.off()

NKselect <- c("NCF2","CCL4","CCL3","LY9","EOMES","CD79A")
ILC1select <- c("IL2RA","IL4R","IL18R1","IL23R","CCR6","CD28")
