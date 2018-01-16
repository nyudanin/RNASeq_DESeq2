## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("ggplot2")
library("GMD")
library("pheatmap")
library("devtools")
library("factoextra")

date <- paste0(Sys.Date())
expNum <- paste0("2-2B-1")

# Convert to matrix----------------------------------------------------
PPGrawcounts <- read.delim("2-2B-1 Spleen Jejunum Raw Counts.txt", row.names=1)
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

## Heatmap Functions and colors ----------------------------------------------------
normheatmap <- function(mtx, cluster_cols=TRUE, title=title, cex=1, h=1, w=1, ...){
  pheatmap(mtx,
           cex = cex,
           cluster_rows=TRUE,
           scale="row",
           breaks = c(seq(-1.5, 1.5, length.out = 256)),
           border_color = NA,
           drop_levels=TRUE,
           color = my_palette,
           show_rownames=TRUE,
           cluster_cols=cluster_cols,
           annotation_col=coldata[c(1,3)],
           annotation_colors = ann_colors,
           annotation_legend=FALSE,
           main= paste0(expNum, " ",title),
           legend=FALSE,
           fontsize= 10,
           treeheight_col = 20,
           treeheight_row = 20,
           height=h*6,
           width=w*6,
           filename= paste0(date," ",expNum," ",title,".pdf"))
}
resheatmap <- function(vsd, genes, samples, cluster_cols=TRUE, title=title, cex=1, h=1, w=1,  ...){
  filtered <- assay(vsd) [genes, samples]
  filtered <- filtered[rowVars(filtered)>0,]
  with(vsd,
       pheatmap(filtered[order( rowVars( assay(vsd)[genes,]), decreasing=TRUE),],
                cex = cex,
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
                main= paste0(expNum, " ",title),
                legend=FALSE,
                fontsize= 10,
                treeheight_col = 10,
                treeheight_row = 20,
                height=15*h,
                width=6*w,
                filename= paste0(date," ",expNum," ",title,".pdf"))
  )
}

my_palette <- colorRampPalette(brewer.pal(11, "RdBu")) (255)
my_palette <- rev(my_palette)

ann_colors = list(
  subset = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(vsd)$subset],
  tissue = c( Spleen="black", Lung="light gray", Jejunum="#707070" )[colData(vsd)$tissue])

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet.
(coldata <- data.frame(row.names=colnames(countdata),tissue,subset))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~tissue+subset+subset:tissue)
dds <- DESeq(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
resultsNames(dds)


# Variance Transformation ----------------------------------------------------
vsd <- varianceStabilizingTransformation(dds)

cairo_pdf("2-2B-1 Jejunum Spleen Subset PCA.pdf", w=6, h=4)
data <- plotPCA(vsd, intgroup=c("subset"), returnData=FALSE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=subset)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

write.csv( pca$scores[,1:5], file="2-2B-1 ILC1 ILC3 PCA Scores.csv")
plot(pca$scores[order(pca$scores[,3], decreasing=TRUE),3])
quantile(rowVars(pca$scores[,-1]))
pca$scores[which(pca$scores[,1] < -20),1]
pca$scores[,3]

res.ILC1.NK <- results(dds, tidy=FALSE, contrast=c("subset","ILC1","NK"))
res.ILC1.NK <- res.ILC1.NK[order(rowVars(assay(vsd)), decreasing=TRUE),]
sdg.ILC1.NK <- subset(res.ILC1.NK, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
write.csv(sdg.ILC1.NK, file="2-2B-1 ILC1.NK-DE.csv")

resheatmap(vsd,
           genes = ILC.NK.genes[[1]],
           samples = -grep("Spleen|ILC3",colnames(vsd)),
           height = 16,
           filename="2-2B-1 ILC1.NK Functional Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = rownames(sdg.ILC1.NK),
           samples= -grep("ILC3",colnames(vsd)),
           height = 11,
           filename="2-2B-1 ILC1.NK DEG Heatmap.pdf")
dev.off()

# Spleen vs. Jejunum effect for each subset----------------------------------------------------
res.ILC1 <- results(dds, contrast=list("tissueSpleen.subsetILC1","tissueJejunum.subsetILC1"), tidy=FALSE)
res.ILC3 <- results(dds, contrast=list("tissueSpleen.subsetILC3","tissueJejunum.subsetILC3"), tidy=FALSE)
res.NK <- results(dds, contrast=list("tissueSpleen.subsetNK","tissueJejunum.subsetNK"), tidy=FALSE)

res.ILC1 <- res.ILC1[order(res.ILC1$padj, decreasing=FALSE),]
res.ILC3 <- res.ILC3[order(res.ILC3$padj, decreasing=FALSE),]
res.NK <- res.NK[order(res.NK$padj, decreasing=FALSE),]

#Filter significant genes by padj<0.01 & LFC >2----------------------------------------------------
sdg.ILC1 <- subset(res.ILC1, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.ILC3 <- subset(res.ILC3, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.NK <- subset(res.NK, padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)

#Order by Log2FC----------------------------------------------------
sdg.ILC1 <- sdg.ILC1[order(sdg.ILC1$log2FoldChange, decreasing=TRUE),]
sdg.ILC3 <- sdg.ILC3[order(sdg.ILC3$log2FoldChange, decreasing=TRUE),]
sdg.NK <- sdg.NK[order(sdg.NK$log2FoldChange, decreasing=TRUE),]

## Write DE results----------------------------------------------------
write.csv(sdg.ILC1,file="2-2B-1 ILC1.Spleen.Jejunum-DE.csv")
write.csv(sdg.ILC3,file="2-2B-1 ILC3.Spleen.Jejunum-DE.csv")
write.csv(sdg.NK,file="2-2B-1 NK.Spleen.Jejunum-DE.csv")


## Tissue Signature Genes ---------------------------------------------------------------------
signatureGenes <- unique(c(rownames(sdg.ILC1),rownames(sdg.ILC3),rownames(sdg.NK)))
signatureGenes <- signatureGenes[which(signatureGenes %in% rownames(assay(vsd)))]
signatureGenes <- signatureGenes[order(rowVars(assay(vsd)[signatureGenes,]), decreasing=TRUE)]
write.csv(counts(dds, normalized=TRUE)[signatureGenes,], file="2-2B-1 Tissue Signature Genes.csv")

resheatmap(vsd, 
           genes = signatureGenes,
           samples= -grep("ILC1",colnames(vsd)),
           height = 8,
           filename="2-2B-1 Tissue Signature Genes Heatmap.pdf")
dev.off()

## Functional Genes ---------------------------------------------------------------------
functionalGenes <- read.table("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/2-2B-1/2-2B-1 DESeq2/Spleen & Jejunum Samples/2-2B-1 ILC1 vs NK Functional Genes.txt", 
                              quote="\"", comment.char="", stringsAsFactors=FALSE)
functionalGenes <- unlist(functionalGenes[,1])
functionalGenes <- subset(functionalGenes, 
                      rowVars(counts(dds)[functionalGenes,])>1)
functionalGenes <- functionalGenes[order( rowVars( assay(vsd)[functionalGenes,]), decreasing=TRUE)]
functionalGenes <- functionalGenes[which(functionalGenes %in% rownames(assay(vsd)))]
functionalGenes <- functionalGenes[order(rowVars(assay(vsd)[functionalGenes,]), decreasing=TRUE)]

resheatmap(vsd, 
           genes = functionalGenes,
           samples= grep("Spleen.*ILC1|Spleen.*NK",colnames(vsd)),
           height = 13,
           title="Spleen ILC1 vs NK Functional Genes Heatmap")
dev.off()

resheatmap(vsd, 
           genes = functionalGenes,
           samples= grep("Jejunum.*ILC1|Jejunum.*NK",colnames(vsd)),
           height = 13,
           title="Jejunum ILC1 vs NK Functional Genes Heatmap")
dev.off()

resheatmap(vsd, 
           genes = functionalGenes,
           samples= grep(".*ILC1|.*NK",colnames(vsd)),
           height = 13,
           title="Total ILC1 vs NK Functional Genes Heatmap")
dev.off()
