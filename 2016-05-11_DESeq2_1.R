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

register(MulticoreParam(4))
setwd("/Volumes/IBD/Yudanin/RNAseq/EC-SM-3509 RNAseq/EC-SM-3509 DESeq2")

## Import data from featureCounts
countdata <- read.delim("/Volumes/IBD/Yudanin/RNAseq/EC-SM-3509 RNAseq/EC-SM-3509 DESeq2/EC-SM-3509 Raw Counts.txt", stringsAsFactors=FALSE)
countdata <- countdata[-(which(duplicated(countdata[,2])==TRUE)),]
row.names(countdata) <- unlist(countdata$Symbol)

## Convert to matrix and remove bad sample #441
countdata <- as.matrix(countdata[,-(1:2)])
countdata <- countdata[,-(grep("X441_NAIVE_KO",colnames(countdata)))]
head(countdata)

## Assign condition 
names <- colnames(countdata)
getsubset <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
subset <- unlist(lapply(names, getsubset))

gettype <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
type <- unlist(lapply(names, gettype))
type <- gsub("CLE","HET",type)

rm(names)

##-------------------------------DESeq2 ----------------------------------------------------
## Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), subset, type))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~subset + type)

### Run DESeq normalization
dds <- DESeq(dds)
write.csv(counts(dds, normalized=TRUE),file="EC-SM-3509 Normalized Counts.csv")


### Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
#rld <- rlogTransformation(dds)

### PCA ----------------------------------------------------
cairo_pdf("EC-SM-3509 PCA.pdf", w=6, h=4)
pca <- plotPCA(vsd, intgroup=c("type","subset"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, shape=subset, color=type)) +
  geom_point(size=4) +
  scale_color_manual(values = c(HET="#4266F6", WT="#269040", KO="#C7302A" )) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

cairo_pdf("EC-SM-3509 WT KO PCA.pdf", w=6, h=4)
pcaWTKO <- plotPCA(, intgroup=c("type","subset"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, shape=subset, color=type)) +
  geom_point(size=4) +
  scale_color_manual(values = c(HET="#4266F6", WT="#269040", KO="#C7302A" )) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
### Sample distance heatmap ---------------------------------------------------------
my_palette <- colorRampPalette(brewer.pal(11, "RdBu")) (255)
my_palette <- rev(my_palette)

sampleDists <- as.matrix(dist(t(assay(vsd))))
nbDists <- as.matrix (dist(t(assay(vsd)[,grep("NB",colnames(vsd))])))
naiveDists <- as.matrix(dist(t(assay(vsd)[,grep("NAIVE",colnames(vsd))])))

cairo_pdf("EC-SM-3509 Sample Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(sampleDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256)),
          ColSideColors = c(HET="#4266F6", WT="#269040", KO="#C7302A" )[colData(vsd)$type],
          RowSideColors = c(NAIVE="#707070", NB="#111111")[colData(vsd)$subset])
dev.off()

cairo_pdf("EC-SM-3509 NB Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(nbDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256)))
dev.off()

cairo_pdf("EC-SM-3509 Naive Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(naiveDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256)))
dev.off()

## Heatmap Function ----------------------------------------------------
resheatmap <- function(vsd, genes, samples, cluster_cols=TRUE, h, w, filename, ...){
  filtered <- assay(vsd) [genes, samples]
  filtered <- filtered[rowVars(filtered)>0,]
  with(vsd,
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
                height=h,
                width=w,
                filename= filename)
  )
}

ann_colors = list(
  type = c(HET="#4266F6", WT="#269040", KO="#C7302A" )[colData(vsd)$type],
  subset = c(NAIVE="#707070", NB="#111111")[colData(vsd)$subset])

### Top Variable Genes Heatmap ----------------------------------------------------
topVarGenes <- head( order( rowVars( assay(vsd)), decreasing=TRUE ), 100 )
nbVarGenes <- head( order( rowVars( assay(vsd)[,grep("NB",colnames(vsd))]), decreasing=TRUE ), 100 )
naiveVarGenes <- head( order( rowVars( assay(vsd)[,grep("NAIVE",colnames(vsd))]), decreasing=TRUE ), 100 )

resheatmap(vsd, 
           genes = topVarGenes,
           samples= colnames(vsd),
           h=15,
           w=6,
           filename="EC-SM-3509 Top Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = nbVarGenes,
           samples= grep("NB",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-SM-3509 NB Top Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = naiveVarGenes,
           samples= grep("NAIVE",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-SM-3509 Naive Top Variable Genes Heatmap.pdf")
dev.off()
### Differential Expression  ----------------------------------------------------

resultsNames(dds)
res.NB.NAIVE <- results(dds, contrast=c("subset","NB","NAIVE"))
res.WT.KO <- results(dds, contrast=c("type","WT","KO"))
res.WT.HET <- results(dds, contrast=c("type","WT","HET"))
res.KO.HET <- results(dds, contrast=c("type","KO","HET"))

#Order by Row Variance
res.NB.NAIVE <- res.NB.NAIVE[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.WT.KO <- res.WT.KO[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.WT.HET <- res.WT.HET[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.KO.HET <- res.KO.HET[order(rowVars(assay(vsd)), decreasing=TRUE),]

#Filter DE genes by mean counts, padj<0.1 & LFC >2----------------------------------------------------
sdg.NB.NAIVE <- subset(res.NB.NAIVE[order(res.NB.NAIVE$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.WT.KO <- subset(res.WT.KO[order(res.WT.KO$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.WT.HET <- subset(res.WT.HET[order(res.WT.HET$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.KO.HET <- subset(res.KO.HET[order(res.KO.HET$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)

# Write DE results
write.csv(res.NB.NAIVE,file="EC-SM-3509 NB vs Naive DEG.csv")
write.csv(res.WT.KO,file="EC-SM-3509 WT vs KO DEG.csv")
write.csv(res.WT.HET,file="EC-SM-3509 WT vs HET DEG.csv")
write.csv(res.KO.HET,file="EC-SM-3509 KO vs HET DEG.csv")

save.image("/Volumes/IBD/Yudanin/RNAseq/EC-SM-3509 RNAseq/EC-SM-3509 DESeq2/EC-SM-3509 DESeq2 Combined.RData")

### Individual Group DESeq ---------------------------------------------------------------------
dds$group <- factor(paste0(dds$subset, dds$type))
design(dds) <- ~ group
dds <- DESeq(dds)
save.image("/Volumes/IBD/Yudanin/RNAseq/EC-SM-3509 RNAseq/EC-SM-3509 DESeq2/EC-SM-3509 DESeq2 Individual.RData")

resultsNames(dds)
write.csv(counts(dds, normalized=TRUE),file="EC-SM-3509 Individual Group Normalized Counts.csv")

vsd <- varianceStabilizingTransformation(dds)

cairo_pdf("EC-SM-3509 Individual Group PCA.pdf", w=6, h=4)
pca <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(pca, aes(PC1, PC2, color=group)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

### NB  Differential Expression ---------------------------------------------------------------------
res.NB.KO.WT <- results(dds, contrast=c("group","NBWT","NBKO"))
res.NB.KO.WT <- res.NB.KO.WT[order(res.NB.KO.WT$baseMean, decreasing = TRUE),]
sdg.NB.KO.WT <- subset(res.NB.KO.WT[order(res.NB.KO.WT$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean > 10)
write.csv(res.NB.KO.WT,file="EC-SM-3509 NB KO.WT DEG.csv")
write.csv(sdg.NB.KO.WT,file="EC-SM-3509 NB KO.WT SigDEG.csv")

res.NB.HET.WT <- results(dds, contrast=c("group","NBWT","NBHET"))
res.NB.HET.WT <- res.NB.HET.WT[order(res.NB.HET.WT$baseMean, decreasing = TRUE),]
sdg.NB.HET.WT <- subset(res.NB.HET.WT[order(res.NB.HET.WT$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean >10)
write.csv(res.NB.HET.WT,file="EC-SM-3509 NB HET.WT DEG.csv")
write.csv(sdg.NB.HET.WT,file="EC-SM-3509 NB HET.WT SigDEG.csv")

res.NB.KO.HET <- results(dds, contrast=c("group","NBKO","NBHET"))
res.NB.KO.HET <- res.NB.KO.HET[order(res.NB.KO.HET$baseMean, decreasing = TRUE),]
sdg.NB.KO.HET <- subset(res.NB.KO.HET[order(res.NB.KO.HET$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean >10)
write.csv(res.NB.KO.HET,file="EC-SM-3509 NB KO.HET DEG.csv")
write.csv(sdg.NB.KO.HET,file="EC-SM-3509 NB KO.HET SigDEG.csv")

resheatmap(vsd, 
           genes = row.names(sdg.NB.KO.WT), 
           samples=grep("NB", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NB.KO.WT)),
           w=4,
           filename="EC-SM-3509 NB KO.WT SigDEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.NB.HET.WT), 
           samples=grep("NB", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NB.HET.WT)),
           w=4,
           filename="EC-SM-3509 NB HET.WT SigDEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.NB.KO.HET), 
           samples=grep("NB", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NB.KO.HET)),
           w=4,
           filename="EC-SM-3509 NB HET.KO SigDEG Heatmap.pdf")
dev.off()

### NAIVE Differential Expression ---------------------------------------------------------------------
res.NAIVE.KO.WT <- results(dds, contrast=c("group","NAIVEWT","NAIVEKO"))
res.NAIVE.KO.WT <- res.NAIVE.KO.WT[order(res.NAIVE.KO.WT$baseMean, decreasing = TRUE),]
sdg.NAIVE.KO.WT <- subset(res.NAIVE.KO.WT[order(res.NAIVE.KO.WT$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean > 10)
write.csv(res.NAIVE.KO.WT,file="EC-SM-3509 NAIVE KO.WT DEG.csv")
write.csv(sdg.NAIVE.KO.WT,file="EC-SM-3509 NAIVE KO.WT SigDEG.csv")

res.NAIVE.HET.WT <- results(dds, contrast=c("group","NAIVEWT","NAIVEHET"))
res.NAIVE.HET.WT <- res.NAIVE.HET.WT[order(res.NAIVE.HET.WT$baseMean, decreasing = TRUE),]
sdg.NAIVE.HET.WT <- subset(res.NAIVE.HET.WT[order(res.NAIVE.HET.WT$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean >10)
write.csv(res.NAIVE.HET.WT,file="EC-SM-3509 NAIVE HET.WT DEG.csv")
write.csv(sdg.NAIVE.HET.WT,file="EC-SM-3509 NAIVE HET.WT SigDEG.csv")

res.NAIVE.KO.HET <- results(dds, contrast=c("group","NAIVEKO","NAIVEHET"))
res.NAIVE.KO.HET <- res.NAIVE.KO.HET[order(res.NAIVE.KO.HET$baseMean, decreasing = TRUE),]
sdg.NAIVE.KO.HET <- subset(res.NAIVE.KO.HET[order(res.NAIVE.KO.HET$baseMean, decreasing = TRUE),] , abs(log2FoldChange) > 1  & baseMean >10)
write.csv(res.NAIVE.KO.HET,file="EC-SM-3509 NAIVE KO.HET DEG.csv")
write.csv(sdg.NAIVE.KO.HET,file="EC-SM-3509 NAIVE KO.HET SigDEG.csv")

resheatmap(vsd, 
           genes = row.names(sdg.NAIVE.KO.WT), 
           samples=grep("NAIVE", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NAIVE.KO.WT)),
           w=4,
           filename="EC-SM-3509 NAIVE KO.WT SigDEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.NAIVE.HET.WT), 
           samples=grep("NAIVE", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NAIVE.HET.WT)),
           w=4,
           filename="EC-SM-3509 NAIVE HET.WT SigDEG Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.NAIVE.KO.HET), 
           samples=grep("NAIVE", colnames(vsd)),
           h= 0.17 * length(row.names(sdg.NAIVE.KO.HET)),
           w=4,
           filename="EC-SM-3509 NAIVE HET.KO SigDEG Heatmap.pdf")
dev.off()
save.image("/Volumes/IBD/Yudanin/RNAseq/EC-SM-3509 RNAseq/EC-SM-3509 DESeq2/EC-SM-3509 DESeq2 Individual.RData")