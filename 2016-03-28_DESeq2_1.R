## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("ggplot2")
library("GMD")
library("pheatmap")
library("vsn")
library("Rtsne")



# Import data from featureCounts
countdata <- read.delim("~/Google Drive/Weill Cornell/Research/Data/RNAseq/Other/EC-LR-3455/EC-LR-3455 Raw Counts.txt", stringsAsFactors=FALSE, header=TRUE)
countdata <- countdata[-(which(duplicated(countdata[,1])==TRUE)),]
row.names(countdata) <- unlist(countdata$Symbol)

# Convert to matrix
countdata <- as.matrix(countdata[,-1])
head(countdata)

# Assign condition 
names <- colnames(countdata)
getsubset <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
subset <- unlist(lapply(names, getsubset))

gettype <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
type <- unlist(lapply(names, gettype))

rm(names)

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), subset, type))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~subset+type+subset:type)

## Run DESeq normalization
dds<-DESeq(dds)
write.csv(counts(dds, normalized=TRUE),file="EC-LR-3455 Normalized Counts.csv")

## Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
#rld <- rlogTransformation(dds)

## PCA ----------------------------------------------------
cairo_pdf("EC-LR-3455 Subset PCA.pdf", w=6, h=4)
data <- plotPCA(vsd, intgroup=c("type","subset", "subset"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=subset)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

## Sample distance heatmap ---------------------------------------------------------
my_palette <- colorRampPalette(brewer.pal(11, "RdBu")) (255)
my_palette <- rev(my_palette)

sampleDists <- as.matrix(dist(t(assay(vsd))))
ILC2Dists <- as.matrix (dist(t(assay(vsd)[,grep("ILC2",colnames(vsd))])))
ILC3Dists <- as.matrix (dist(t(assay(vsd)[,grep("ILC3",colnames(vsd))])))
GFDists <- as.matrix (dist(t(assay(vsd)[,grep("GF",colnames(vsd))])))
SPFDists <- as.matrix (dist(t(assay(vsd)[,grep("SPF",colnames(vsd))])))

cairo_pdf("EC-LR-3455 Sample Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(sampleDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256)),
          ColSideColors = c(ILC2="#E7A626", ILC3="#4266F6")[colData(vsd)$subset],
          RowSideColors = c(GF="#269040", SPF="#C7302A")[colData(vsd)$type])
dev.off()

cairo_pdf("EC-LR-3455 ILC2 Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(ILC2Dists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256))
          )
dev.off()

cairo_pdf("EC-LR-3455 ILC3 Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(ILC3Dists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256))
)
dev.off()

cairo_pdf("EC-LR-3455 GF Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(GFDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256))
)
dev.off()

cairo_pdf("EC-LR-3455 SPF Distance Heatmap.pdf", w=6, h=6)
heatmap.2(as.matrix(SPFDists), 
          scale="row",
          key=FALSE, 
          trace="none",
          col=my_palette,
          margin=c(10, 10),
          breaks = c(seq(-1.5, 1.5, length.out = 256))
)
dev.off()

# Heatmap Function ----------------------------------------------------
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
  type = c(GF="#269040", SPF="#C7302A")[colData(vsd)$type],
  subset = c(ILC2="#E7A626", ILC3="#4266F6")[colData(vsd)$subset])

## Top Variable Genes Heatmap ----------------------------------------------------
topVarGenes <- head( order( rowVars( assay(vsd)), decreasing=TRUE ), 100 )
ILC2VarGenes <- head( order( rowVars( assay(vsd)[,grep("ILC2",colnames(vsd))]), decreasing=TRUE ), 100 )
ILC3VarGenes <- head( order( rowVars( assay(vsd)[,grep("ILC3",colnames(vsd))]), decreasing=TRUE ), 100 )
SPFVarGenes <- head( order( rowVars( assay(vsd)[,grep("SPF",colnames(vsd))]), decreasing=TRUE ), 100 )
GFVarGenes <- head( order( rowVars( assay(vsd)[,grep("GF",colnames(vsd))]), decreasing=TRUE ), 100 )

resheatmap(vsd, 
           genes = topVarGenes,
           samples= colnames(vsd),
           h=15,
           w=6,
           filename="EC-LR-3455 Top Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = ILC2VarGenes,
           samples= grep("ILC2",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-LR-3455 ILC2 Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = ILC3VarGenes,
           samples= grep("ILC3",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-LR-3455 ILC3 Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = SPFVarGenes,
           samples= grep("SPF",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-LR-3455 SPF Variable Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = GFVarGenes,
           samples= grep("GF",colnames(vsd)),
           h=15,
           w=4,
           filename="EC-LR-3455 GF Variable Genes Heatmap.pdf")
dev.off()

## Differential Expression  ----------------------------------------------------
res.GF.SPF <- results(dds, contrast=c("type","GF","SPF"))
res.ILC2.ILC3 <- results(dds, contrast=c("subset","ILC2","ILC3"))

#Order by Row Variance
res.GF.SPF <- res.GF.SPF[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.ILC2.ILC3 <- res.ILC2.ILC3[order(rowVars(assay(vsd)), decreasing=TRUE),]

#Filter DE genes by mean counts, padj<0.1 & LFC >2----------------------------------------------------
sdg.GF.SPF <- subset(res.GF.SPF[order(res.GF.SPF$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.ILC2.ILC3 <- subset(res.ILC2.ILC3[order(res.ILC2.ILC3$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)

## Write DE results
write.csv(res.GF.SPF,file="EC-LR-3455 GF.SPF-DE.csv")
write.csv(res.ILC2.ILC3,file="EC-LR-3455 ILC2.ILC3-DE.csv")

## Subset DE Genes Heatmap ---------------------------------------------------------------------
sigDEGenes <- row.names(sdg.ILC2.ILC3)
resheatmap(vsd, 
          genes = sigDEGenes[1:100], 
          samples=colnames(vsd),
          h=15,
          w=6,
          filename="EC-LR-3455 ILC2.ILC3 DE Genes Heatmap.pdf")
dev.off()

## Individual Group DESeq ---------------------------------------------------------------------
dds$group <- factor(paste0(dds$subset, dds$type))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

write.csv(counts(dds, normalized=TRUE),file="EC-LR-3455 Individual Group Normalized Counts.csv")
vsd <- varianceStabilizingTransformation(dds)

cairo_pdf("EC-LR-3455 Individual Group PCA.pdf", w=6, h=4)
data <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=group)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

res.ILC2 <- results(dds, contrast=c("group","ILC2GF","ILC2SPF"))
res.ILC2 <- res.ILC2[order(rowVars(assay(vsd)), decreasing=TRUE),]
sdg.ILC2 <- subset(res.ILC2[order(res.ILC2$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
write.csv(sdg.ILC2,file="EC-LR-3455 ILC2-DE.csv")

res.ILC3 <- results(dds, contrast=c("group","ILC3GF","ILC3SPF"))
res.ILC3 <- res.ILC3[order(rowVars(assay(vsd)), decreasing=TRUE),]
sdg.ILC3 <- subset(res.ILC3[order(res.ILC3$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
write.csv(sdg.ILC3,file="EC-LR-3455 ILC3-DE.csv")

res.GF<- results(dds, contrast=c("group","ILC2GF","ILC3GF"))
res.GF <- res.GF[order(rowVars(assay(vsd)), decreasing=TRUE),]
sdg.GF <- subset(res.GF[order(res.GF$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
write.csv(sdg.GF,file="EC-LR-3455 GF-DE.csv")

res.SPF<- results(dds, contrast=c("group","ILC2SPF","ILC3SPF"))
res.SPF <- res.SPF[order(rowVars(assay(vsd)), decreasing=TRUE),]
sdg.SPF <- subset(res.SPF[order(res.SPF$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
write.csv(sdg.SPF,file="EC-LR-3455 SPF-DE.csv")


resheatmap(vsd, 
           genes = row.names(sdg.ILC2), 
           samples=colnames(vsd),
           h=10,
           w=6,
           filename="EC-LR-3455 ILC2 DE Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.ILC3), 
           samples=colnames(vsd),
           h=10,
           w=6,
           filename="EC-LR-3455 ILC3 DE Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.GF)[1:100], 
           samples=colnames(vsd),
           h=15,
           w=6,
           filename="EC-LR-3455 GF DE Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = row.names(sdg.SPF)[1:100], 
           samples=colnames(vsd),
           h=15,
           w=6,
           filename="EC-LR-3455 SPF DE Genes Heatmap.pdf")
dev.off()
