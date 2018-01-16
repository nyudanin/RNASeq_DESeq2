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
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/2-1A-4/2-1A-4 DESeq2")
date <- paste0(Sys.Date())
expNum <- paste0("2-1A-4")
colorder <- c(grep("Colon", colnames(dds)), grep("SI", colnames(dds)), grep("MWAT", colnames(dds)), grep("Lung", colnames(dds)))

# Import data from featureCounts
countdata <- read.table("2-1A-4 Raw Counts.txt", header=TRUE, row.names=1)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition 
names <- colnames(countdata)
getday <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
day <- unlist(lapply(names, getday))

gettissue <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
tissue <- unlist(lapply(names, gettissue))

getreplicate <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
replicate <- unlist(lapply(names, getreplicate))
rm(names)

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), tissue, replicate, day))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~tissue+replicate+day)

## Run DESeq normalization
dds<-DESeq(dds)
write.csv(counts(dds, normalized=TRUE)[,colorder],file= paste0(date," ",expNum," ", "Normalized Counts.csv"))
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

## Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

## PCA ----------------------------------------------------
pca <- plotPCA(vsd, intgroup=c("tissue","day", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

tissuefill <- c(Colon="#E7A626", SI="#4266F6", MWAT="#269040", Lung="#C7302A" ) [pca$tissue]
dayshapes <- c(D1=21, D2=22, D3=25) [pca$day]
replicatecolor <- c(A="#212121", B="#000000") [pca$replicate]

cairo_pdf(paste0(date," ",expNum," Tissue, Day & Replicate PCA.pdf"), w=6, h=4)
ggplot(pca, aes(PC1, PC2, shape=day, color=replicate)) +
  geom_point(size=4, stroke = 0.5, aes(shape=day, color=replicate, fill=tissue)) +
  scale_shape_manual(values = dayshapes) +
  scale_color_manual(values = replicatecolor) +
  scale_fill_manual(values = tissuefill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste0(expNum," Tissue, Day & Replicate PCA.pdf")) 
dev.off()

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
           annotation_col=coldata[1:3],
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
resheatmap <- function(vsd, genes, samples, cluster_cols=TRUE, title=title, cex, h=1, w=1,  ...){
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
                annotation_col=coldata[1:3],
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
  day = c(D1="#C7302A", D2="#707070", D3="#4266F6" )[colData(vsd)$day],
  replicate = c(A="#959595", B="#000000")[colData(vsd)$replicate],
  tissue = c(Colon="#E7A626", SI="#4266F6", MWAT="#269040", Lung="#C7302A" )[colData(vsd)$tissue] )

## Sample distance heatmap ---------------------------------------------------------
sampleDists <- as.matrix(dist(t(assay(vsd))))

normheatmap(sampleDists,
            title= "Sample Distance Heatmap")

## tSNE ----------------------------------------------------
rtsnesample <- Rtsne(sampleDists, is_distance = TRUE, perplexity = 5)

tissuefill <- c(Colon="#E7A626", SI="#4266F6", MWAT="#269040", Lung="#C7302A" ) [colData(vsd)$tissue]
dayshapes <- c(D1=21, D2=22, D3=25) [colData(vsd)$day]
replicatecolor <- c(A="#212121", B="#000000") [colData(vsd)$replicate]

cairo_pdf(paste0(date," ",expNum," Tissue, Day & Replicate tSNE.pdf"), w=6, h=4)
plot(rtsnesample$Y, col=replicatecolor, bg=tissuefill, pch=dayshapes, cex=1, xlab="tSNE 1", ylab="tSNE 2", main="2-1A-4 Tissue, Day & Replicate tSNE")
dev.off()

## Top Variable Genes Heatmap ----------------------------------------------------
minNOTzero <- rownames(assay(vsd)[which(rowMin(counts(dds))>30),])
minNOTzero <- minNOTzero[which(rowVars(assay(vsd)[minNOTzero,])>1)]

topVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero,]), decreasing=TRUE)], 100)
colonVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("Colon",colnames(vsd))]), decreasing=TRUE)], 100)
siVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("SI",colnames(vsd))]), decreasing=TRUE)], 100)
mwatVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("MWAT",colnames(vsd))]), decreasing=TRUE)], 100)
lungVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("Lung",colnames(vsd))]), decreasing=TRUE)], 100)
allVarGenes <- unique(c(topVarGenes, colonVarGenes, siVarGenes, mwatVarGenes, lungVarGenes))

write.csv(counts(dds, normalized=TRUE)[allVarGenes,colorder],file= paste0(date," ",expNum," ", "Top Variable Genes.csv"))

resheatmap(vsd, 
           genes = topVarGenes,
           samples= colnames(vsd),
           cex= 1,
           title= "Top Variable Genes Heatmap")
dev.off()

resheatmap(vsd, 
           genes = allVarGenes,
           samples= colnames(vsd),
           cex= 0.9,
           title= "All Top Variable Genes Heatmap")
dev.off()

### Group Design DESeq Analysis ---------------------------------------------------------------------
dds$group <- factor(paste0(dds$day,"_",dds$tissue,"_",dds$replicate))
design(dds) <- ~ group
ddsGroup <- DESeq(dds)
vsdGroup <- varianceStabilizingTransformation(ddsGroup)
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

## Differential Expression  ----------------------------------------------------
res.SI.Colon <- results(dds, contrast=c("tissue","SI","Colon"))
res.SI.Lung <- results(dds, contrast=c("tissue","SI","Lung"))
res.SI.MWAT <- results(dds, contrast=c("tissue","SI","MWAT"))
res.Colon.Lung <- results(dds, contrast=c("tissue","Colon","Lung"))
res.Colon.MWAT <- results(dds, contrast=c("tissue","Colon","MWAT"))
res.Lung.MWAT <- results(dds, contrast=c("tissue","Lung","MWAT"))

#Order by Row Variance
res.SI.Colon <- res.SI.Colon[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.SI.Lung <- res.SI.Lung[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.SI.MWAT <- res.SI.MWAT[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.Colon.Lung <- res.Colon.Lung[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.Colon.MWAT <- res.Colon.MWAT[order(rowVars(assay(vsd)), decreasing=TRUE),]
res.Lung.MWAT <- res.Lung.MWAT[order(rowVars(assay(vsd)), decreasing=TRUE),]

#Filter DE genes by mean counts, padj<0.1 & LFC >2----------------------------------------------------
sdg.SI.Colon <- subset(res.SI.Colon[order(res.SI.Colon$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)
sdg.SI.Lung <- subset(res.SI.Lung[order(res.SI.Lung$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)
sdg.SI.MWAT <- subset(res.SI.MWAT[order(res.SI.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)
sdg.Colon.Lung <- subset(res.Colon.Lung[order(res.Colon.Lung$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)
sdg.Colon.MWAT <- subset(res.Colon.MWAT[order(res.Colon.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)
sdg.Lung.MWAT <- subset(res.Lung.MWAT[order(res.Lung.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >100)


## Tissue DE Genes Heatmap ---------------------------------------------------------------------
sigDEGenes <- unique(c(rownames(sdg.SI.Colon),
                           rownames(sdg.SI.Lung),
                           rownames(sdg.SI.MWAT),
                           rownames(sdg.Colon.Lung),
                           rownames(sdg.Colon.MWAT),
                           rownames(sdg.Lung.MWAT)))
sigDEGenes <- subset(sigDEGenes, 
                     rowMin(counts(dds)[sigDEGenes,])>50)
sigDEGenes <- sigDEGenes[order(rowVars(assay(vsd)[sigDEGenes, ]), decreasing=TRUE)]
write.csv(counts(dds, normalized=TRUE)[sigDEGenes,colorder],file=paste0(date," ",expNum," Significant Differentially Expressed Genes.csv" ))

resheatmap(vsd, 
           genes = sigDEGenes,
           samples= colnames(vsd),
           cex= 0.9,
           title= "Significant DE Genes Heatmap")
dev.off()

## Comparison to 3-1A-3 --------------------- 
signaturegenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 DAVID/3-1A-3 signaturegenes.txt", 
                             quote="\"", comment.char="", na.strings="", stringsAsFactors=FALSE)
signaturegenes <- signaturegenes[[1]]
signaturegenes <- signaturegenes[which(signaturegenes %in% rownames(counts(dds)))]
signaturegenes <- signaturegenes[order( rowVars( assay(vsd)[signaturegenes,]), decreasing=TRUE)]
signaturegenes <- subset(signaturegenes, 
                     rowMin(counts(dds)[signaturegenes,])>10)
signaturegenes <- signaturegenes[order(rowVars(assay(vsd)[signaturegenes, ]), decreasing=TRUE)]

resheatmap(vsd, 
           genes = signaturegenes,
           samples= colnames(vsd),
           cex= 0.9,
           title= "3-1A-3 Comparison Heatmap")
dev.off()




