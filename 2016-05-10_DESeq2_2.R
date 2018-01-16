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
dir <- setwd("/Volumes/IBD/Yudanin/RNAseq/2-2B-1 RNAseq/2-2B-1 DESeq2")
date <- paste0(Sys.Date())
expNum <- paste0("2-2B-1")


# Import data from featureCounts
countdata <- read.table("2-2B-1 Raw Counts.txt", header=TRUE, row.names=1)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition 
names <- colnames(countdata)
getdonor <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
donor <- unlist(lapply(names, getdonor))

gettissue <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
tissue <- unlist(lapply(names, gettissue))

getsubset <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
subset <- unlist(lapply(names, getsubset))
rm(names)

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), tissue, donor, subset))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~tissue+donor+subset)
colorder <- c(grep("NK", colnames(dds), grep("ILC1", colnames(dds)), grep("ILC2", colnames(dds)), grep("ILC3", colnames(dds))))

## Run DESeq normalization
dds<-DESeq(dds)
write.csv(counts(dds, normalized=TRUE)[,colorder],file= paste0(date," ",expNum," ", "Normalized Counts.csv"))
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

## Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

## PCA ----------------------------------------------------
pca <- plotPCA(vsd, intgroup=c("tissue","donor", "subset"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

subsetfill <- c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[pca$subset]
donorshapes <- c("178"=21, "201"=22, "206"=23, "217"=24, "223"=25) [pca$donor]
tissuecolor <- c(Jejunum="#212121", Colon="#000000", Lung="#717171") [pca$tissue]

cairo_pdf(paste0(date," ",expNum," Tissue, Subset & Donor PCA.pdf"), w=6, h=4)
ggplot(pca, aes(PC1, PC2, shape=donor, color=tissue)) +
  geom_point(size=4, stroke = 0.5, aes(shape=donor, color=tissue, fill=subset)) +
  scale_shape_manual(values = donorshapes) +
  scale_color_manual(values = tissuecolor) +
  scale_fill_manual(values = subsetfill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste0(expNum," Tissue, Donor & Subset PCA.pdf")) 
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
  subset = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(vsd)$subset],
  tissue = c( Spleen="black", Lung="light gray", Jejunum="#707070" )[colData(vsd)$tissue])

## Sample distance heatmap ---------------------------------------------------------
sampleDists <- as.matrix(dist(t(assay(vsd))))

normheatmap(sampleDists,
            cex=0.5,
            title= "Sample Distance Heatmap")

## tSNE ----------------------------------------------------
rtsnesample <- Rtsne(sampleDists, is_distance = TRUE, perplexity = 5)

subsetfill <- c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(vsd)$subset]
donorshapes <- c("178"=21, "201"=22, "206"=23, "217"=24, "223"=25) [colData(vsd)$donor]
tissuecolor <- c(Jejunum="#212121", Colon="#000000", Lung="#717171") [colData(vsd)$tissue]

cairo_pdf(paste0(date," ",expNum," Tissue, Donor & Subset tSNE.pdf"), w=6, h=4)
plot(rtsnesample$Y, col=tissuecolor, bg=subsetfill, pch=donorshapes, cex=1, xlab="tSNE 1", ylab="tSNE 2", main="2-1B-1 Tissue, Donor & Subset tSNE")
dev.off()

## Top Variable Genes Heatmap ----------------------------------------------------
minNOTzero <- rownames(assay(vsd)[which(rowMin(counts(dds))>30),])
minNOTzero <- minNOTzero[which(rowVars(assay(vsd)[minNOTzero,])>1)]

topVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero,]), decreasing=TRUE)], 100)
write.csv(counts(dds, normalized=TRUE)[topVarGenes,colorder],file= paste0(date," ",expNum," ", "Top Variable Genes.csv"))

resheatmap(vsd, 
           genes = topVarGenes,
           samples= colnames(vsd),
           cex= 1,
           title= "Top Variable Genes Heatmap")
dev.off()

## Selected Genes Heatmaps ----------------------------------------------------

ilc2genes <- c("Gata3", "Il1rl1", "Il4", "Il5", "Il13", "Areg", "Ccl5", "Ccl1", "Ccr2", "Cxcl1", "Cd28", "Alox5", "Hspa1b", "Hspa1a", "Acbd7", "C5ar2")

resheatmap(vsd, 
           genes = ilc2genes,
           samples= -grep("MWAT",colnames(vsd)),
           h=4,
           w=5,
           filename="2-1A-4 ILC2 Genes Heatmap.pdf")
dev.off()

adiposegenes <- read.table("2-1A-4 Adipose Genes.txt", stringsAsFactors = FALSE)
adiposegenes <- as.character(adiposegenes[[1]])
uniquegenes <- unique(c(adiposegenes,ilc2genes))
badgenes <- c("Ppargc1a","Ppargc1b", "Pparg", "Il7r", "Ppargc1a", "Il21", "Il21r", "Xbp1", "Il33", "Batf")
uniquegenes <- uniquegenes[!uniquegenes %in% badgenes]

resheatmap(vsd, 
           genes = uniquegenes,
           samples= -grep("Lung",colnames(vsd)),
           h=5,
           w=4,
           filename="2-1A-4 Adipose ILC2 Genes Heatmap.pdf")
dev.off()

receptorgenes <- read.table("receptors.txt", stringsAsFactors = FALSE)
receptorgenes <- as.character(receptorgenes[[1]])
resheatmap(vsd, 
           genes = unique(receptorgenes),
           samples= colnames(vsd),
           h=8,
           w=4,
           filename="2-1A-4 Tissue Receptor Genes Heatmap.pdf")
dev.off()

signaturegenes <- unique(c(receptorgenes,uniquegenes))
signaturegenes <- unique(c(signaturegenes,"Il10ra"))

resheatmap(vsd, 
           genes = signaturegenes,
           samples= colnames(vsd),
           h=11,
           w=5,
           filename="2-1A-4 Tissue Signature Genes Heatmap.pdf")
dev.off()

functiongenes <-read.delim("genefunctions.txt", stringsAsFactors=FALSE)

lunggenes <- as.character(functiongenes[[1]])
lunggenes <- lunggenes[which(lunggenes %in% rownames(assay(vsd)))]
lunggenes <- lunggenes[which(rowMin(counts(dds)[lunggenes, ])>10, useNames=TRUE)]
lunggenes <- lunggenes[order(rowVars(assay(vsd)[lunggenes,]), decreasing=FALSE)]
resheatmap(vsd, 
           genes = lunggenes,
           samples= colnames(vsd),
           h=15,
           w=4,
           filename="2-1A-4 Lung Genes Heatmap.pdf")
dev.off()


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
sdg.SI.Colon <- subset(res.SI.Colon[order(res.SI.Colon$log2FoldChange),] , padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.SI.Lung <- subset(res.SI.Lung[order(res.SI.Lung$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.SI.MWAT <- subset(res.SI.MWAT[order(res.SI.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.Colon.Lung <- subset(res.Colon.Lung[order(res.Colon.Lung$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.Colon.MWAT <- subset(res.Colon.MWAT[order(res.Colon.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)
sdg.Lung.MWAT <- subset(res.Lung.MWAT[order(res.Lung.MWAT$log2FoldChange),], padj < 0.1 & abs(log2FoldChange)>2 & baseMean >10)

## Write DE results
write.csv(sdg.SI.Colon,file="2-1A-4 SI.Colon-DE.csv")
write.csv(sdg.SI.Lung,file="2-1A-4 SI.Lung-DE.csv")
write.csv(sdg.SI.MWAT,file="2-1A-4 SI.MWAT-DE.csv")
write.csv(sdg.Colon.Lung,file="2-1A-4 Colon.Lung-DE.csv")
write.csv(sdg.Colon.MWAT,file="2-1A-4 Colon.MWAT-DE.csv")
write.csv(sdg.Lung.MWAT,file="2-1A-4 Lung.MWAT-DE.csv")

## Tissue DE Genes Heatmap ---------------------------------------------------------------------
sigDEGenes <- unique(c(rownames(sdg.SI.Colon),
                           rownames(sdg.SI.Lung),
                           rownames(sdg.SI.MWAT),
                           rownames(sdg.Colon.Lung),
                           rownames(sdg.Colon.MWAT),
                           rownames(sdg.Lung.MWAT)))
sigDEGenes <- sigDEGenes[which(rowMeans(counts(dds)[sigDEGenes,])>10, useNames=TRUE)]
sigDEGenes <- sigDEGenes[order(rowVars(assay(vsd)[sigDEGenes, ]), decreasing=TRUE)]
write.csv(counts(dds, normalized=TRUE)[sigDEGenes,], file="2-1A-4 Tissue DE Genes.csv")

resheatmap(vsd, 
          genes = sigDEGenes[1:100], 
          samples=colnames(vsd),
          h=15,
          w=6,
          filename="2-1A-4 Tissue DE Genes Heatmap.pdf")
dev.off()

filteredGenes <- rownames(filteredDEgenes)
resheatmap(vsd, 
           genes = filteredGenes[1:100], 
           samples=colnames(vsd),
           h=15,
           w=6,
           filename="2-1A-4 Filtered DE Genes Heatmap.pdf")
dev.off()

## Tissue Functional Genes Heatmaps ---------------------------------------------------------------------
tissuegenes <- read.delim("~/Google Drive/Weill Cornell/Research/Data Analysis/RNAseq/2-1A-4 RNAseq/tissuegenes.txt", 
                          na.strings="",
                          stringsAsFactors=FALSE)

mwatGenes <- as.character(na.omit(tissuegenes[[1]]))
mwatGenes <- mwatGenes[which(mwatGenes %in% rownames(assay(vsd)))]
mwatGenes <- mwatGenes[which(rowMeans(counts(dds)[mwatGenes,-grep("Lung", colnames(vsd))])>20, useNames=TRUE)]
mwatGenes <- mwatGenes[which(rowMin(counts(dds)[mwatGenes,-grep("Lung", colnames(vsd))])>0, useNames=TRUE)]
mwatGenes <- mwatGenes[order(rowVars(assay(vsd)[mwatGenes, ]), decreasing=TRUE)]

colonGenes <- as.character(na.omit(tissuegenes[[2]]))
colonGenes <- colonGenes[which(colonGenes %in% rownames(assay(vsd)))]
colonGenes <- colonGenes[which(rowMeans(counts(dds)[colonGenes,-grep("Lung", colnames(vsd))])>20, useNames=TRUE)]
colonGenes <- colonGenes[which(rowMin(counts(dds)[colonGenes,-grep("Lung", colnames(vsd))])>0, useNames=TRUE)]
colonGenes <- colonGenes[order(rowVars(assay(vsd)[colonGenes, ]), decreasing=TRUE)]

siGenes <- as.character(na.omit(tissuegenes[[3]]))
siGenes <- siGenes[which(siGenes %in% rownames(assay(vsd)))]
siGenes <- siGenes[which(rowMeans(counts(dds)[siGenes,-grep("Lung", colnames(vsd))])>20, useNames=TRUE)]
siGenes <- siGenes[which(rowMin(counts(dds)[siGenes,-grep("Lung", colnames(vsd))])>0, useNames=TRUE)]
siGenes <- siGenes[order(rowVars(assay(vsd)[siGenes, ]), decreasing=TRUE)]  

lungGenes <- as.character(na.omit(tissuegenes[[4]]))
lungGenes <- lungGenes[which(lungGenes %in% rownames(assay(vsd)))]
lungGenes <- lungGenes[which(rowMeans(counts(dds)[lungGenes,-grep("Lung", colnames(vsd))])>20, useNames=TRUE)]
lungGenes <- lungGenes[which(rowMin(counts(dds)[lungGenes, -grep("Lung", colnames(vsd))])>0, useNames=TRUE)]
lungGenes <- lungGenes[order(rowVars(assay(vsd)[lungGenes, ]), decreasing=TRUE)]

tissueGenes <- unique (c(mwatGenes, colonGenes, siGenes, lungGenes))
tissueGenes <- tissueGenes[order(rowVars(assay(vsd)[tissueGenes, ]), decreasing=TRUE)]

resheatmap(vsd, 
           genes = tissueGenes [1:200],
           samples=colnames(vsd),
           h=30,
           w=6,
           filename="2-1A-4 Tissue Functional Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = mwatGenes,
           samples=-grep("D2", colnames(vsd)),
           h=6,
           w=6,
           filename="2-1A-4 MWAT Functional Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = colonGenes, 
           samples=-grep("D2_Lung|D3_MWAT",colnames(vsd)),
           h=10,
           w=6,
           filename="2-1A-4 Colon Functional Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = siGenes, 
           samples=-grep("D2_Lung|D3_MWAT",colnames(vsd)),
           h=35,
           w=6,
           filename="2-1A-4 SI Functional Genes Heatmap.pdf")
dev.off()

resheatmap(vsd, 
           genes = lungGenes, 
           samples=colnames(vsd),
           h=40,
           w=6,
           filename="2-1A-4 Lung Functional Genes Heatmap.pdf")
dev.off()

## Violin Plots --------------------- 
library(ggplot2)
library(gridExtra)
library(reshape2)

signaturecounts <- t(counts(dds, normalized=TRUE)[signaturegenes,])
signaturecounts <- cbind2(coldata, signaturecounts)
meltedSC <- melt(signaturecounts, variable.name = "gene", value.name = "count" )

cairo_pdf(filename="2-1A-4 Tissue Signature Gene Plots.pdf", w=20, h=15)
ggplot(data = meltedSC,  aes(x = tissue, y = log(count), fill = tissue)) +
  geom_violin(scale="width") +
  facet_wrap(~gene, nrow=6 ) +
  scale_fill_manual(values = c(Colon="#E7A626", SI="#4266F6", MWAT="#269040", Lung="#C7302A" )) +
  ggtitle("Tissue Signature Genes") +
  theme(axis.title.x = element_blank(),
        legend.position='none') +
  ylab("Log Normalized Counts")
dev.off()







