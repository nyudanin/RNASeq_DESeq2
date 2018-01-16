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
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/2-2B-1/2-2B-1 DESeq2")
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
colorder <- c(grep("NK", colnames(dds)), grep("ILC1", colnames(dds)), grep("ILC2", colnames(dds)), grep("ILC3", colnames(dds)))

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

pcaTissueDonor <- plotPCA(vsd, intgroup=c("tissue","donor"), returnData=TRUE)
percentVarTissueDonor <- round(100 * attr(pcaTissueDonor, "percentVar"))

pcaSubsetDonor <- plotPCA(vsd, intgroup=c("subset","donor"), returnData=TRUE)
percentVarSubsetDonor <- round(100 * attr(pcaSubsetDonor, "percentVar"))

pcaSubsetTissue <- plotPCA(vsd, intgroup=c("subset","tissue"), returnData=TRUE)
percentVarSubsetTissue <- round(100 * attr(pcaSubsetTissue, "percentVar"))

subsetfill <- c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[pca$subset]
donorshapes <- c("178"=21, "201"=22, "216"=23, "217"=24, "223"=25) [pca$donor]
tissuecolor <- c(Jejunum="#717171", Spleen="#000000", Lung="#FFFFFF") [pca$tissue]

cairo_pdf(paste0(date," ",expNum," Tissue & Donor PCA.pdf"), w=8, h=6)
ggplot(pcaTissueDonor, aes(PC1, PC2, shape=donor, color=tissue)) +
  geom_point(size=3, stroke = 0.8, aes(shape=donor, color=tissue, fill=tissue)) +
  scale_shape_manual(values = donorshapes) +
  scale_color_manual(values = tissuecolor) +
  scale_fill_manual(values = tissuecolor) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVarTissueDonor[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarTissueDonor[2],"% variance")) +
  ggtitle(paste0(expNum," Tissue & Donor PCA.pdf")) 
dev.off()

cairo_pdf(paste0(date," ",expNum," Subset & Donor PCA.pdf"), w=8, h=6)
ggplot(pcaTissueDonor, aes(PC1, PC2, shape=donor, color=tissue)) +
  geom_point(size=3, stroke = 0.8, aes(shape=donor, color=tissue, fill=subset)) +
  scale_shape_manual(values = donorshapes) +
  scale_color_manual(values = tissuecolor) +
  scale_fill_manual(values = subsetfill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVarSubsetDonor[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarSubsetDonor[2],"% variance")) +
  ggtitle(paste0(expNum," Subset & Donor PCA.pdf")) 
dev.off()

cairo_pdf(paste0(date," ",expNum," Tissue, Subset & Donor PCA.pdf"), w=8, h=6)
ggplot(pca, aes(PC1, PC2, shape=donor, color=subset)) +
  geom_point(size=3, stroke = 1, aes(shape=donor, color=subset, fill=tissue)) +
  scale_shape_manual(values = donorshapes) +
  scale_color_manual(values = subsetfill) +
  scale_fill_manual(values = tissuecolor) +
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
                annotation_col=coldata[c(1,3)],
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
donorshapes <- c("178"=21, "201"=22, "216"=23, "217"=24, "223"=25) [colData(vsd)$donor]
tissuecolor <- c(Jejunum="#212121", Spleen="#000000", Lung="#717171") [colData(vsd)$tissue]

cairo_pdf(paste0(date," ",expNum," Tissue, Donor & Subset tSNE.pdf"), w=6, h=4)
plot(rtsnesample$Y, col=tissuecolor, bg=subsetfill, pch=donorshapes, cex=1, xlab="tSNE 1", ylab="tSNE 2", main="2-1B-1 Tissue, Donor & Subset tSNE")
dev.off()

## Top Variable Genes Heatmap ----------------------------------------------------
minNOTzero <- rownames(assay(vsd)[which(rowMin(counts(dds))>1),])
minNOTzero <- minNOTzero[which(rowVars(assay(vsd)[minNOTzero,])>1)]

topVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero,]), decreasing=TRUE)], 100)
ilc1VarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_ILC1",colnames(vsd))]), decreasing=TRUE)], 100)
ilc2VarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_ILC2",colnames(vsd))]), decreasing=TRUE)], 100)
ilc3VarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_ILC3",colnames(vsd))]), decreasing=TRUE)], 100)
nkVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_NK",colnames(vsd))]), decreasing=TRUE)], 100)

allVarGenes <- unique(c(topVarGenes,
                        ilc3VarGenes,
                        ilc2VarGenes,
                        ilc1VarGenes,
                        nkVarGenes))

write.csv(counts(dds, normalized=TRUE)[topVarGenes,colorder],file= paste0(date," ",expNum," ", "Top Variable Genes.csv"))
write.csv(counts(dds, normalized=TRUE)[allVarGenes,colorder],file= paste0(date," ",expNum," ", "All Variable Genes.csv"))

resheatmap(vsd, 
           genes = topVarGenes,
           samples= colnames(vsd),
           cex= 0.7,
           w=1.2,
           title= "Top Variable Genes Heatmap")
dev.off()

resheatmap(vsd, 
           genes = allVarGenes,
           samples= colnames(vsd),
           cex= 0.7,
           w=1.2,
           title= "All Variable Genes Heatmap")
dev.off()

remove <- read.table("~/Desktop/Remove.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
remove <- remove[[1]]
remove <- remove[which(remove %in% rownames(counts(dds)))]
newVarGenes <- unique(c(allVarGenes,remove))

resheatmap(vsd, 
           genes = newVarGenes,
           samples= colnames(vsd),
           cex= 0.7,
           w=1.2,
           title= "Top Variable Genes Heatmap")
dev.off()

### Differential Expression Functions---------------------------------------------------------------------
deresults <- function (dds, fcThresh=2, baseMeanThresh=100, padjThresh=0.1){
  reslist <- list()
  resnames <- resultsNames(dds)[-1]
  for (num in 2:length(resnames)){
    for (den in 1:(num-1)){
      print(paste0(resnames[num],".", resnames[den]))
      res <- results(dds, contrast=list(resnames[num], resnames[den]))
      res <- subset(res, abs(log2FoldChange) > fcThresh  & baseMean > baseMeanThresh & padj < padjThresh)
      reslist[[paste0(resnames[num],".", resnames[den])]] <- res
    }
  }
  return(reslist)
}

### Differential Expression---------------------------------------------------------------------
sdg <- deresults(dds)

## Comparison to 3-1A-3 --------------------- 
signaturegenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 DAVID/3-1A-3 signaturegenes.txt", 
                             quote="\"", comment.char="", na.strings="", stringsAsFactors=FALSE)
signaturegenes <- toupper(signaturegenes[[1]])
signaturegenes <- signaturegenes[which(signaturegenes %in% rownames(counts(dds)))]
signaturegenes <- signaturegenes[order( rowVars( assay(vsd)[signaturegenes,]), decreasing=TRUE)]
signaturegenes <- subset(signaturegenes, 
                         rowMin(counts(dds)[signaturegenes,])>0)
signaturegenes <- signaturegenes[order(rowVars(assay(vsd)[signaturegenes, ]), decreasing=TRUE)]

resheatmap(vsd, 
           genes = signaturegenes,
           samples= colnames(vsd),
           cex= 0.9,
           h=0.7,
           w=1.2,
           title= "3-1A-3 Comparison Heatmap")
dev.off()

## Comparison to 2-1A-4 --------------------- 
tissuesignaturegenes <- read.csv("/Volumes/IBD/Yudanin/RNAseq/2-1A-4 RNAseq/2-1A-4 DESeq2/2-1A-4 Tissue Signature Genes.csv",
                                 header=FALSE, stringsAsFactors=FALSE)
tissuesignaturegenes <- toupper(tissuesignaturegenes[[1]])
tissuesignaturegenes <- tissuesignaturegenes[which(tissuesignaturegenes %in% rownames(counts(dds)))]
tissuesignaturegenes <- subset(tissuesignaturegenes, 
                         rowVars(counts(dds)[tissuesignaturegenes,])>1)
tissuesignaturegenes <- tissuesignaturegenes[order( rowVars( assay(vsd)[tissuesignaturegenes,]), decreasing=TRUE)]

resheatmap(vsd, 
           genes = tissuesignaturegenes,
           samples= colnames(vsd),
           cex= 0.9,
           h=0.7,
           w=1.2,
           title= "2-1A-4 Comparison Heatmap")
dev.off()


## Tissue & Subset Signature Genes --------------------- 
tissueGenes <- read.csv("/Volumes/IBD/Yudanin/RNAseq/2-2B-1 RNAseq/2-2B-1 DESeq2/All Samples/2-2B-1 Tissue Signature Genes.csv",
                        header=FALSE, stringsAsFactors=FALSE)
tissueGenes <- toupper(tissueGenes[[1]])
tissueGenes <- tissueGenes[which(tissueGenes %in% rownames(counts(dds)))]
tissueGenes <- subset(tissueGenes, 
                               rowVars(counts(dds)[tissueGenes,])>1)
tissueGenes <- tissueGenes[order( rowVars( assay(vsd)[tissueGenes,]), decreasing=TRUE)]

subsetGenes <- read.csv("/Volumes/IBD/Yudanin/RNAseq/2-2B-1 RNAseq/2-2B-1 DESeq2/All Samples/2-2B-1 Subset Signature Genes.csv",
                        header=FALSE, stringsAsFactors=FALSE)
subsetGenes <- toupper(subsetGenes[[1]])
subsetGenes <- subsetGenes[which(subsetGenes %in% rownames(counts(dds)))]
subsetGenes <- subset(subsetGenes, 
                      rowMin(counts(dds)[subsetGenes,])>0)
subsetGenes <- subsetGenes[order( rowVars( assay(vsd)[subsetGenes,]), decreasing=TRUE)]

subsettissueGenes <- unique(c(tissueGenes, subsetGenes))
subsettissueGenes <- subsettissueGenes[order( rowVars( assay(vsd)[subsettissueGenes,]), decreasing=TRUE)]
resheatmap(vsd, 
           genes = subsettissueGenes,
           samples= colnames(vsd),
           cex= 0.8,
           h=0.6,
           w=1.2,
           title= "Tissue & Subset Signature Genes Heatmap")
dev.off()

ilc1NKgenes <- unlist(ilc1NKgenes)
ilc1NKgenes <- read.table("~/Desktop/GeneList.txt", quote="\"", comment.char="")
ilc1NKgenes <- toupper(ilc1NKgenes[1])
ilc1NKgenes <- ilc1NKgenes[which(ilc1NKgenes %in% rownames(counts(dds)))]
ilc1NKgenes <- subset(ilc1NKgenes, 
                      rowMin(counts(dds)[ilc1NKgenes,])>0)
ilc1NKgenes <- ilc1NKgenes[order( rowVars( assay(vsd)[ilc1NKgenes,]), decreasing=TRUE)]

resheatmap(dds, 
           genes = ilc1NKgenes,
           samples= colnames(dds),
           cex= 0.8,
           h=0.6,
           w=1.2,
           title= "ILC1 vs NK Genes Heatmap")
dev.off()

ilc1NKgenes <- as.character(ilc1NKgenes)
char


write.csv(counts(dds, normalized=TRUE)[ilc1NKgenes,],file= paste0(date," ",expNum," ", "ILC1 vs NK Genes.csv"))
