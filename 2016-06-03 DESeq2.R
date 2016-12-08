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
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/3-1A-3/3-1A-3 DESeq2")
date <- paste0(Sys.Date())
expNum <- paste0("3-1A-3")
colorder <- c(7,1,12,4,8,2,10,5,9,3,11,6)

## DESeq2 ----------
## Import data from featureCounts
countdata <- read.delim("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 DESeq2/3-1A-3 Raw Counts.txt", stringsAsFactors=FALSE)
countdata <- countdata[-(which(duplicated(countdata[,1])==TRUE)),]
row.names(countdata) <- unlist(countdata$Symbol)

## Convert to matrix 
countdata <- as.matrix(countdata[,-1])
countdata <- countdata[,colorder]
## Assign condition 
names <- colnames(countdata)
getdiet <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
diet <- unlist(lapply(names, getdiet))
diet <- gsub("C","Control",diet)
diet <- gsub("F","HiFat",diet)

gettissue <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
tissue <- unlist(lapply(names, gettissue))
tissue <- gsub("E","EWAT",tissue)
tissue <- gsub("M","MWAT",tissue)

getsubset <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
subset <- unlist(lapply(names, getsubset))

rm(names)

group <- paste0(diet," ",subset," ",tissue)

## Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

(coldata <- data.frame(row.names=colnames(countdata), diet, subset, tissue, group))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~subset + diet + tissue)

### Run DESeq normalization
dds <- DESeq(dds, parallel = TRUE)

write.csv(counts(dds, normalized=TRUE),file= paste0(date," ",expNum," ", "Normalized Counts.csv"))
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

### Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
rld <- rlogTransformation(dds)

### PCA ----------------------------------------------------
pca <- plotPCA(vsd, intgroup=c("tissue","subset", "diet"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

subsetfill <- c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6" ) [pca$subset]
tissueshapes <- c(EWAT=21, MWAT=22) [pca$tissue]
dietcolor <- c(Control="#212121", HiFat="#269040") [pca$diet]

cairo_pdf(paste0(date," ",expNum," Subset & Diet PCA.pdf"), w=6, h=4)
ggplot(pca, aes(PC1, PC2, shape=diet, color=subset)) +
  geom_point(size=4, stroke = 1, aes(shape=tissue, color=diet, fill=subset)) +
  scale_shape_manual(values = tissueshapes) +
  scale_color_manual(values = dietcolor) +
  scale_fill_manual(values = subsetfill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste0(expNum," Subset & Diet PCA")) 
dev.off()

pca2 <- prcomp(assay(dds))
write.csv(pca2$x, file = "3-1A-3 PCA Genes.csv")
write.csv(pca2$rotation, file = "3-1A-3 PCA Samples.csv")
plot(pca2$rotation[,c(3,1)], bg=subsetfill, pch=tissueshapes, col=dietcolor, cex=2)

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
  subset = c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6" )[colData(vsd)$subset],
  diet = c(HiFat="#269040", Control="#212121")[colData(vsd)$diet],
  tissue = c(EWAT="#E7A626", MWAT="#9C27B0")[colData(vsd)$tissue] )

### Sample distance heatmaps ---------------------------------------------------------
sampleDists <- as.matrix(dist(t(assay(vsd))))
EWATdists <- as.matrix(dist(t(assay(vsd)[,grep("_E_", colnames(vsd))])))
MWATdists <- as.matrix(dist(t(assay(vsd)[,grep("_M_", colnames(vsd))])))
HiFatdists <- as.matrix(dist(t(assay(vsd)[,grep("_F_", colnames(vsd))])))
Controldists <- as.matrix(dist(t(assay(vsd)[,grep("_C_", colnames(vsd))])))

normheatmap(sampleDists,
            title= "Sample Distance Heatmap")
normheatmap(EWATdists,
            title = "EWAT Sample Distance Heatmap")
normheatmap(MWATdists,
            title = "MWAT Sample Distance Heatmap")
normheatmap(HiFatdists,
            title = "HiFat Sample Distance Heatmap")
normheatmap(Controldists,
            title = "Control Sample Distance Heatmap")

### Top Variable Genes Heatmap ----------------------------------------------------
minNOTzero <- rownames(assay(rld)[which(rowMin(counts(dds))>50),])
minNOTzero <- minNOTzero[which(rowVars(assay(rld)[minNOTzero,])>0)]

topVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero,]), decreasing=TRUE)], 100)
ewatVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_E_",colnames(rld))]), decreasing=TRUE)], 100)
mwatVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_M_",colnames(rld))]), decreasing=TRUE)], 100)
hifatVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_F_",colnames(rld))]), decreasing=TRUE)], 100)
controlVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_C_",colnames(rld))]), decreasing=TRUE)], 100)
tregVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_Treg",colnames(rld))]), decreasing=TRUE)], 100)
ilc2VarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_ILC2",colnames(rld))]), decreasing=TRUE)], 100)
cd4VarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, grep("_CD4",colnames(rld))]), decreasing=TRUE)], 100)

ilc2tregVarGenes <- head( minNOTzero[order( rowVars( assay(rld)[minNOTzero, -grep("_CD4",colnames(vsd))]), decreasing=TRUE)], 100)

allVarGenes <- unique(c(topVarGenes,
                        ewatVarGenes,
                        mwatVarGenes,
                        hifatVarGenes,
                        controlVarGenes,
                        tregVarGenes,
                        ilc2VarGenes,
                        cd4VarGenes))

write.csv(counts(dds, normalized=TRUE)[allVarGenes,colorder],file= paste0(date," ",expNum," ", "Top Variable Genes.csv"))

resheatmap(rld, 
           genes = ilc2tregVarGenes,
           samples= -grep("_CD4|_C_E_2_ILC2|_M_",colnames(rld)),
           cex= 0.9,
           w= 0.67,
           title= "ILC2 Treg EWAT Variable Genes Heatmap"
)
dev.off()

resheatmap(rld, 
           genes = topVarGenes,
           samples= colnames(rld),
           cex= 1,
           title= "Top Variable Genes Heatmap"
           )
dev.off()

resheatmap(vsd, 
           genes = ewatVarGenes,
           samples= grep("_E_",colnames(vsd)),
           cex= 0.9,
           w= 0.67,
           title= "EWAT Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd, 
           genes = mwatVarGenes,
           samples= grep("_M_",colnames(vsd)),
           cex= 0.9,
           w= 0.67,
           title= "MWAT Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd, 
           genes = hifatVarGenes,
           samples= grep("_F_",colnames(vsd)),
           cex= 0.9,
           w= 0.67,
           title= "HiFat Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd, 
           genes = controlVarGenes,
           samples= grep("_C_",colnames(vsd)),
           cex= 0.9,
           w= 0.67,
           title= "Control Top Variable Genes Heatmap"
)
dev.off()

### Group Design DESeq Analysis ---------------------------------------------------------------------
dds$group <- factor(paste0(dds$diet,"_",dds$tissue,"_",dds$subset))
design(dds) <- ~ group
ddsGroup <- DESeq(dds)
vsdGroup <- varianceStabilizingTransformation(ddsGroup)
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

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

#results <- deresults(ddsGroup)

sdgnames <- function (sdg, dds=ddsGroup, vsd=vsdGroup, rowMin=30, rowVars=2,...){
  filtered <- unique(rownames(sdg))
  filtered <- filtered[order( rowVars( assay(vsd)[filtered,]), decreasing=TRUE)]
  subset(filtered,
         rowMin(counts(dds)[filtered,]) > rowMin &
         rowVars(counts(dds)[filtered,]) > rowVars)
}  

### Tissue Differential Expression ---------------------------------------------------------------------
res.EWAT_HiFat_CD4.MWAT_HiFat_CD4 <- results(ddsGroup, contrast=c("group","HiFat_EWAT_CD4","HiFat_MWAT_CD4"))
res.EWAT_HiFat_Treg.MWAT_HiFat_Treg <- results(ddsGroup, contrast=c("group","HiFat_EWAT_Treg","HiFat_MWAT_Treg"))
res.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2 <- results(ddsGroup, contrast=c("group","HiFat_EWAT_ILC2","HiFat_MWAT_ILC2"))
res.EWAT_Control_CD4.MWAT_Control_CD4 <- results(ddsGroup, contrast=c("group","Control_EWAT_CD4","Control_MWAT_CD4"))
res.EWAT_Control_Treg.MWAT_Control_Treg <- results(ddsGroup, contrast=c("group","Control_EWAT_Treg","Control_MWAT_Treg"))
res.EWAT_Control_ILC2.MWAT_Control_ILC2 <- results(ddsGroup, contrast=c("group","Control_EWAT_ILC2","Control_MWAT_ILC2"))

sdg.EWAT_HiFat_CD4.MWAT_HiFat_CD4 <- sdg (res.EWAT_HiFat_CD4.MWAT_HiFat_CD4) 
sdg.EWAT_HiFat_Treg.MWAT_HiFat_Treg <- sdg (res.EWAT_HiFat_Treg.MWAT_HiFat_Treg)
sdg.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2 <- sdg (res.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2)
sdg.EWAT_Control_CD4.MWAT_Control_CD4 <- sdg (res.EWAT_Control_CD4.MWAT_Control_CD4)
sdg.EWAT_Control_Treg.MWAT_Control_Treg <- sdg (res.EWAT_Control_Treg.MWAT_Control_Treg)
sdg.EWAT_Control_ILC2.MWAT_Control_ILC2 <- sdg (res.EWAT_Control_ILC2.MWAT_Control_ILC2)

sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4 <- sdgnames (sdg.EWAT_HiFat_CD4.MWAT_HiFat_CD4) 
sdgnames.EWAT_HiFat_Treg.MWAT_HiFat_Treg <- sdgnames (sdg.EWAT_HiFat_Treg.MWAT_HiFat_Treg)
sdgnames.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2 <- sdgnames (sdg.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2)
sdgnames.EWAT_Control_CD4.MWAT_Control_CD4 <- sdgnames (sdg.EWAT_Control_CD4.MWAT_Control_CD4)
sdgnames.EWAT_Control_Treg.MWAT_Control_Treg <- sdgnames (sdg.EWAT_Control_Treg.MWAT_Control_Treg)
sdgnames.EWAT_Control_ILC2.MWAT_Control_ILC2 <- sdgnames (sdg.EWAT_Control_ILC2.MWAT_Control_ILC2)

tissueSDG <- unique(c(sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4,
                   sdgnames.EWAT_HiFat_Treg.MWAT_HiFat_Treg,
                   sdgnames.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2,
                   sdgnames.EWAT_Control_CD4.MWAT_Control_CD4,
                   sdgnames.EWAT_Control_Treg.MWAT_Control_Treg,
                   sdgnames.EWAT_Control_ILC2.MWAT_Control_ILC2))

write.csv(counts(ddsGroup, normalized=TRUE)[tissueSDG,colorder],file=paste0(date," ",expNum," Tissue Differentially Expressed Genes.csv" ))
tissueSDGlist <- list(sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4,
                      sdgnames.EWAT_HiFat_Treg.MWAT_HiFat_Treg,
                      sdgnames.EWAT_HiFat_ILC2.MWAT_HiFat_ILC2,
                      sdgnames.EWAT_Control_CD4.MWAT_Control_CD4,
                      sdgnames.EWAT_Control_Treg.MWAT_Control_Treg,
                      sdgnames.EWAT_Control_ILC2.MWAT_Control_ILC2)

CD4tissueSDG <- c(tissueSDGlist[[1]], tissueSDGlist[[4]])
CD4tissueSDG  <- CD4tissueSDG [order( rowMeans( assay(vsd)[CD4tissueSDG,]), decreasing=TRUE)]
TregtissueSDG <- c(tissueSDGlist[[2]], tissueSDGlist[[5]])
TregtissueSDG  <- TregtissueSDG [order( rowMeans( assay(vsd)[TregtissueSDG,]), decreasing=TRUE)]
ILC2tissueSDG <- c(tissueSDGlist[[3]], tissueSDGlist[[6]])
ILC2tissueSDG  <- ILC2tissueSDG[order( rowMeans( assay(vsd)[ILC2tissueSDG,]), decreasing=TRUE)]

resheatmap(vsdGroup, 
           genes = CD4tissueSDG[1:100],
           samples= grep("_CD4",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "CD4 Tissue SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = TregtissueSDG[1:100],
           samples= grep("_Treg",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "Treg Tissue SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = ILC2tissueSDG[1:100],
           samples= grep("_ILC2",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "ILC2 Tissue SDE Genes Heatmap")
dev.off()

### Diet Differential Expression ---------------------------------------------------------------------
res.EWAT_HiFat_CD4.EWAT_Control_CD4 <- results(ddsGroup, contrast=c("group","HiFat_EWAT_CD4","Control_EWAT_CD4"))
res.EWAT_HiFat_Treg.EWAT_Control_Treg <- results(ddsGroup, contrast=c("group","HiFat_EWAT_Treg","Control_EWAT_Treg"))
res.EWAT_HiFat_ILC2.EWAT_Control_ILC2 <- results(ddsGroup, contrast=c("group","HiFat_EWAT_ILC2","Control_EWAT_ILC2"))
res.MWAT_HiFat_CD4.MWAT_Control_CD4 <- results(ddsGroup, contrast=c("group","HiFat_MWAT_CD4","Control_MWAT_CD4"))
res.MWAT_HiFat_Treg.MWAT_Control_Treg <- results(ddsGroup, contrast=c("group","HiFat_MWAT_Treg","Control_MWAT_Treg"))
res.MWAT_HiFat_ILC2.MWAT_Control_ILC2 <- results(ddsGroup, contrast=c("group","HiFat_MWAT_ILC2","Control_MWAT_ILC2"))

sdg.EWAT_HiFat_CD4.EWAT_Control_CD4 <- sdg (res.EWAT_HiFat_CD4.EWAT_Control_CD4) 
sdg.EWAT_HiFat_Treg.EWAT_Control_Treg <- sdg (res.EWAT_HiFat_Treg.EWAT_Control_Treg)
sdg.EWAT_HiFat_ILC2.EWAT_Control_ILC2 <- sdg (res.EWAT_HiFat_ILC2.EWAT_Control_ILC2)
sdg.MWAT_HiFat_CD4.MWAT_Control_CD4 <- sdg (res.MWAT_HiFat_CD4.MWAT_Control_CD4) 
sdg.MWAT_HiFat_Treg.MWAT_Control_Treg <- sdg (res.MWAT_HiFat_Treg.MWAT_Control_Treg)
sdg.MWAT_HiFat_ILC2.MWAT_Control_ILC2 <- sdg (res.MWAT_HiFat_ILC2.MWAT_Control_ILC2)

sdgnames.EWAT_HiFat_CD4.EWAT_Control_CD4 <- sdgnames (sdg.EWAT_HiFat_CD4.EWAT_Control_CD4) 
sdgnames.EWAT_HiFat_Treg.EWAT_Control_Treg <- sdgnames (sdg.EWAT_HiFat_Treg.EWAT_Control_Treg)
sdgnames.EWAT_HiFat_ILC2.EWAT_Control_ILC2 <- sdgnames (sdg.EWAT_HiFat_ILC2.EWAT_Control_ILC2)
sdgnames.MWAT_HiFat_CD4.MWAT_Control_CD4 <- sdgnames (sdg.MWAT_HiFat_CD4.MWAT_Control_CD4) 
sdgnames.MWAT_HiFat_Treg.MWAT_Control_Treg <- sdgnames (sdg.MWAT_HiFat_Treg.MWAT_Control_Treg)
sdgnames.MWAT_HiFat_ILC2.MWAT_Control_ILC2 <- sdgnames (sdg.MWAT_HiFat_ILC2.MWAT_Control_ILC2)

 
dietSDG <- unique(c(sdgnames.EWAT_HiFat_CD4.EWAT_Control_CD4,
                sdgnames.EWAT_HiFat_Treg.EWAT_Control_Treg,
                sdgnames.EWAT_HiFat_ILC2.EWAT_Control_ILC2,
                sdgnames.MWAT_HiFat_CD4.MWAT_Control_CD4,
                sdgnames.MWAT_HiFat_Treg.MWAT_Control_Treg,
                sdgnames.MWAT_HiFat_ILC2.MWAT_Control_ILC2))
write.csv(counts(ddsGroup, normalized=TRUE)[dietSDG,colorder],file=paste0(date," ",expNum," Diet Differentially Expressed Genes.csv" ))

dietSDGlist <- list(sdgnames.EWAT_HiFat_CD4.EWAT_Control_CD4,
                    sdgnames.EWAT_HiFat_Treg.EWAT_Control_Treg,
                    sdgnames.EWAT_HiFat_ILC2.EWAT_Control_ILC2,
                    sdgnames.MWAT_HiFat_CD4.MWAT_Control_CD4,
                    sdgnames.MWAT_HiFat_Treg.MWAT_Control_Treg,
                    sdgnames.MWAT_HiFat_ILC2.MWAT_Control_ILC2)

CD4dietSDG <- c(dietSDGlist[[1]], dietSDGlist[[4]])
CD4dietSDG  <- CD4dietSDG [order( rowMeans( assay(vsd)[CD4dietSDG,]), decreasing=TRUE)]
TregdietSDG <- c(dietSDGlist[[2]], dietSDGlist[[5]])
TregdietSDG  <- TregdietSDG [order( rowMeans( assay(vsd)[TregdietSDG,]), decreasing=TRUE)]
ILC2dietSDG <- c(dietSDGlist[[3]], dietSDGlist[[6]])
ILC2dietSDG  <- ILC2dietSDG[order( rowMeans( assay(vsd)[ILC2dietSDG,]), decreasing=TRUE)]

resheatmap(vsdGroup, 
           genes = CD4dietSDG[1:100],
           samples= grep("_CD4",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "CD4 Diet SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = TregdietSDG[1:100],
           samples= grep("_Treg",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "Treg Diet SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = ILC2dietSDG[1:100],
           samples= grep("_ILC2",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "ILC2 Diet SDE Genes Heatmap")
dev.off()

## Diet & Tissue SDG---------------------------
CD4sdg <- unique(c(CD4dietSDG, CD4tissueSDG))
CD4sdg <- CD4sdg [order( rowMeans( assay(vsd)[CD4sdg,]), decreasing=TRUE)]
Tregsdg <- unique(c(TregdietSDG, TregtissueSDG))
Tregsdg <- Tregsdg[order( rowMeans( assay(vsd)[Tregsdg,]), decreasing=TRUE)]
ILC2sdg <- unique(c(ILC2dietSDG, ILC2tissueSDG))
ILC2sdg <- ILC2sdg[order( rowMeans( assay(vsd)[ILC2sdg,]), decreasing=TRUE)]

resheatmap(vsdGroup, 
           genes = CD4sdg[1:100],
           samples= grep("_CD4",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.67,
           title= "CD4 Diet & Tissue SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = Tregsdg[1:100],
           samples= grep("_Treg",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.67,
           title= "Treg Diet & Tissue SDE Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = ILC2sdg[1:100],
           samples= grep("_ILC2",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.68,
           title= "ILC2 Diet & Tissue SDE Genes Heatmap")
dev.off()

diettissueSDG <- unique(c(dietSDG,tissueSDG))
diettissueSDG <- diettissueSDG[order( rowVars( assay(dds)[diettissueSDG,]), decreasing=TRUE)]
write.csv(counts(ddsGroup, normalized=TRUE)[diettissueSDG,colorder],file=paste0(date," ",expNum," Diet & Tissue Differentially Expressed Genes.csv" ))

## DEG STRING Analysis ---------------------------

stringGenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 STRING/stringGenes.txt", 
                          quote="\"", comment.char="", stringsAsFactors=FALSE)
stringGenes <- stringGenes[[1]]
stringGenes <- stringGenes[which(stringGenes %in% rownames(counts(dds)))]
stringGenes <- stringGenes[order( rowVars( assay(vsd)[stringGenes,]), decreasing=TRUE)]

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= colnames(vsdGroup),
           cex= 0.9,
           w= 0.8,
           title= "STRING Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= grep("_CD4",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "CD4 STRING Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= grep("_Treg",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "Treg STRING Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= grep("_ILC2",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "ILC2 STRING Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= grep("_E_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "EWAT STRING Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = stringGenes[1:100],
           samples= grep("_M_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "MWAT STRING Genes Heatmap")
dev.off()

## DEG DAVID Analysis ---------------------------
DAVIDgenes <- read.table("~/Desktop/DAVIDgenes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

DAVIDgenes <- DAVIDgenes[[1]]
DAVIDgenes <- DAVIDgenes[which(DAVIDgenes %in% rownames(assay(vsd)))]
DAVIDgenes <- DAVIDgenes[which(rowVars( assay(vsd)[DAVIDgenes,])>1.5 & rowMin(counts(dds)[DAVIDgenes,])>50)]
DAVIDgenes <- DAVIDgenes[order( rowVars( assay(vsd)[DAVIDgenes,]), decreasing=TRUE)]

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= -grep("CD4|_E_", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "Treg ILC2 MWAT DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= -grep("CD4|_M_", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "Treg ILC2 EWAT DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_CD4|",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "CD4 DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_Treg",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "Treg DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_ILC2",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "ILC2 DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_E_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "EWAT DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_M_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "MWAT DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_F_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "HiFat DAVID Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = DAVIDgenes,
           samples= grep("_C_",colnames(vsdGroup)),
           cex= 0.9,
           w= 0.6,
           title= "Control DAVID Genes Heatmap")
dev.off()

## DEG Signature Genes ----------
signaturegenes <- read.table("~/Desktop/3-1A-3 signaturegenes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
signaturegenes <- signaturegenes[[1]]

signaturegenes <- signaturegenes[which(signaturegenes %in% rownames(assay(vsd)))]
signaturegenes <- signaturegenes[which(rowVars( assay(vsd)[signaturegenes,])>0 & rowMin(counts(dds)[signaturegenes,])>50)]
signaturegenes <- signaturegenes[order( rowVars( assay(vsd)[signaturegenes,]), decreasing=TRUE)]
signaturegenes <- unique(signaturegenes)

write.csv(counts(dds, normalized=TRUE)[signaturegenes,colorder],file= paste0(date," ",expNum," ", "Signature Genes.csv"))

resheatmap(vsd, 
           genes = signaturegenes,
           samples= colnames(vsdGroup),
           cex= 0.9,
           w= 0.6,
           title= "Signature Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = signaturegenes,
           samples= -grep("_E_", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "MWAT Signature Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = signaturegenes,
           samples= -grep("_M_", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "EWAT Signature Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = signaturegenes,
           samples= grep("_CD4", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "CD4 Signature Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = signaturegenes,
           samples= grep("_ILC2", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "ILC2 Signature Genes Heatmap")
dev.off()

resheatmap(vsdGroup, 
           genes = signaturegenes,
           samples= grep("_Treg", colnames(vsdGroup)),
           cex= 0.8,
           w= 0.5,
           title= "Treg Signature Genes Heatmap")
dev.off()

## Signature Gene Bar Plots --------------------- 
library(ggplot2)
library(gridExtra)
library(reshape2)

signaturecounts <- t(counts(dds, normalized=TRUE)[signaturegenes,])
signaturecounts <- cbind2(coldata, signaturecounts)
meltedSC <- melt(signaturecounts, variable.name = "gene", value.name = "count" )

subsetfill <- c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6" ) [meltedSC$subset]
tissueshapes <- c(EWAT=21, MWAT=22) [meltedSC$tissue]
dietcolor <- c(Control="#212121", HiFat="#269040") [meltedSC$diet]

cairo_pdf(paste0(date," ",expNum," Signature Genes.pdf"), w=15, h=10)
ggplot(data=meltedSC, aes(x=subset, y=log2(count), color=diet, fill=subset)) +
  facet_wrap(tissue ~ gene, scales="free_y", ncol = 10) +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.75) +
  scale_fill_manual(values = c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6")) +
  scale_color_manual(values = c(Control="#212121", HiFat="#269040")) +
  labs(title=paste0(expNum," Signature Genes"), x="Subset", y= "Log (counts)") +
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.background = element_rect(color="#E5E5E5", linetype="solid", fill=NA),
        panel.grid.major = element_line(color="#E5E5E5", linetype="dotted"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = rel(2)))
dev.off()



mTGFBgenes <- read.table("~/Desktop/mTGFBgenes.txt", quote="\"", comment.char="", na.strings="", stringsAsFactors=FALSE)
mTGFBgenes <- mTGFBgenes[[1]]
mTGFBgenes <- mTGFBgenes[which(mTGFBgenes %in% rownames(counts(dds)))]
mTGFBgenes <- mTGFBgenes[which(rowMin(counts(dds)[mTGFBgenes,])>0)]


mTGFBgenecounts <- t(counts(dds, normalized=TRUE)[mTGFBgenes,])
mTGFBgenecounts <- cbind2(coldata, mTGFBgenecounts)
meltedTGFB <- melt(mTGFBgenecounts, variable.name = "gene", value.name = "count" )

subsetfill <- c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6" ) [meltedTGFB$subset]
tissueshapes <- c(EWAT=21, MWAT=22) [meltedTGFB$tissue]
dietcolor <- c(Control="#212121", HiFat="#269040") [meltedTGFB$diet]

cairo_pdf(paste0(date," ",expNum," TGFb Genes.pdf"), w=5, h=30)
ggplot(data=meltedTGFB, aes(x=tissue, y=log2(count), color=diet, fill=subset)) +
  facet_wrap(  ~ gene+subset, scales="free_y", ncol = 6) +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.75) +
  scale_fill_manual(values = c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6")) +
  scale_color_manual(values = c(Control="#212121", HiFat="#269040")) +
  labs(title=paste0(expNum," TGFb Genes"), x="Subset", y= "Log (counts)") +
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.background = element_rect(color="#E5E5E5", linetype="solid", fill=NA),
        panel.grid.major = element_line(color="#E5E5E5", linetype="dotted"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = rel(2)))
dev.off()

