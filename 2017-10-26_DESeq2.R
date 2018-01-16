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
dir <- setwd("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/3-1E-1/3-1E-1 DESeq2")
date <- paste0(Sys.Date())
expNum <- paste0("3-1E-1")


## DESeq2 ----------
## Import data from featureCounts
countdata <- read.delim("~/Google Drive/Weill Cornell/Research/Analysis/RNAseq/3-1E-1/3-1E-1 DESeq2/3-1E-1 Raw Counts.txt", row.names=1)
countdata <- countdata[-(which(duplicated(countdata[,1])==TRUE)),]
row.names(countdata) <- unlist(countdata$Symbol)

## Convert to matrix
countdata <- as.matrix(countdata[,-(1)])
colorder <- c(grep("ILC",colnames(countdata)), grep("Treg",colnames(countdata)))
countdata <- countdata[,colorder]

## Assign condition
names <- colnames(countdata)
getdiet1 <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
diet1 <- unlist(lapply(names, getdiet1))
diet1 <- gsub("C","Control",diet1)
diet1 <- gsub("F","HiFat",diet1)

getdiet2 <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
diet2 <- unlist(lapply(names, getdiet2))
diet2 <- gsub("C","Control",diet2)
diet2 <- gsub("F","HiFat",diet2)

getreplicate <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
replicate <- unlist(lapply(names, getreplicate))
replicate <- gsub( "1", "R1", replicate)
replicate <- gsub( "2", "R2", replicate)
replicate <- gsub( "3", "R3", replicate)

getsubset <- function (colname){
  strsplit(colname,"_")[[1]][4]
}
subset <- unlist(lapply(names, getsubset))

rm(names)

dietgroup <- paste0(diet1," ",diet2)
samplegroup <- paste0(diet1," ",diet2," ",subset)

## Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

(coldata <- data.frame(row.names=colnames(countdata), diet1, diet2, subset, dietgroup, samplegroup))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~dietgroup + subset)

### Run DESeq normalization
dds <- DESeq(dds, parallel = TRUE)

write.csv(counts(dds, normalized=TRUE),file= paste0(date," ",expNum," ", "Normalized Counts.csv"))
save.image(paste0(dir,"/", date," ",expNum, " DESeq2.RData"))

### Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)
rld <- rlogTransformation(dds)

### PCA ----------------------------------------------------
pca <- plotPCA(vsd, intgroup=c("diet1", "subset"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

dietgroupfill <- c("Control Control" = "#0057E7", "Control HiFat" = "#BBDEFB", "HiFat Control" = "#C8E6C9", "HiFat HiFat" = "#008744") [pca$diet1]
subsetshapes <- c(Treg=21, ILC2=22) [pca$subset]
subsetcolors <- c(Treg="#D62D20", ILC2="#212121") [pca$subset]

pdf(paste0(date," ",expNum," Subset & Diet PCA.pdf"), w=6, h=6)
ggplot(pca, aes(PC1, PC2, shape=subset, color=diet1)) +
  geom_point(size=4, stroke = 1, aes(shape=subset, color=subset, fill=diet1)) +
  scale_shape_manual(values = subsetshapes) +
  scale_color_manual(values = subsetcolors) +
  scale_fill_manual(values = diet1fill) +
  theme(legend.title = element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste0(expNum," Subset & Diet PCA"))
dev.off()

pca2 <- prcomp(assay(dds))
write.csv(pca2$x, file = "3-1E-1 PCA Genes.csv")
write.csv(pca2$rotation, file = "3-1E-1 PCA Samples.csv")

pdf(paste0(date," ",expNum," Subset & Diet PCA2.pdf"), w=7, h=7)
plot(pca2$rotation[,c(1,2)], bg=dietgroupfill, pch=subsetshapes, col=subsetcolors, cex=2)
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
  subset = c(Treg="#C7302A", ILC2="#707070")[colData(vsd)$subset],
  dietgroup = c("Control Control" = "#0057E7", "Control HiFat" = "#BBDEFB", "HiFat Control" = "#C8E6C9", "HiFat HiFat" = "#008744")[colData(vsd)$dietgroup])

### Sample distance heatmaps ---------------------------------------------------------
sampleDists <- as.matrix(dist(t(assay(vsd))))
ILC2dists <- as.matrix(dist(t(assay(vsd)[,grep("_ILC", colnames(vsd))])))
Tregdists <- as.matrix(dist(t(assay(vsd)[,grep("_Treg", colnames(vsd))])))


normheatmap(sampleDists,
            title= "Sample Distance Heatmap")
normheatmap(ILC2dists,
            title = "ILC2 Sample Distance Heatmap")
normheatmap(Tregdists,
            title = "Treg Sample Distance Heatmap")



### Top Variable Genes Heatmap ----------------------------------------------------
minNOTzero <- rownames(assay(vsd)[which(rowMin(counts(dds))>100),-grep("_3",colnames(vsd))])
minNOTzero <- minNOTzero[which(rowVars(assay(vsd)[minNOTzero,])>1)]

topVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero,]), decreasing=TRUE)], 500)
tregVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_Treg",colnames(vsd))]), decreasing=TRUE)], 500)
ilc2VarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("_ILC2",colnames(vsd))]), decreasing=TRUE)], 500)

FatVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("F_",colnames(vsd))]), decreasing=TRUE)], 500)
ConVarGenes <- head( minNOTzero[order( rowVars( assay(vsd)[minNOTzero, grep("C_",colnames(vsd))]), decreasing=TRUE)], 500)

allVarGenes <- unique(c(topVarGenes,
                        tregVarGenes,
                        ilc2VarGenes,
                        FatVarGenes,
                        ConVarGenes))

write.csv(counts(dds, normalized=TRUE)[allVarGenes,colorder],file= paste0(date," ",expNum," ", "Top Variable Genes.csv"))

resheatmap(vsd,
           genes = allVarGenes[1:250],
           samples= -grep("_3", colnames(vsd)),
           cex= 0.5,
           w= 0.67,
           title= "Combined Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd,
           genes = topVarGenes,
           samples= colnames(vsd),
           cex= 1,
           title= "Top Variable Genes Heatmap"
           )
dev.off()

resheatmap(vsd,
           genes = topVarGenes,
           samples= -grep("C_F_2|F_F_3|F_C_3|_Treg|C_C_1", colnames(vsd)),
           cex= 0.4,
           w= 0.6,
           title= "ILC2 Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd,
           genes = tregVarGenes,
           samples= colnames(vsd),
           cex= 1,
           title= "Treg Top Variable Genes Heatmap"
)
dev.off()
resheatmap(vsd,
           genes = FatVarGenes,
           samples= colnames(vsd),
           cex= 0.5,
           w= 0.67,
           title= "Fat Top Variable Genes Heatmap"
)
dev.off()

resheatmap(vsd,
           genes = ConVarGenes,
           samples= colnames(vsd),
           cex= 0.5,
           w= 0.67,
           title= "Control Top Variable Genes Heatmap"
)
dev.off()
### Group Design DESeq Analysis ---------------------------------------------------------------------
dds$group <- factor(paste0(dds$diet1,"_",dds$diet2,"_",dds$subset))
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

results <- deresults(dds)

sdgnames <- function (sdg, dds=dds, vsd=vsd, rowMin=30, rowVars=2,...){
  filtered <- unique(rownames(sdg))
  filtered <- filtered[order( rowVars( assay(vsd)[filtered,]), decreasing=TRUE)]
  subset(filtered,
         rowMin(counts(dds)[filtered,]) > rowMin &
         rowVars(counts(dds)[filtered,]) > rowVars)
}

### Tissue Differential Expression ---------------------------------------------------------------------
res.EWAT_HiFat_CD4.MWAT_HiFat_CD4 <- results(ddsGroup, contrast=c("group","HiFat_EWAT_CD4","HiFat_MWAT_CD4"))

sdg.EWAT_Control_ILC2.MWAT_Control_ILC2 <- sdg (res.EWAT_Control_ILC2.MWAT_Control_ILC2)

sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4 <- sdgnames (sdg.EWAT_HiFat_CD4.MWAT_HiFat_CD4)


dietSDG <- unique(c(sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4,
                 ))

write.csv(counts(ddsGroup, normalized=TRUE)[dietSDG,colorder],file=paste0(date," ",expNum," Tissue Differentially Expressed Genes.csv" ))
dietSDGlist <- list(sdgnames.EWAT_HiFat_CD4.MWAT_HiFat_CD4,
                     )


ILC2tissueSDG <- c(tissueSDGlist[[3]], tissueSDGlist[[6]])
ILC2tissueSDG  <- ILC2tissueSDG[order( rowMeans( assay(vsd)[ILC2tissueSDG,]), decreasing=TRUE)]


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

stringGenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1E-1 RNAseq/3-1E-1 STRING/stringGenes.txt",
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
TGFBgenes <- read.table("~/Desktop/tgfb.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

TGFBgenes <- TGFBgenes[[1]]
TGFBgenes <- TGFBgenes[which(TGFBgenes %in% rownames(assay(vsd)))]
TGFBgenes <- TGFBgenes[which(rowMin(counts(dds)[TGFBgenes,])>50)]
TGFBgenes <- TGFBgenes[order( rowVars( assay(vsd)[TGFBgenes,]), decreasing=TRUE)]
sharedTGFBgenes <- TGFBgenes[which(TGFBgenes %in% topVarGenes)]
sharedTGFBgenes
resheatmap(vsdGroup,
           genes = TGFBgenes,
           samples= -grep("C_F_2_I", colnames(vsdGroup)),
           cex= 0.8,
           title= "TGFB Genes Heatmap")
dev.off()
resheatmap(vsd,
           genes = TGFBgenes,
           samples= colnames(vsd),
           cex= 0.8,
           title= "vsd TGFB Genes Heatmap")
dev.off()
resheatmap(vsd,
           genes = TGFBgenes,
           samples= -grep("T", colnames(vsd)),
           cex= 0.8,
           title= "ILC2 TGFB Genes Heatmap")
dev.off()

resheatmap(vsd,
           genes = TGFBgenes,
           samples= -grep("I", colnames(vsd)),
           cex= 0.8,
           title= "Treg TGFB Genes Heatmap")
dev.off()
## DEG Signature Genes ----------
signaturegenes <- read.table("~/Desktop/3-1E-1 signaturegenes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
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

