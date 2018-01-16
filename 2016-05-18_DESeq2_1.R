## Load required libraries-----
library("BiocParallel")
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("ggplot2")
library("gplots")
library("pheatmap")
save.image(paste0(dir,"/", date," ",expNum, " TGFb Genes.RData"))

## Set environment variables---------
register(MulticoreParam(4))
dir <- setwd("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 Functional Analysis/")
date <- paste0(Sys.Date())
expNum <- paste0("3-1A-3")
colorder <- c(4,7,16,19,1,10,13,5,8,17,20,2,11,14,6,9,18,3,12,15)

## Define functions and color palettes
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
resheatmap <- function(vsd, genes, samples, cluster_cols=TRUE, title=title, cex, h=1, w=1, scale="row",  ...){
  filtered <- assay(vsd) [genes, samples]
  filtered <- filtered[rowVars(filtered)>0,]
  with(vsd,
       pheatmap(filtered[order( rowVars( assay(vsd)[genes,]), decreasing=TRUE),],
                cex = cex,
                cluster_rows=TRUE,
                scale=scale,
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
  subset = c(Treg="#C7302A", ILC2="#707070", CD4="#4266F6" )[colData(dds)$subset],
  diet = c(HiFat="#269040", Control="#212121")[colData(dds)$diet],
  tissue = c(EWAT="#E7A626", MWAT="#9C27B0")[colData(dds)$tissue] )

## Import data-----------
countdata <- read.delim("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 DESeq2/3-1A-3 Raw Counts.txt", stringsAsFactors=FALSE)
countdata <- countdata[-(which(duplicated(countdata[,2])==TRUE)),]
row.names(countdata) <- unlist(countdata$Symbol)

## Convert to matrix 
countdata <- as.matrix(countdata[,-(1:2)])

## Assign conditions and create coldata metadataframe-----------
names <- colnames(countdata)
getdiet <- function (colname){
  strsplit(colname,"_")[[1]][2]
}
diet <- unlist(lapply(names, getdiet))
diet <- gsub("C","Control",diet)
diet <- gsub("F","HiFat",diet)

gettissue <- function (colname){
  strsplit(colname,"_")[[1]][3]
}
tissue <- unlist(lapply(names, gettissue))
tissue <- gsub("E","EWAT",tissue)
tissue <- gsub("M","MWAT",tissue)

getsubset <- function (colname){
  strsplit(colname,"_")[[1]][5]
}
subset <- unlist(lapply(names, getsubset))
rm(names)

group <- paste0(subset," ",tissue," ",diet)

(coldata <- data.frame(row.names=colnames(countdata), diet, subset, tissue, group))
coldata <- coldata[order(coldata$subset,coldata$tissue, coldata$diet),]

## Instantiate DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~subset + diet + tissue)
dds <- DESeq(dds, parallel=TRUE)

## Import gene list of interest and export corresponding counts-----------
ILC2genes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/ILC2_genes.txt", 
                          quote="\"", comment.char="", stringsAsFactors=FALSE)
ILC2genes <- ILC2genes[[1]]
ILC2genes <- ILC2genes[which(ILC2genes %in% rownames(counts(dds)))]
ILC2genes <- ILC2genes[order( rowMeans( assay(dds)[ILC2genes,]), decreasing=TRUE)]

write.csv(counts(dds, normalized=TRUE)[ILC2genes,colorder],
          file=paste0(date," ",expNum," ILC2 Genes Normalized Counts.csv" ))
write.csv(counts(dds, normalized=FALSE)[ILC2genes,colorder],
          file=paste0(date," ",expNum," ILC2 Genes Raw Counts.csv" ))

TGFbgenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 Functional Analysis/Genes.txt", 
                        quote="\"", comment.char="", stringsAsFactors=FALSE)
TGFbgenes <- TGFbgenes[[1]]
TGFbgenes <- unique(TGFbgenes)
TGFbgenes <- TGFbgenes[which(TGFbgenes %in% rownames(counts(dds)))]
TGFbgenes <- TGFbgenes[which(rowMin(counts(dds)[TGFbgenes,])>100)]
TGFbgenes <- TGFbgenes[order( rowVars( assay(dds)[TGFbgenes,]), decreasing=TRUE)]
TGFbgenes <- TGFbgenes[which(rowVars(counts(dds)[TGFbgenes,])>10)]
### Selected Genes Heatmaps ----------------------------------------------------
ILC2normcounts <- counts(dds, normalized=TRUE)[ILC2genes[1:30],]
TGFbnormcounts <- counts(dds, normalized=TRUE)[TGFbgenes,]

normheatmap(ILC2normcounts[,grep("ILC2",colnames(dds))], 
           cex= 1,
           h=1.5,
           w=0.5,
           title= "ILC2 Genes Heatmap"
)
dev.off()

normheatmap(TGFbnormcounts[,grep("ILC2",colnames(dds))], 
            cex= 1,
            h=1.5,
            w=0.5,
            title= "ILC2 TGFb Genes Heatmap"
)
dev.off()

## Selected Gene Bar Plots --------------------- 
library(ggplot2)
library(gridExtra)
library(reshape2)

TGFbcounts <- t(counts(dds, normalized=TRUE)[TGFbgenes,])
TGFbcounts <- cbind2(coldata, TGFbcounts)
meltedSC <- melt(TGFbcounts, variable.name = "gene", value.name = "count" )
meltedSC$group <- paste0(meltedSC$subset, " ",meltedSC$tissue, " ",meltedSC$diet)
meltedSC$subsetdiet <- paste0(meltedSC$subset,"_",meltedSC$diet)

subsetdietcolors <- c(Treg_Control="#E39794", 
  ILC2_Control="#d4d4d4", 
  CD4_Control="#b3c1fb",
  Treg_HiFat="#C7302A", 
  ILC2_HiFat="#707070", 
  CD4_HiFat="#4266F6")[meltedSC$subsetdiet]
  
  
subsetcolors <- c(Treg="#4f1310", 
  ILC2="#212121", 
  CD4="#131e49")[meltedSC$subset]

cairo_pdf(paste0(date," ",expNum," ILC2 TGFb Genes.pdf"), w=10, h=13)
ggplot(data=meltedSC, aes(x=group, y=log2(count+1), fill=subsetdiet, color=subset)) +
  facet_wrap( ~ gene, scales="fixed", ncol = 6, dir="h") +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.75) +
  scale_fill_manual(values = subsetdietcolors) +
  scale_color_manual(values = subsetcolors ) +
  labs(title=paste0(expNum," ILC2 TGFb Genes")) +
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.background = element_rect(color="#707070", linetype="solid", fill=NA, size=0.5),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = rel(2)))
dev.off()


