## Load required libraries-----
library("BiocParallel")
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("ggplot2")
library("gplots")
library("pheatmap")
save.image(paste0(dir,"/", date," ",expNum, " ILC2 Genes.RData"))

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

group <- paste0(diet," ",subset," ",tissue)

(coldata <- data.frame(row.names=colnames(countdata), diet, subset, tissue, group))
coldata <- coldata[order(coldata$subset,coldata$tissue, coldata$diet),]

## Instantiate DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~subset + diet + tissue)
dds <- DESeq(dds, parallel=TRUE)

## Import gene list of interest and export corresponding counts-----------
Selectedgenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/Genes.txt", 
                        quote="\"", comment.char="", stringsAsFactors=FALSE)
Selectedgenes <- Selectedgenes[[1]]
Selectedgenes <- Selectedgenes[which(Selectedgenes %in% rownames(counts(dds)))]
Selectedgenes <- Selectedgenes[order( rowMeans( assay(dds)[Selectedgenes,]), decreasing=TRUE)]

write.csv(counts(dds, normalized=TRUE)[Selectedgenes,colorder],
          file=paste0(date," ",expNum," Selected Genes Normalized Counts.csv" ))

Dupegenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/GenesDupes.txt", 
                            quote="\"", comment.char="", stringsAsFactors=FALSE)
Dupegenes <- Dupegenes[[1]]
Dupegenes <- Dupegenes[which(Dupegenes %in% rownames(counts(dds)))]
Dupegenes <- Dupegenes[order( rowMeans( assay(dds)[Dupegenes,]), decreasing=TRUE)]

write.csv(counts(dds, normalized=TRUE)[Dupegenes,colorder],
          file=paste0(date," ",expNum," Duplicated Selected Genes Normalized Counts.csv" ))

### Selected Genes Heatmaps ----------------------------------------------------
normcounts <- counts(dds, normalized=TRUE)[Selectedgenes, -grep("_CD4", colnames(dds))]
normheatmap(normcounts, 
           cex= 1,
           h=3,
           w=0.6,
           title= "Selected Genes Heatmap"
)
dev.off()

## Selected Gene Bar Plots --------------------- 
library(ggplot2)
library(gridExtra)
library(reshape2)

SelectCounts <- t(counts(dds, normalized=TRUE)[Selectedgenes,])
SelectCounts <- cbind2(coldata, SelectCounts)
meltedSC <- melt(SelectCounts, variable.name = "gene", value.name = "count" )
meltedSC$group <- paste0(meltedSC$subset, " ",meltedSC$tissue, " ",meltedSC$diet)
meltedSC$subsetdiet <- paste0(meltedSC$subset,"_",meltedSC$diet)
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

cairo_pdf(paste0(date," ",expNum," ", title, ".pdf"), w=8, h=20)
ggplot(data=meltedSC, aes(x=tissue, y=log2(count+1), fill=subsetdiet, color=subset)) +
  facet_wrap( ~ gene, scales="fixed", ncol = 6, dir="h") +
  stat_summary (fun.y = "mean", geom = "bar", position = position_dodge(0.9), width=0.75) +
  scale_fill_manual(values = subsetdietcolors) +
  scale_color_manual(values = subsetcolors ) +
  labs(title=paste0(expNum," ", title)) +
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.background = element_rect(color="#707070", linetype="solid", fill=NA, size=0.5),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = rel(2)))
dev.off()

title <- "Selected Genes"


dcast(meltedSC, )
