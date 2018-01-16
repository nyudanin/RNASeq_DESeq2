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
dir <- setwd("/Volumes/IBD/Yudanin/RNAseq/2-1A-4 RNAseq")
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
save.image(paste0(dir,"/", date," ",expNum, " ILC2 TGFb Genes.RData"))

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), tissue, replicate, day))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~tissue+replicate+day)

## Run DESeq normalization
dds<-DESeq(dds)



## Variance transformation for clustering/heatmaps, etc
vsd <- varianceStabilizingTransformation(dds)

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
tissuefill <- c(Colon="#E7A626", SI="#4266F6", MWAT="#269040", Lung="#C7302A" ) [colData(vsd)$tissue]
dayshapes <- c(D1=21, D2=22, D3=25) [colData(vsd)$day]
replicatecolor <- c(A="#212121", B="#000000") [colData(vsd)$replicate]

## Comparison to ILC2 TGFb Genes  --------------------- 
TGFbgenes <- read.table("/Volumes/IBD/Yudanin/RNAseq/3-1A-3 RNAseq/3-1A-3 Functional Analysis/ILC2 TGFb genes.txt", 
                        quote="\"", comment.char="", stringsAsFactors=FALSE)
TGFbgenes <- TGFbgenes[[1]]
TGFbgenes <- TGFbgenes[which(TGFbgenes %in% rownames(assay(vsd)))]
TGFbgenes <- TGFbgenes[order( rowMeans( assay(vsd)[TGFbgenes,]), decreasing=TRUE)]

### Selected Genes Heatmaps ----------------------------------------------------
TGFbnormcounts <- (counts(dds, normalized=TRUE)[TGFbgenes,])
TGFbnormcounts <- as.matrix(TGFbnormcounts)

normheatmap(TGFbnormcounts, 
            cex= 1,
            h=1.3,
            w=0.8,
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




