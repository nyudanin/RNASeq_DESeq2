# Assign conditions ----------------------------------------------------
names <- colnames(countdata)
getfield <- function (colname){
  strsplit(colname,"_")[[1]][1]
}
field <- unlist(lapply(names, getfield))

## DESeq Heatmaps ---------------------------------------------------------------------
ann_colors = list(
  subset = c( ILC1="#C7302A", ILC2="#4266F6", ILC3="#269040", NK="#E7A626" )[colData(vsd)$subset ],
  tissue = c( Spleen="black", Lung="light gray", Jejunum="#707070" )[colData(vsd)$tissue ])

resheatmap <- function(rld, genes, samples, cluster_cols=TRUE, height = 6, width = 6, filename, ...){
  filtered <- assay(rld) [genes, samples]
  filtered <- filtered[rowVars(filtered)>0,]
  with(rld,
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
                height=height,
                width=width,
                filename= filename)
  )
}

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

## DESeq MA Plot ----------------------------------------------------
maplot <- function (res, lfcthresh=2, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log10(baseMean), log2FoldChange, pch=20, cex=.5, ...))
  with(subset(res, padj<thresh & abs(log2FoldChange)>lfcthresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh & abs(log2FoldChange)>lfcthresh), textxy(baseMean, log2FoldChange, labs=row, cex=textcx, col=2))
  }
}

# Volcano plot with "significant" genes labeled----------------------------------------------------
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) 
{
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=row, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
