## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("ggplot2")
library("GMD")
library("pheatmap")

# Convert to matrix----------------------------------------------------
PPGrawcounts <- read.delim("2-2B-1 Raw Counts.txt", row.names=1)
countdata <- as.matrix(PPGrawcounts)
head(countdata)

# Assign condition ----------------------------------------------------
tissue <- c("Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Jejunum",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Lung",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen",
            "Spleen")
subset <- c("ILC1",
            "ILC1",
            "ILC1",
            "NK",
            "NK",
            "NK",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC2",
            "ILC2",
            "ILC2",
            "ILC3",
            "ILC3",
            "NK",
            "ILC2",
            "ILC2",
            "ILC2",
            "ILC1",
            "ILC1",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "ILC2",
            "ILC2",
            "ILC2",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC1",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK",
            "ILC1",
            "ILC1",
            "ILC2",
            "ILC2",
            "ILC3",
            "ILC3",
            "ILC3",
            "NK",
            "NK",
            "NK")

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet.
(coldata <- data.frame(row.names=colnames(countdata),tissue,subset))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~tissue+subset)

# Run DESeq normalization ----------------------------------------------------
dds<-DESeq(dds)
write.csv(counts(dds, normalized=TRUE),file="2-2B-1 Normalized Counts.csv")


## Regularized log and variance stabilizing transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds,fast=TRUE)
vsd <- varianceStabilizingTransformation(dds)
write.csv(assay(rld),file="2-2B-1 Log Transformed Normalized Counts.csv")


## PCA (subset only)
cairo_pdf("RO1 Subset PCA.pdf", w=6, h=4)
data <- plotPCA(rld, intgroup=c("subset"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=subset)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

#Key genes heatmaps-----------------
resheatmap <- function(rld, genes, samples, filename){
  with(rld,
       pheatmap(assay(rld) [genes, samples],
                cluster_rows=TRUE,
                scale="row",
                breaks = c(seq(-2, 2, length.out = 256)),
                border_color = NA,
                drop_levels=TRUE,
                color = my_palette,
                show_rownames=TRUE,
                cluster_cols=TRUE,
                annotation_col=coldata,
                annotation_colors = ann_colors,
                annotation_legend=FALSE,
                height=11,
                width=8.5,
                na.rm=TRUE,
                filename= filename)
       )
}
ann_colors = list(
  tissue = c(Spleen="orange", Lung="black", Jejunum="purple" )[coldata$tissue],
  subset = c(ILC1="red", ILC2="blue", ILC3="green", NK="gray" )[coldata$subset ])

resheatmap(rld, 
           genes = c("KLRB1",
                     "NCR1",
                     "TBX21",
                     "RORC",
                     "IFNG",
                     "CCR6",
                     "IL7R",
                     "IL23R",
                     "KIT",
                     "PTGDR2",
                     "IL1RL1",
                     "IL17RB",
                     "IL13"),                    
           samples=grep("Spleen_223_NK|Spleen_223_ILC3|Spleen_223_ILC1|Lung_223",colnames(rld)), 
           filename="Lung RO1 Key Genes Heatmap.pdf")
dev.off()

colnames(coldata)

