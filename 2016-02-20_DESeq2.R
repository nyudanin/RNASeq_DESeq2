## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("ggplot2")
library("GMD")
library("pheatmap")
library("vsn")
library("Rtsne")


data_TM <- data_CM
data_TS <- data_CS

rm(data_CM)
rm(data_CS)

colnames(data_IM) <- paste("IM_", seq.int(1:length(data_IM)), sep="")
colnames(data_IS) <- paste("IS_", seq.int(1:length(data_IS)), sep="")
colnames(data_TM) <- paste("TM_", seq.int(1:length(data_TM)), sep="")
colnames(data_TS) <- paste("TS_", seq.int(1:length(data_TS)), sep="")


plot(colSums(data_IS)[1:200])
plot(colSums(data_TS)[1:200])
plot(colSums(data_IM)[1:200])
plot(colSums(data_TM)[1:200])

data_TS <- as.matrix(data_TS[,order(colSums(data_TS), decreasing=TRUE)])
data_IS <- as.matrix(data_IS[,order(colSums(data_IS), decreasing=TRUE)])
data_TM <- as.matrix(data_TM[,order(colSums(data_TM), decreasing=TRUE)])
data_IM <- as.matrix(data_IM[,order(colSums(data_IM), decreasing=TRUE)])
write.csv(data_IM, file = "800 MLN ILCs.csv")
write.csv(data_IS, file = "1900 Spleen ILCs.csv")
write.csv(data_TM, file = "1100 MLN T Cells.csv")
write.csv(data_IM, file = "3700 Spleen T Cells.csv")


quantile(colSums(data_IS, na.rm=TRUE))
quantile(colSums(data_TM, na.rm=TRUE))
quantile(colSums(data_TS, na.rm=TRUE))
quantile(colSums(data_IM, na.rm=TRUE))

quantile(rowSums(data_IS, na.rm=TRUE))
quantile(rowSums(data_TM, na.rm=TRUE))
quantile(rowSums(data_TS, na.rm=TRUE))
quantile(rowSums(data_IM, na.rm=TRUE))

plot(data_IM ["Il7r",1:200])
plot(data_IS ["Il7r",1:200])
plot(data_TM [c("Cd3e"),1:200])
plot(data_TS [c("Cd3e"),1:200])

MLNcounts <- merge.data.frame(data_TM[,1:300], data_IM[,1:300], by.x = 0, by.y = 0, all = FALSE)
row.names(MLNcounts) <- as.character(MLNcounts [,1])
MLNcounts <- MLNcounts [,-1]
MLNcounts <- MLNcounts[which(rowVars(MLNcounts)>0),which(MLNcounts["Actb",]>1)]
MLNcounts <- MLNcounts[,order(colSums(MLNcounts), decreasing=TRUE)]

SPLcounts <- merge.data.frame(data_TS[,1:300], data_IS[,1:300], by.x = 0, by.y = 0, all = FALSE)
row.names(SPLcounts) <- as.character(SPLcounts [,1])
SPLcounts <- SPLcounts [,-1]
SPLcounts <- SPLcounts[which(rowVars(SPLcounts)>0), which(SPLcounts["Actb",]>1)]
SPLcounts <- SPLcounts[,order(colSums(SPLcounts), decreasing=TRUE)]

ALLcounts <- merge.data.frame(MLNcounts, SPLcounts, by.x = 0, by.y = 0, all = FALSE)
row.names(ALLcounts) <- as.character(ALLcounts [,1])
ALLcounts <- ALLcounts [,-1]
ALLcounts <- ALLcounts[,order(colSums(ALLcounts), decreasing=TRUE)]

plot(colSums(ALLcounts, na.rm=TRUE))
plot(rowSums(ALLcounts, na.rm=TRUE))
quantile(colSums(ALLcounts, na.rm=TRUE))

# Assign conditions & colors-----------
subset <- unlist(substring(colnames(ALLcounts), 1, 1))
tissue <- unlist(substring(colnames(ALLcounts), 2, 2))

MLNsubset <- unlist(substring(colnames(MLNcounts), 1, 1))
SPLsubset <- unlist(substring(colnames(SPLcounts), 1, 1))

MLNcolor <- c(I="#C7302A", T="#4266F6")[MLNsubset]
SPLcolor <- c(I="#C7302A", T="#4266F6")[SPLsubset]

tissueColor <- c(M="#707070", S="#000000")[tissue]
subsetShapes <- c(I=21, T=22)[subset]
subsetColor <- c(I="#707070", T="#000000")[subset]
tissueShapes <- c(M=21, S=25)[tissue]

ALLgenecolor <- ifelse((unlist(ALLcounts["Il2rg",])) > 0,
                       "#C7302A",
                       ifelse((unlist(ALLcounts["Cd3g",])) > 0,
                              "#4266F6",
                              ifelse((unlist(ALLcounts["Ly6e",])) > 0,
                                     "#269040",
                                     ifelse((unlist(ALLcounts["Ccl5",])) > 0,
                                            "#E7A626",
                                            "#9E9E9E"))))

# tSNE -------------
rtsneALL <- Rtsne(unique(t(logALLcounts)), is_distance = FALSE)
rtsneMLN <- Rtsne(unique(t(logMLNcounts)), is_distance = FALSE)
rtsneSPL <- Rtsne(unique(t(logSPLcounts)), is_distance = FALSE)
ALLVarGenes <- ALLVarGenes[order(logALLcounts[ALLVarGenes,], decreasing=TRUE),]
cairo_pdf("5-1A-2 All Subset tSNE.pdf", w=6, h=6)
plot(rtsneALL$Y, 
     col=ALLgenecolor, 
     pch=tissueShapes, 
     bg=subsetColor, 
     cex=0.6, 
     xlab="tSNE 1", ylab="tSNE 2", main="5-1A-2 All Subset tSNE")
legend("bottomright",legend=c("Ccl5","Cd3g","Ly6e","Il7r"),fill=c("#C7302A","#4266F6", "#269040", "#E7A626"))
dev.off()

## Top Variable Genes ---------------
ALLcounts <- ALLcounts[,order(ALLcounts["Actb",], decreasing=TRUE)]
logALLcounts <- log10(ALLcounts+1)
logALLcounts <- logALLcounts[order(rowVars(logALLcounts), decreasing=TRUE),]
ALLVarGenes <- rownames(logALLcounts[(1:200),])

SPLcounts <- SPLcounts[,order(SPLcounts["Actb",], decreasing=TRUE)]
logSPLcounts <- log10(SPLcounts+1)
logSPLcounts <- logSPLcounts[order(rowVars(logSPLcounts), decreasing=TRUE),]
SPLVarGenes <- rownames(logSPLcounts[(1:200),])

MLNcounts <- MLNcounts[,order(MLNcounts["Actb",], decreasing=TRUE)]
logMLNcounts <- log10(MLNcounts+1)
logMLNcounts <- logMLNcounts[order(rowVars(logMLNcounts), decreasing=TRUE),]
MLNVarGenes <- rownames(logMLNcounts[(1:200),])

ALLpca <- prcomp(t(logALLcounts))
SPLpca <- prcomp(t(logSPLcounts))
MLNpca <- prcomp(t(logMLNcounts))

write.csv(logALLcounts[topvargenes,], file="logALLcounts Top Variable Genes.csv")
logALLcounts <- logALLcounts[,order(colSums(logALLcounts[topvargenes,]), decreasing=TRUE)]
