# Silene vulgaris - population Stranska Skala - root tissue
# prepare 1 big matrix (table) from all featureCounts outputs - here called: myTab.csv
# myTab.csv contains 7 columns: contig names, ctl1, ctl2, ctl3, azi1, azi2, azi3
# in rstudio or R
# install packages: DESeq2 and pheatmap

library(DESeq2)
library(pheatmap)
setwd("path_to_your_files")

## ROOTS - SS - control condition versus copper treatment 

tab <- read.table("myTab.csv", sep="\t", header = TRUE)
tab$sum <- apply(tab[,2:7],1, sum)
tab$mean <- apply(tab[,2:7],1, mean)
tab$med <- apply(tab[,2:7],1, median)
par(mfrow=c(1,3))
hist(log2(tab$sum))
hist(log2(tab$med))
hist(log2(tab$mean))
# based on distriutions, prefilter your data a bit 
# I do not use extremely covered transcripts and transcripts with less than 10 mapped reads 
tab_red <- tab[which(tab$med > 9 & (tab$sum > 15 | tab$sum < 10000)),]
head(tab_red)
mat <- as.matrix(tab_red[,c(2,3,4,5,6,7)])
rownames(mat) <- tab_red[,1]
head(mat)
conditions <- c(rep("SS_wt",3),rep("SS_cu",3))
Batch <- rep(1,6)
myColData <- data.frame(condition = conditions, Batch = Batch)
rownames(myColData) <-  c("R_SS1_wt","R_SS2_wt","R_SS3_wt","R_SS1_cu","R_SS2_cu","R_SS3_cu")
myColData
dds <- DESeqDataSetFromMatrix(countData = mat, colData = myColData, design =  ~ condition)
## be sure that you know your reference 
dds$condition <- relevel(dds$condition, ref="SS_wt")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.table(res, "SS_refWT_queriedCopper_results.csv", sep="\t")
strictCond <- res[which(res$padj <= 0.001 & res$baseMean > 10),]
strictCond
write.table(res[which(res$padj <= 0.001 & res$baseMean > 10),], "SS_refWT_queriedCopper_results_padj1e-3.csv", sep="\t")



