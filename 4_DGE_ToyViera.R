# Silene vulgaris - population Stranska Skala - root tissue
# prepare 1 big matrix (table) from all featureCounts outputs - here called: myTab.csv
# myTab.csv contains 7 columns: contig names, ctl1, ctl2, ctl3, azi1, azi2, azi3
# in rstudio or R
# install packages: DESeq2 and pheatmap

library(DESeq2)
library(pheatmap)
setwd("/home/isabel/Documents/Doktorarbeit/RNA_Seb/CDS/")

## Ngo -- 

tab <- read.table("Ngo_CDS_Counts.csv", sep=" ", header = TRUE)
tab$sum <- apply(tab[,2:7],1, sum)
tab$mean <- apply(tab[,2:7],1, mean)
tab$med <- apply(tab[,2:7],1, median)

par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(log2(tab$sum), ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(log2(tab$med), ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(log2(tab$mean), ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")

par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(tab$sum, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(tab$med, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(tab$mean, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")

# based on distributions, prefilter your data a bit 
# I do not use extremely covered transcripts and transcripts with less than 10 mapped reads 
tab_red <- tab[which(tab$med > 9 & (tab$sum > 15 | tab$sum < 10000)),]
head(tab_red)
mat <- as.matrix(tab_red[,c(2,3,4,5,6,7)])
rownames(mat) <- tab_red[,1]
head(mat)
conditions <- c(rep("SS_ctl",3),rep("SS_azi",3))
Batch <- rep(1,6)
myColData <- data.frame(condition = conditions, Batch = Batch)
rownames(myColData) <-  c("R_SS1_ctl","R_SS2_ctl","R_SS3_ctl","R_SS1_azi","R_SS2_azi","R_SS3_azi")
myColData
dds <- DESeqDataSetFromMatrix(countData = mat, colData = myColData, design =  ~ condition)
## be sure that you know your reference 
dds$condition <- relevel(dds$condition, ref="SS_ctl")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.table(res, "SS_refWT_queriedCopper_results.csv", sep="\t")
strictCond <- res[which(res$padj <= 0.001 & res$baseMean > 10),]
strictCond
write.table(res[which(res$padj <= 0.001 & res$baseMean > 10),], "SS_refWT_queriedCopper_results_padj1e-3.csv", sep="\t")
