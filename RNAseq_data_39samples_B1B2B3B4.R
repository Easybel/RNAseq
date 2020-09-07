setwd("/home/veve/Documents/Bsubt_HGT/improvedRNAseq4")
library(DESeq2)
library(ggplot2)
library(ggdendro)
library(ape)
library(ff)
library(ggfortify)
library(cluster)
library(ggrepel)
library(pheatmap)

tab <- read.table("Bsubt_39samples.counts", header = TRUE, sep="\t")
colnames(tab)

myMat <- as.matrix(tab[,c(2:40)])
rownames(myMat) <- tab[,1]
head(myMat)
colnames(myMat)
sumSamples <- apply(myMat, 2, sum)
par(mar=c(7,2.5,2,0.5), mgp=c(1.5,0.5,0), cex.lab=0.85, cex.axis=0.85, cex.main=0.9)
plot(x=c(1:39), y=sumSamples, type="b", xlab="", ylab="raw read counts", main="Raw reads, 39 samples", xaxt="n")
axis(1, at=c(1:39), labels=c("ancestral_a_1", "ancestral_a_2", "ancestral_a_3","ancestral_b_1", "ancestral_b_2", "ancestral_b_3", "c21_noDNA_a_1", "c21_noDNA_a_2", "c21_noDNA_b_1", "c21_noDNA_b_2", "c21_noDNA_c_1", "c21_noDNA_c_2", "c21_DNA_W1_1", "c21_DNA_W1_2", "c21_DNA_W2_1", "c21_DNA_W2_2", "c21_DNA_W4_1", "c21_DNA_W4_2", "c21_DNA_W4_3",  "c21_DNA_W5_1", "c21_DNA_W5_2", "c21_DNA_W5_3", "c21_DNA_W5_4", "c21_DNA_W5_5", "c21_DNA_W6_1", "c21_DNA_W6_2", "c21_DNA_W7_1", "c21_DNA_W7_2", "c21_DNA_W8_1", "c21_DNA_W8_2", "c21_DNA_W1_3", "c21_DNA_W1_4", "c21_DNA_W1_5", "c21_noDNA_b_3", "c21_noDNA_b_4", "c21_noDNA_a_3", "c21_noDNA_a_4", "c21_noDNA_c_3", "c21_noDNA_c_4"), las=2)
# png file name: 39samples_sum_RawReads
# we have ancestrals only in the batch 1 !!! I can not wisely remove batch effects from them. 
# here, we see differences in library see - we need to normalized for them

### Normalization, DESeq2 pck. ###

condition <- c(rep("ancestral_a",3),rep("ancestral_b",3), "c21_noDNA_a", "c21_noDNA_a", "c21_noDNA_b", "c21_noDNA_b", "c21_noDNA_c", "c21_noDNA_c","W1","W1","W2","W2",rep("W4",3),rep("W5",5),"W6","W6","W7","W7", "W8","W8", rep("W1",3), "c21_noDNA_b", "c21_noDNA_b", "c21_noDNA_a", "c21_noDNA_a", "c21_noDNA_c", "c21_noDNA_c")
Batch <- c(rep("B1",6), rep("B2",2), rep("B3",4), rep("B2",4), rep("B1",6), rep("B2",8), rep("B4", 9))
myColData <- data.frame(condition = condition, Batch = Batch)
rownames(myColData) <- colnames(myMat)
myColData
# you can save the table with data descriptions: 
# write.table(myColData, "Data_description.csv", sep="\t")

dds <- DESeqDataSetFromMatrix(countData = myMat, colData = myColData, design =  ~ condition)
# normalization - I am using rlog function
# there is also another normalization type called "vst" in DEseq2
# using the mean of gene-wise disperion estimates as the fitted value: fitType="mean"
# the local regression does better, as it has more degrees of freedom compared to the parametric with only two parameters: fitType="local"
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#local-or-mean-dispersion-fit
myRLOG <- rlog(dds)
head(assay(myRLOG))
# write.table(assay(myRLOG), "Bsubt_39samples_rlogs.csv", sep="\t")
rlogMat <- assay(myRLOG)
head(rlogMat)

### tests ###

# R-R plots (replicates)
head(rlogMat[,c(7,36)])
# the same sample which differs in batch
par(mfrow=c(2,1), mar=c(2.5,2.5,2,0.5), mgp=c(1.5,0.5,0), cex.lab=0.75, cex.axis=0.75, cex.main=0.85)
plot(myMat[,7], myMat[,36], xlab="N17_B1", ylab="N17_B4", main="raw counts")
cor.test(myMat[,7], myMat[,36])
text(x=10000, y=100000, "r=0.697")
cor.test(rlogMat[,7], rlogMat[,36])
plot(rlogMat[,7], rlogMat[,36], xlab="N17_B1", ylab="N17_B4", main="normalized values")
text(x=0, y=12, "r=0.96")

# heatmaps (dendrograms)
pheatmap(myMat, cluster_rows = FALSE, labels_col = paste(condition, Batch, sep="_"), show_rownames = FALSE, cellheight = 0, legend = FALSE, main="raw counts")
# even if we use raw counts, we see that biological replicates cluster together
# clustering does not follow conditions
# png file: Heatmap_rowCounts_39samples
pheatmap(rlogMat, cluster_rows = FALSE, labels_col = paste(condition, Batch, sep="_"), show_rownames = FALSE, cellheight = 0, legend = FALSE, main="normalized values")
# clustering improved after rlog normalization but stil follow batch numbers
# png file: Heatmap_normalizedValues_39samples

# the plotPCA function is implemented within DESeq2 pckg. and uses spectral decomposition 
plotPCA(myRLOG, intgroup=c("Batch")) + ggtitle("PCA - Bsubtilis - rlog normalized counts")
#  png file: PCA_rlogs_normalizedValues_byBatch


### BATCH CORRECTIOn ###
# https://rdrr.io/bioc/sva/man/ComBat.html
# https://academic.oup.com/biostatistics/article/17/1/29/1744261
# The input data are assumed to be cleaned and normalized before batch effect removal. 
library(sva)
myCol <- c(rep("darksalmon",6), rep("forestgreen",2), rep("steelblue", 4), rep("forestgreen",4), rep("darksalmon",6), rep("forestgreen",8), rep("maroon2",9))
# remove genes with 
rlogMatNo0 <- rlogMat[-unlist(which(apply(rlogMat,1,sum) == 0)),]
# 6 genes have sum of rlog values zero; these must be removed
# mod: Model matrix for outcome of interest and other covariates besides batch
# par.prior = FALSE: non-parametric empirical Bayes frameworks for adjusting data for batch effects
# rior.plots = TRUE: give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
# par.prior = TRUE:  indicates parametric adjustments will be used
rlogComb <- ComBat(rlogMatNo0, batch=Batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
# write.table(rlogComb, "Bsubt_39samples_rlogs_BatchEffAdj_ParametricCombat.csv", sep="\t")

# PCA
autoplot(prcomp(t(rlogMat)), x=1, y=2, label = FALSE) + theme_classic() + ggtitle("39 samples - normalized values - PC1 and PC2 - coloured by batches") + geom_point(color = myCol, size = 1) + theme(plot.title = element_text(size=10,face="bold")) + geom_text_repel(aes(label = paste(condition, Batch, sep="_")), size=2.5, nudge_y=0.03, direction = "x", angle="90", segment.size = 0, show.legend = FALSE, colour=myCol)
# png file: SVD_normalizedValues_39samples
autoplot(prcomp(t(rlogComb)), x=1, y=2, label = FALSE) + theme_classic() + ggtitle("39 samples - normalized values, removed batch eff. - PC1 and PC2 \ncoloured by batches") + geom_point(color = myCol, size = 1) + theme(plot.title = element_text(size=10,face="bold")) + geom_text_repel(aes(label = paste(condition, Batch, sep="_")), size=2.5, nudge_y=0.03, direction = "x", angle="90", segment.size = 0, show.legend = FALSE, colour=myCol)
# png file: SVD_normalizedValues_BatchEffCorrected_39samples

# heatmap
pheatmap(rlogComb, cluster_rows = FALSE, labels_col = paste(condition, Batch, sep="_"), show_rownames = FALSE, cellheight = 0, legend = FALSE, main="normalized values after the correction for batch effects")
# clustering looks good

# NEED TO BE REWRITTEN
library(vioplot)
par(mfrow=c(2,1), mar=c(1.5,5,2,0.5), cex.main=1, cex.lab=0.75, cex.axis=0.75,  mgp=c(1.5, 0.5, 0), las=1)
vioplot(apply(rlogMat[,1:6], 1, mean), apply(myRlogs[,7:8], 1, mean), apply(myRlogs[,9:12], 1, mean), names=c("c0 noDNA B1", "c21 noDNA B2", "c21 noDNA B3"),  col=c("cornsilk2", "salmon", "steelblue"), main="Rlogs with batch effects", range=2, horizontal=TRUE)
vioplot(apply(rlogComb[,7:8], 1, mean), apply(matNoBatchEff[,7:8], 1, mean), apply(matNoBatchEff[,9:12], 1, mean), names=c("c0 noDNA B1", "c21 noDNA B2", "c21 noDNA B3"),  col=c("cornsilk2", "salmon", "steelblue"), main="Rlogs with NO batch effects; nonparametric combat", range=2, horizontal=TRUE)

### Differences btw c21 noDNA and c21 +DNA ?
# linear regression; because we do not have count data anymore, we need to use an approach for continous data (not implemented in DESeq2)
colnames(rlogComb)
allPval <- c()
L2FCH <- c()
for (i in seq(1:dim(rlogComb)[1])){
  # lm is used to fit linear models. It can be used to carry out regression.
  # conditions: nonUV versus UV
  M <- lm(rlogComb[i,c(7:12,13:33,34:39)] ~ c(rep("ref",6),rep("treated",21),rep("ref",6)))
  # to test differences between "reference" group and "treated" group, I use anova with Fisher test
  allPval <- c(allPval, anova(M, test = "F")[[5]][1] )
  A <- mean(rlogComb[i,c(7:12,34:39)])
  B <- mean(rlogComb[i,c(13:33)])
  # log2fch = treated/reference = tr/wt
  L2FCH <- c(L2FCH, B-A)
}
# p-values will be adjusted
# here, I do the bonferroni multiple-testing correction
# to give strong control of the family-wise error rate
# rlogsTab$Pval_day0noDNA_vs_day21 <-  p.adjust(allPval, method="bonferroni")
# we can also do the fdr correction prefered by Adreas Beyer's group
PvalsFDRs <-  p.adjust(allPval, method="fdr")
noDNAs_vs_HGT <- as.data.frame(cbind(A, B, allPval,PvalsFDRs, L2FCH))
rownames(noDNAs_vs_HGT) <- rownames(rlogComb)
colnames(noDNAs_vs_HGT) <- c("meanRef","meanTreatm","pvalues","padjFDR","log2FCH")
head(noDNAs_vs_HGT)
# write.table(noDNAs_vs_HGT, "GeneTranscr_c21_noDNA_vs_HGT_rlogs_combat.csv", sep="\t")

length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & ( as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 | as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 )))
# 17
length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 ))
# 17
length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 ))
# 0
DEG <- noDNAs_vs_HGT[which( (noDNAs_vs_HGT$log2FCH >= 1 | noDNAs_vs_HGT$log2FCH <= -1) & noDNAs_vs_HGT$padjFDR < 0.05 ),]
# write.table(DEG, "DEgenes_c21_noDNA_vs_HGT_rlogs_combat.csv", sep="\t")
pheatmap(rlogComb[which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & ( as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 | as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 )), c(7:12,13:33,34:39)], labels_col = paste(condition, Batch, sep="_")[c(7:12,13:33,34:39)], fontsize_row = 7, clustering_distance_rows = "euclidean")

### Differences btw c0 noDNA (ancestrals) and c21 +DNA ?
# linear regression; because we do not have count data anymore, we need to use an approach for continous data (not implemented in DESeq2)
colnames(rlogComb)
allPval <- c()
L2FCH <- c()
for (i in seq(1:dim(rlogComb)[1])){
  # lm is used to fit linear models. It can be used to carry out regression.
  # conditions: nonUV versus UV
  M <- lm(rlogComb[i,c(1:6,13:33)] ~ c(rep("ref",6),rep("treated",21)))
  # to test differences between "reference" group and "treated" group, I use anova with Fisher test
  allPval <- c(allPval, anova(M, test = "F")[[5]][1] )
  A <- mean(rlogComb[i,c(1:6)])
  B <- mean(rlogComb[i,c(13:33)])
  # log2fch = treated/reference = tr/wt
  L2FCH <- c(L2FCH, B-A)
}
# p-values will be adjusted
# here, I do the bonferroni multiple-testing correction
# to give strong control of the family-wise error rate
# rlogsTab$Pval_day0noDNA_vs_day21 <-  p.adjust(allPval, method="bonferroni")
# we can also do the fdr correction prefered by Adreas Beyer's group
PvalsFDRs <-  p.adjust(allPval, method="fdr")
noDNAs_vs_HGT <- as.data.frame(cbind(A, B, allPval,PvalsFDRs, L2FCH))
rownames(noDNAs_vs_HGT) <- rownames(rlogComb)
colnames(noDNAs_vs_HGT) <- c("meanRef","meanTreatm","pvalues","padjFDR","log2FCH")
head(noDNAs_vs_HGT)
# write.table(noDNAs_vs_HGT, "GeneTranscr_c0_noDNA_vs_HGT_rlogs_combat.csv", sep="\t")

length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & ( as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 | as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 )))
# 150
length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 ))
# 131
length(which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 ))
# 19
# why such low number of upregulated genes in treated samples?
DEG <- noDNAs_vs_HGT[which( (noDNAs_vs_HGT$log2FCH >= 1 | noDNAs_vs_HGT$log2FCH <= -1) & noDNAs_vs_HGT$padjFDR < 0.05 ),]
# write.table(DEG, "DEgenes_c0_noDNA_vs_HGT_rlogs_combat.csv", sep="\t")
pheatmap(rlogComb[which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & ( as.double(as.character(noDNAs_vs_HGT$log2FCH))  <= -1 | as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 )), c(1:6,13:33)], labels_col = paste(condition, Batch, sep="_")[c(1:6,13:33)], show_rownames = FALSE, clustering_distance_rows = "euclidean")
pheatmap(rlogComb[which(as.double(as.character(noDNAs_vs_HGT$padjFDR)) < 0.05 & as.double(as.character(noDNAs_vs_HGT$log2FCH)) >= 1 ), c(1:6,13:33)], labels_col = paste(condition, Batch, sep="_")[c(1:6,13:33)], clustering_distance_rows = "euclidean" )
                                                                                  
                        