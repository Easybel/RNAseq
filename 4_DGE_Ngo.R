## Neisseria gonorrhoeae transcriptome project
# from featureCounts output 1 big matrix (table) is generated - here called: Ngo_CDS_Counts.csv
# - this table only contains the CDS counts
# - it contains 7 columns: gene ID, ctl1, ctl2, ctl3, azi1, azi2, azi3

### install packages and load them
library(DESeq2)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
setwd("/home/isabel/Dokumente/ExpEvol/Neisseria/Analysis/CDS/")

#####################
## Ngo -- here the data is loaded and columns are added with sum, mean, median
tab <- read.table("Ngo_CDS_Counts.csv", sep=" ", header = TRUE)
tab$sum <- apply(tab[,2:7],1, sum)
tab$mean <- apply(tab[,2:7],1, mean)
tab$med <- apply(tab[,2:7],1, median)
tab$var <- apply(tab[,2:7],1, var)

######################### plots prior to analysis
## a plot is opened in which the log2 distributions are plotted
par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(log2(tab$sum), ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(log2(tab$med), ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(log2(tab$mean), ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")
# a plot is opened in which the raw data distributions are plotted
par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(tab$sum, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(tab$med, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(tab$mean, breaks=1000, xlim=c(0,1000),
     ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")

################## Handle the data
# Pre-filter: Based on distributions the data is pre-filtered 
tab_red <- tab[which(tab$med > 9 & (tab$sum > 15 | tab$sum < 10000)),] #?? is this necessary??

# convert data in table to matrix
mat <- as.matrix(tab_red[,c(2,3,4,5,6,7)]) #here we have information on row names
head(mat)
rownames(mat) <- tab_red[,1] # but with this function we can add the row names to the matrix

# give information on samples, conditions, batch
conditions <- c(rep("ctl",3),rep("azi",3))
Batch <- rep(1,6)
myColData <- data.frame(condition = conditions, Batch = Batch)
rownames(myColData) <-  c("ctl1","ctl2","ctl3","azi1","azi2","azi3")
myColData
dds <- DESeqDataSetFromMatrix(countData = mat, colData = myColData, design =  ~ condition)
## be sure that you know your reference 
dds$condition <- relevel(dds$condition, ref="ctl")

################### Run DESeq2 and get the results
# run DESeq2: this wraps up the most important analysis steps
dds <- DESeq(dds)
#?? should I also implement lfcShrink and why do I have a lfcSE column in my resutls?
resLFC <- lfcShrink(dds, coef= "condition_azi_vs_ctl", type="apeglm") # other options: ashr, normal (default)

# generate the results tables - this function already does independant filtering, alpha is the FDR cutoff (default would be 0.1)
#?? or should I use independant hypothesis weighing? library("IHW"), resIHW <- results(dds, filterFun=ihw)
res <- results(dds, alpha=0.05) # ?? I changed this to 0.05, right?
summary(res)
mcols(res)$description # what means what in the results??

resOrdered <- res[order(res$padj),]
resOrdered_Mat <- as.matrix(res_sort[,c(1:6)]) # convert data from res to matrix for later use

# write data to files
write.table(resOrdered, "refctl_queriedAzi_results.csv", sep="\t")
strictCond <- resOrdered[which(res$padj <= 0.001 & res$baseMean > 10),]
write.table(strictCond, "refctl_queriedazi_results_padj1e-3.csv", sep="\t")

################## PLOT results
# MA-plot: shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
par(mfrow=c(1,2)); plotMA(res, ylim=c(-3,3));  title("MA-plot with LFC"); plotMA(resLFC, ylim=c(-3,3));  title("MA-plot with shrunken LFC");

# what the shrinkage does 
par(mfrow=c(1,1)); plot(res$log2FoldChange,resLFC$log2FoldChange); title("What the shrinkage does")

# plot the 9 most extreme 
par(mfrow=c(3,3))
for (i in c(1:9)){
plotCounts(dds, gene=rownames(resOrdered_Mat)[i], intgroup="condition")}

# Make a basic volcano plot
par(mfrow=c(1,1)); with(resOrdered, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resOrdered, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resOrdered, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

############### data transformation and visualization
# for DGE we used raw counts - but for visualization and clustering it can be useful to use transformed count data
# the two options are rlog (regularized logarithm) and vst (variance stabilizing transformations)
# look at the PCA
rld <- rlog(dds, blind=FALSE) #?? why false and why rlog?
head(assay(rld), 3)
par(mfrow=c(1,1)); plotPCA(rld, intgroup="condition") 

# look at the heatmap 
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition","Batch")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
