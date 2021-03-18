## Neisseria gonorrhoeae transcriptome project

# This script takes the count table and analyses if everything is ok with the samples
# -- do they cluster?
# -- Do we need batch correction?

# input: count table
#  ----> this is the output from featureCounts
#  ----> this table only contains the CDS counts; 
#       for example: gene ID, ctl1, ctl2, ctl3, azi1, azi2, azi3

# we have 3 (4) different  
# -- ** condition: treatment (antibiotica: Azi or Cef) of no treatment (DMSO) **
# -- duration of treatment/ treatement time: 15 or 60 min 
# -- day of treatment: for each treatment + duration we have 3 replicates
#       -- these were made on different days (day 1,2, and 3) 
# -- Batch: we assume that all samples are from the same batch --> but are they?

### install packages and load them
library(DESeq2)
library(pheatmap)
#library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(ggdendro)
library(ape)
library(ff)
library(ggfortify)
library(cluster)
library(ggrepel) 

setwd("/home/isabel/sciebo/P5_ExpEvol_Ngo/RNA/2021Feb_Run2/1_Result_MapCount/")
inName = "NgoG4_Run2_CDS_allCounts.csv"
# setwd("/home/isabel/sciebo/P5_ExpEvol_Ngo/RNA/2020_Run1/data/")
# inName = "Ngo_CDS_Counts.csv"

#####################
## Ngo -- here the data is loaded and columns are added with sum, mean, median
data_raw <- read.table(inName, sep=" ", header = TRUE)
ID2geneInfo <- read.csv("/home/isabel/sciebo/P5_ExpEvol_Ngo/dictionaries/MS11_ncbi_IR/Manual_Lists/NgoBlast2_NmenANDFA1090_onlyCDS.csv",
                        sep="\t", header = TRUE)
data_samnum <- length(data_raw) -1

# # information about the data
# # all initial infos that apply to all are capitalized to distinguish from later use
Conditions <- c(rep("DMSO",6),rep("Azi",6),rep("Cef",6))                                #set
TreatTimes <- c(rep("15",3),rep("60",3),rep("15",3),rep("60",3),rep("15",3),rep("60",3)) #set
Days <- rep(c(1,2,3),data_samnum/3)                                                     #set
Batch <- rep(1,data_samnum)                                                             #set
                                                          

# automatically generates the Sample Names
SamNames <- paste(paste(Conditions,TreatTimes,sep = "_"),Days,sep = "min_")

# which treatment durations, days and conditions are you interested in? 
# if more than 1, add them with |
cond_oi <- c("Cef")
ctl_oi <- c("DMSO")
dur_oi <- c("60")
day_oi <- c("1|2|3")

# search pattern
idx_dur_oi <- grep(dur_oi,TreatTimes)
idx_cond_oi <- grep(cond_oi,Conditions)
idx_day_oi <- grep(day_oi,Days)
idx_ctl_oi <- grep(ctl_oi,Conditions)

idx_ctl     <- intersect(intersect(idx_dur_oi,idx_ctl_oi),idx_day_oi)
idx_samples <- intersect(intersect(idx_dur_oi,idx_cond_oi),idx_day_oi)
idx_all     <- sort(c(idx_ctl,idx_samples))

outName <- paste(paste(paste("NgoG4run2_",cond_oi,dur_oi,"min",sep = ""),paste(ctl_oi,dur_oi,"min",".tab",sep = ""),sep = "_vs_"))
########################################################################################
##### Look at means and sums for all genes

data_raw$sum  <- apply(data_raw[,2:data_samnum],1, sum)
data_raw$mean <- apply(data_raw[,2:data_samnum],1, mean)
data_raw$med  <- apply(data_raw[,2:data_samnum],1, median)
data_raw$var  <- apply(data_raw[,2:data_samnum],1, var)

######################### plots prior to analysis
## a plot is opened in which the log2 distributions are plotted
par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(log2(data_raw$sum), ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(log2(data_raw$med), ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(log2(data_raw$mean), ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")
# a plot is opened in which the raw data distributions are plotted
par(mfrow=c(1,3)) #combines plots, in 1 line and 3 columns
hist(data_raw$sum, breaks=1000, xlim=c(0,20000),
     ylab="Genes", xlab="log of sum count, all conditions", main="Distribution of counts")
hist(data_raw$med, breaks=1000, xlim=c(0,20000),
     ylab="Genes", xlab="log of median count, all conditions", main="Distribution of counts")
hist(data_raw$mean, breaks=1000, xlim=c(0,20000),
     ylab="Genes", xlab="log of mean count, all conditions", main="Distribution of counts")

#############################################################################################################
################## Process the data

# Pre-filter: Based on distributions the data is pre-filtered 
data_filt <- data_raw #[which(data_raw$med > 9 & (data_raw$sum > 15 | data_raw$sum < 10000)),]

# convert data in table to matrix
data_mat <- as.matrix(data_filt[,(1:data_samnum)+1]) 

rownames(data_mat) <- data_filt[,1] 
colnames(data_mat) <- SamNames

################## Create the dds object for the analysis

# give information on samples, conditions, batch

data_mat_oi <- data_mat[,idx_all]

conds <- as.factor(Conditions[idx_all])
batch <- Batch[idx_all]
days <- as.character(Days[idx_all])
dur <- as.character(TreatTimes[idx_all])

myColData <- data.frame(conds = conds, batch = batch, days = days, dur = dur)
rownames(myColData) <-  SamNames[idx_all]

dds <- DESeqDataSetFromMatrix(countData = data_mat_oi, colData = myColData, design =  ~ conds)

## be sure that you know your reference 
dds$conds <- relevel(dds$conds, ref="DMSO")

##############################################################################################
################### Run DESeq2 and get the results
# run DESeq2: this wraps up the most important analysis steps
dds_Final <- DESeq(dds)
#?? should I also implement lfcShrink and why do I have a lfcSE column in my resutls?
#resLFC <- lfcShrink(dds_Final, coef= "condition_azi_vs_ctl", type="apeglm") # other options: ashr, normal (default)

# generate the results tables - this function already does independant filtering, alpha is the FDR cutoff (default would be 0.1)
#?? or should I use independant hypothesis weighing? library("IHW"), resIHW <- results(dds, filterFun=ihw)
res <- results(dds_Final) # default alpha=0.1
summary(res)
mcols(res)$description # what means what in the results??

resNamed <- data.frame(res[,]) # convert data from res to matrix for later use
resNamed$GN_in_Ngo <- ID2geneInfo[,2]
resNamed$GN_in_Nmen <- ID2geneInfo[,5]
resNamed$GN_in_FA1090 <- ID2geneInfo[,6]
resNamed$Product <- ID2geneInfo[,4]
resNamed_sort     <- resNamed[order(resNamed$padj),]

# write data to files
write.table(resNamed_sort, paste("../2b_Results_DGE/",outName,sep = ""), sep="\t")

################## PLOT results
# MA-plot: shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
par(mfrow=c(1,2)); plotMA(res, ylim=c(-3,3));  title("MA-plot with LFC"); plotMA(resLFC, ylim=c(-3,3));  title("MA-plot with shrunken LFC");

# what the shrinkage does 
par(mfrow=c(1,1)); plot(res$log2FoldChange,resLFC$log2FoldChange); title("What the shrinkage does")

# plot the 9 most extreme 
par(mfrow=c(3,3))
for (i in c(1:9)){
  plotCounts(dds, gene=rownames(resNamed_sort)[i], intgroup="condition")}

# Make a basic volcano plot
par(mfrow=c(1,1)); with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resNamed_sort, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resNamed_sort, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

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

