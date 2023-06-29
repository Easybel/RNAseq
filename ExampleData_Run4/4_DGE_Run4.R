## Neisseria gonorrhoeae transcriptome project

# This script takes the count table and analyses if everything is ok with the samples
# -- do they cluster?
# -- Do we need batch correction?

# input: count table
#  ----> this is the output from featureCounts
#  ----> this table only contains the CDS counts; 
#       for example: gene ID, NgoG4_1, _2, _3, ...

# we have 7 different phenotypes: 
# -- we want to compare the planktonic against biofilm **
# -- for both conditions we have 3 replicates 
#       -- these were made on the same day
# -- Batch: we assume that all samples are from the same batch --> but are they?

### install packages and load them
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggdendro)
library(ape)
library(ff)
library(ggfortify)
library(cluster)
library(ggrepel) 

setwd("/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/Map2_MS11Ref/")
inName = "Run4Ngo_CDS_Counts_G4_DpilE_NG17_NG24_NG32.csv"

#####################
## Ngo -- here the data is loaded and columns are added with sum, mean, median
data_raw <- read.table(inName, sep=" ", header = TRUE)
ID2geneInfo <- read.csv("/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/NgoBlast2_NmenANDFA1090_onlyCDS.csv",
                        sep="\t", header = TRUE)
data_samnum <- length(data_raw) -1

# information about the data
Conditions <- c(rep("G4wt",3),rep("G4pilEKO",3),rep("G4pilENG17",3),rep("G4pilENG24",3),rep("G4pilENG32",3))                                #set
replicate <- rep(c(1,2,3),data_samnum/3)                                                     #set
Batch <- rep(1,data_samnum)      

# set which conds is tested against which?
conds1 <- "G4pilENG32"          # condition that is tested against conds2 
conds2 <- "G4wt"    # reference 

# automatically generates the Sample Names
SamNames <- paste(Conditions,replicate,sep = "_")

# which conditions and biological replicates are you interested in? 
# if more than 1, add them with |
# all samples that should be taken into consideration in the analysis
cond_oi <- c("G4wt|G4pil") 
replicate_oi <- c("1|2|3")

# later on, all data is read into the DESeq2 function, and then the contrast is called

# search pattern
idx_cond_oi <- grep(cond_oi,Conditions)
idx_replicate_oi <- grep(replicate_oi,replicate)

idx_all     <- intersect(idx_cond_oi,idx_replicate_oi)

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
#batch <- Batch[idx_all]
replicate <- as.character(replicate[idx_all])

myColData <- data.frame(conds = conds, replicate = replicate)
rownames(myColData) <-  SamNames[idx_all]

dds <- DESeqDataSetFromMatrix(countData = data_mat_oi, colData = myColData, design =  ~  replicate + conds)

## be sure that you know your reference 
dds$conds <- relevel(dds$conds, ref=conds2)

##############################################################################################
################### Run DESeq2 and get the results
# run DESeq2: this wraps up the most important analysis steps
dds_Final <- DESeq(dds)

# generate the results tables - this function already does independant filtering, alpha is the FDR cutoff (default would be 0.1)
#?? or should I use independant hypothesis weighing? library("IHW"), resIHW <- results(dds, filterFun=ihw)
# The text: condition treated vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).
res <- results(dds_Final, contrast = c("conds",conds1,conds2), alpha = 0.05) #"conds",cond_oi,ctl_oi)) # default alpha=0.1summary(res)
summary(res)
mcols(res)$description # what means what in the results??

#?? should I also implement lfcShrink and why do I have a lfcSE column in my resutls?
resLFC <- lfcShrink(dds_Final, coef= paste("conds_",conds1,"_vs_",conds2,sep = ""), type="apeglm") # other options: ashr, normal (default)

resNamed <- data.frame(res[,]) # convert data from res to matrix for later use

# Test if annotation fits 
if(isFALSE(identical(rownames(resNamed),paste("gene-",ID2geneInfo[,1],sep = "")))){
  stop("The annotation does not fit to the imported data!! Please check")
}

resNamed$GN_in_Ngo <- ID2geneInfo[,2]
resNamed$GN_in_Nmen <- ID2geneInfo[,6]
resNamed$GN_in_FA1090 <- ID2geneInfo[,8]
resNamed$product <- ID2geneInfo[,4]
resNamed$locustag_in_Nmen <- ID2geneInfo[,5]
resNamed$locustag_in_FA1090 <- ID2geneInfo[,7]
resNamed_sort     <- resNamed[order(resNamed$padj),]

# write data to files
outName <- paste(paste("NgoG4Run4_",conds1,"_vs_",conds2,sep = ""),".tab",sep = "")
#write.table(resNamed_sort, paste("/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/RNA/2022March_Run4/1b_Analysis_naiv_woMS11_DpilT/2b_Results_DGE/",Sys.Date(),outName,sep = ""), sep="\t", quote = FALSE, na = "")

################## PLOT results
# MA-plot: shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
par(mfrow=c(1,2)); plotMA(res, ylim=c(-3,3));  title("MA-plot with LFC"); plotMA(resLFC, ylim=c(-3,3));  title("MA-plot with shrunken LFC");

# what the shrinkage does 
par(mfrow=c(1,1)); plot(res$log2FoldChange,resLFC$log2FoldChange); title("What the shrinkage does")

# plot the 9 most extreme 
par(mfrow=c(3,3))
for (i in c(1:9)){
  plotCounts(dds_Final, gene=rownames(resNamed_sort)[i], intgroup="conds")}

# Make a basic volcano plot
par(mfrow=c(1,1)); with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3) ,ylim =c(0,25)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resNamed_sort, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resNamed_sort, padj<.1 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

############### data transformation and visualization
# for DGE we used raw counts - but for visualization and clustering it can be useful to use transformed count data
# the two options are rlog (regularized logarithm) and vst (variance stabilizing transformations)
# look at the PCA
rld <- rlog(dds,blind=True) #?? why false and why rlog?
head(assay(rld), 3)
par(mfrow=c(1,1)); plotPCA(rld, intgroup="conds") 

# look at the heatmap 
# select <- order(rowMeans(counts(dds_Final,normalized=TRUE)),
#                 decreasing=TRUE)[1:100]
# df <- as.data.frame(colData(dds_Final)[,c("conds","batch")])
# pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)

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

