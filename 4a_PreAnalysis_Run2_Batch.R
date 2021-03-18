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
# Conditions <- c(rep("DMSO",6),rep("Azi",6),rep("Cef",6))                                #set
# TreatTimes <- c(rep("15",3),rep("60",3),rep("15",3),rep("60",3),rep("15",3),rep("60",3)) #set
# Days <- rep(c(1,2,3),data_samnum/3)                                                     #set
# Batch <- rep(1,data_samnum)                                                             #set
# 
# # automatically generates the Sample Names
# SamNames <- paste(paste(Conditions,TreatTimes,sep = "_"),Days,sep = "min_")
# 
# # which treatment durations, days and conditions are you interested in? 
# # if more than 1, add them with |
# cond_oi <- c("Azi")
# ctl_oi <- c("DMSO")
# dur_oi <- c("15|60")
# day_oi <- c("1|2|3")

# information about the data
# all initial infos that apply to all are capitalized to distinguish from later use
Conditions <- c(rep("DMSO",6),rep("Azi",6),rep("Cef",6))                                #set
TreatTimes <- c(rep("15",3),rep("60",3),rep("15",3),rep("60",3),rep("15",3),rep("60",3)) #set
Days <- rep(c(1,2,3),data_samnum/3)                                                     #set
Batch <- rep(1,data_samnum)                                                             #set

# automatically generates the Sample Names
SamNames <- paste(paste(Conditions,TreatTimes,sep = "_"),Days,sep = "min_")

# which treatment durations, days and conditions are you interested in? 
# if more than 1, add them with |
cond_oi <- c("Azi|Cef")
ctl_oi <- c("DMSO")
dur_oi <- c("15|60")
day_oi <- c("1|2|3")

# search pattern
idx_dur_oi <- grep(dur_oi,TreatTimes)
idx_cond_oi <- grep(cond_oi,Conditions)
idx_day_oi <- grep(day_oi,Days)
idx_ctl_oi <- grep(ctl_oi,Conditions)

# these are the samples that don't have super low mRNA fraction
# idx_ctl <- c(1,3,5,6)
# idx_oi <- c(7,9,10,11,14,15,17,18)

idx_ctl     <- intersect(intersect(idx_dur_oi,idx_ctl_oi),idx_day_oi)
idx_samples <- intersect(intersect(idx_dur_oi,idx_cond_oi),idx_day_oi)
idx_all     <- sort(c(idx_ctl,idx_samples))

#######################################################################################
############### check the number of counts per sample -- analyse if there are connections

counts_per_sample <- apply(data_raw[,(1:data_samnum)+1],2,sum)
counts <- data.frame(counts_per_sample,TreatTimes,Conditions,Days)
counts$Days <- as.factor(counts$Days)
counts$Conditions <- as.factor(counts$Conditions)
counts$TreatTime <- as.factor(counts$TreatTime)


# Basic box plot
counts_depCond <-ggplot(counts, aes(x=Conditions, y=counts_per_sample, fill = Conditions)) +
  geom_boxplot(notch=FALSE)
counts_depCond

counts_depDur <-ggplot(counts, aes(x=Days, y=counts_per_sample, fill = TreatTime)) +
  geom_boxplot(notch=FALSE) 

counts_depDur

#############################################################################################################
################## Process the data

# Pre-filter: Based on distributions the data is pre-filtered 
data_filt <- data_raw #[which(tab$med > 9 & (tab$sum > 15 | tab$sum < 10000)),]

# convert data in table to matrix
data_mat <- as.matrix(data_filt[,(1:data_samnum)+1]) 

rownames(data_mat) <- data_filt[,1] 
colnames(data_mat) <- SamNames

################## Create the dds object for the analysis

# give information on samples, conditions, batch

data_mat_oi <- data_mat[,idx_all]

conds <- as.factor(Conditions[idx_all])
batch <- Batch[idx_all]
days <- as.factor(Days[idx_all])
dur <- as.factor(TreatTimes[idx_all])

myColData <- data.frame(conds = conds, batch = batch, days = days, dur = dur)
rownames(myColData) <-  SamNames[idx_all]

dds <- DESeqDataSetFromMatrix(countData = data_mat_oi, colData = myColData, design =  ~ conds)

## be sure that you know your reference 
#dds$condition <- relevel(dds$condition, ref="DMSO")

##############################################################################################
########################## transformed values
# (Viera) normalization - I am using rlog function
# there is also another normalization type called "vst" in DEseq2

myRLOG <- rlog(dds)

# with assay, the matrix is extracted
rlogMat <- assay(myRLOG)
#head(rlogMat)

### tests ###
# R-R plots (replicates)
which_rep <- c(8,9) # set
sample1 <- rownames(myColData)[which_rep[1]]
sample2 <- rownames(myColData)[which_rep[2]]
# the same sample which differs in batch
par(mfrow=c(2,1), mar=c(2.5,2.5,2,0.5), mgp=c(1.5,0.5,0), cex.lab=0.75, cex.axis=0.75, cex.main=0.85)
plot(data_mat_oi[,which_rep[1]], data_mat_oi[,which_rep[2]], xlab=sample1, ylab=sample2, main="raw counts")
r1 <- cor.test(data_mat_oi[,which_rep[1]], data_mat_oi[,which_rep[2]])
text(x=10000, y=100000, as.character(r1$estimate))
r2 <- cor.test(rlogMat[,which_rep[1]], rlogMat[,which_rep[2]])
plot(rlogMat[,which_rep[1]], rlogMat[,which_rep[2]],  xlab=sample1, ylab=sample2, main="normalized values")
text(x=1, y=10, as.character(r2$estimate))

# heatmaps (dendrograms)
pheatmap(data_mat_oi, cluster_rows = FALSE, labels_col = SamNames[idx_all], show_rownames = FALSE,
         cellheight = 0, legend = FALSE, main="raw counts")
# even if we use raw counts, we see that biological replicates cluster together
# clustering does not follow conditions

pheatmap(rlogMat, cluster_rows = FALSE, labels_col = SamNames[idx_all], show_rownames = FALSE,
         cellheight = 0, legend = FALSE, main="normalized values")
# clustering improved after rlog normalization but stil follow batch numbers
# png file: Heatmap_normalizedValues_39samples

# the plotPCA function is implemented within DESeq2 pckg. and uses spectral decomposition 
par(mfrow=c(2,2))
plotPCA(myRLOG, intgroup=c("days")) + ggtitle("PCA - rlog normalized counts - per day")
plotPCA(myRLOG, intgroup=c("conds")) + ggtitle("PCA - rlog normalized counts - per condition")
plotPCA(myRLOG, intgroup=c("dur")) + ggtitle("PCA - rlog normalized counts - per Time of treatment")


pcaData <- plotPCA(myRLOG, intgroup=c("conds", "dur", "days"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conds, shape = days, size = dur)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

### try to plot PCA plot that include all infomtion
rv <- rowVars(assay(myRLOG))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds)[select,]))

PCcollect <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2], PC3=pca$x[,3],PC4=pca$x[,4],PC5=pca$x[,5],
                        conds = conds, days = days, dur = dur)
PCcollect$dur <- as.factor(PCcollect$dur)
ggplot(PCcollect, aes(x=PC3, y=PC4, color=conds, shape = days, size = dur)) + geom_point()
