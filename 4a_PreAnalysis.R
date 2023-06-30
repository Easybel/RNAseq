## Neisseria gonorrhoeae transcriptome project

# This script takes the count table and analyses if everything is ok with the samples
# -- do they cluster biologically?
# -- Do we need batch correction?

# input: count table
#  ----> this is the output from featureCounts
#  ----> this table only contains the CDS counts; 
#       for example: gene ID, NgoG4_1, _2, _3, ...

# we have 7 different phenotypes: 
# -- we want to compare the planktonic against biofilm **
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

# do you want to exclude multi mapper?
excMM = "YES"

################ LOAD DATA #######################
if (excMM == "YES") {
List_MM_file <- read.csv("/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/MULTIMAPPER_Ms11toMs11/MS11_2_MS11_allGenes_woRNA.csv",
                        sep="\t", header = TRUE)
List_MM      <-  paste("gene-",List_MM_file[,1],sep="")
exc_pilS_E = paste("gene-",c("NGFG_01821","NGFG_02431","NGFG_02484","NGFG_02482","NGFG_02405","NGFG_02481","NGFG_02253","NGFG_00014","NGFG_02485","NGFG_02487","NGFG_01819","NGFG_01818"),sep="")

List_MM <- c(List_MM,exc_pilS_E)}
####
## Ngo -- here the data is loaded
data_raw <- read.table(inName, sep=" ", header = TRUE)
ID2geneInfo <- read.csv("/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1A_MS11_ManualLists/NgoBlast2_NmenANDFA1090_onlyCDS.csv",
                        sep="\t", header = TRUE)
data_samnum <- length(data_raw) -1

############################# EXTRA INFORMATION ################################
# information about the data
Conditions <- c(rep("G4wt",3),rep("G4pilEKO",3),rep("G4pilENG17",3),rep("G4pilENG24",3),rep("G4pilENG32",3))   #set
replicate <- rep(c(1,2,3),data_samnum/3)                                                                       #set

# automatically generates the Sample Names
SamNames <- paste(Conditions,replicate,sep = "_")

# which replicates and conditions (here: genotypes) are you interested in? 
# if more than 1, add them with |
cond_oi <- c("G4pilEKO|G4pilENG17|G4pilENG24|G4pilENG32") 
ctl_oi <- c("G4wt")
replicate_oi <- c("1|2|3")

# search pattern
idx_cond_oi <- grep(cond_oi,Conditions)
idx_replicate_oi <- grep(replicate_oi,replicate)
idx_ctl_oi <- grep(ctl_oi,Conditions)

idx_ctl     <- intersect(idx_ctl_oi,idx_replicate_oi)
idx_samples <- intersect(idx_cond_oi,idx_replicate_oi)
idx_all     <- sort(c(idx_ctl,idx_samples))

#######################################################################################
############### check the number of counts per sample -- analyse if there are connections

counts_per_sample <- apply(data_raw[,(1:data_samnum)+1],2,sum)
counts <- data.frame(counts_per_sample,Conditions,replicate)
counts$replicate <- as.factor(counts$replicate)
counts$Conditions <- as.factor(counts$Conditions)

# Basic box plot
counts_depCond <-ggplot(counts, aes(x=Conditions, y=counts_per_sample, fill = Conditions)) +
  geom_boxplot(notch=FALSE)
counts_depCond

#############################################################################################################
################## Process the data

# Pre-filter: optionally, data is pre-filtered based on distributions
data_filt <- data_raw #[which(tab$med > 9 & (tab$sum > 15 | tab$sum < 10000)),]

# convert data in table to matrix
data_mat <- as.matrix(data_filt[,(1:data_samnum)+1]) 

rownames(data_mat) <- data_filt[,1] 
colnames(data_mat) <- SamNames

#exclude the multimapper
if (excMM == "YES") {
idx_exc <- which(data_filt[,1] %in% List_MM)
data_mat <- data_mat[-(idx_exc),]}

################## Create the dds object for the analysis
# give information on samples, conditions, batch

data_mat_oi <- data_mat[,idx_all]

conds <- as.factor(Conditions[idx_all])
replicate <- as.factor(replicate[idx_all])

myColData <- data.frame(conds = conds, replicate = replicate)
rownames(myColData) <-  SamNames[idx_all]

dds <- DESeqDataSetFromMatrix(countData = data_mat_oi, colData = myColData, design = ~ replicate + conds)

## be sure that you know your reference 
# dds$condition <- relevel(dds$condition, ref="biofilm")

##############################################################################################
########################## transformed values
# (Viera) normalization - I am using rlog function
# there is also another normalization type called "vst" in DEseq2

myRLOG <- rlog(dds,blind=TRUE)

# with assay, the matrix is extracted
rlogMat <- assay(myRLOG)
rlogMat_cent <- rlogMat - apply(rlogMat,1,mean)
#head(rlogMat)

# heatmaps (dendrograms)
pheatmap(data_mat_oi, cluster_rows = FALSE, labels_col = SamNames[idx_all], show_rownames = FALSE,
         cellheight = 0, legend = FALSE, main="raw counts")
# clustering improved after rlog normalization - fits to sample type!
pheatmap(rlogMat_cent, cluster_rows = FALSE, labels_col = SamNames[idx_all], show_rownames = FALSE,
         cellheight = 0, legend = FALSE, main="normalized values")

# the plotPCA function is implemented within DESeq2 pckg. and uses spectral decomposition 
par(mfrow=c(2,2))
plotPCA(myRLOG, intgroup=c("replicate"),ntop=1900) + ggtitle("PCA - rlog normalized counts - per day")
plotPCA(myRLOG, intgroup=c("conds"),ntop=1900) + ggtitle("PCA - rlog normalized counts - per condition")

pcaData <- plotPCA(myRLOG, intgroup=c("conds", "replicate"), returnData=TRUE) #, 
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conds)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

### try to plot PCA plot that include all information
rv <- rowVars(assay(myRLOG))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds)[select,]))

PCcollect <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2], PC3=pca$x[,3],PC4=pca$x[,4],PC5=pca$x[,5],
                        conds = conds, replicate = replicate)
ggplot(PCcollect, aes(x=PC3, y=PC4, color=conds, shape = replicate)) + 
  geom_point(size=3)
