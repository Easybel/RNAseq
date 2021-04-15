  #topGO muss geladen werden damit wir es später benutzen können.
library(topGO)

# read in the gene 2 GO mapping that is done with the Nmen and FA1090 species, as Ngo is not annotated
geneID2GO <- readMappings(file ="/home/isabel/Dokumente/P5_ExpEvol_Ngo/dictionaries/Nmen_FA1090_GO_sorted.txt")
str(head(geneID2GO))


# here all gene names are saved -- this is the gene universe


# read in the diff. expression data
DiffExprResults <- read.delim("/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/output/ExpDataNgo_wIDsFromBlast_onlyGOHit.txt",sep = "\t", header = TRUE)

geneNames <- as.character(DiffExprResults$info2GO)
head(geneNames)

myInterestingGenes <- as.character(DiffExprResults$info2GO[DiffExprResults$l2f >= 0.5 & DiffExprResults$pvalue <0.05])
HowmanyIntGenes <- sum(as.integer(geneNames %in% myInterestingGenes))

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)


#Mit den topGO-Funktionen können wir jetzt alle GO-Terme testen. 
ergebnisseClassic <- lapply(c("BP","MF","CC"),FUN=function(ont){ 
  GOdata <- new("topGOdata",
                ontology = ont,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO,
                nodeSize=10) #Nur Terme mit mindestens 10 Genen (egal ob sie signifikant sind) werden geteste
  
  
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(object=GOdata, pvalue = resultFisher,topNodes=length(resultFisher@score))
allRes$pvalue <- resultFisher@score[allRes$GO.ID] #Leider wurden die p-Werte formatiert. Wir ersetzen sie mit den echten Zahlen.
rownames(allRes) <- allRes[,1] #Wir benutzen die GO.IDs als Zeilennamen.
return(allRes)
})

head(ergebnisseClassic[[1]])

# save the data
write.table(ergebnisseClassic[[1]], file = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/Ontology/BP_absl2fhigher1_pvaluelower0.05.csv", append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
write.table(ergebnisseClassic[[2]], file = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/Ontology/MF_absl2fhigher1_pvaluelower0.05.csv", append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
write.table(ergebnisseClassic[[3]], file = "/home/isabel/Dokumente/P5_ExpEvol_Ngo/RNA/Ontology/CC_absl2fhigher1_pvaluelower0.05.csv", append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
