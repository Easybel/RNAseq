#Ist topGO installiert? 
# Wenn nicht, wird es jetzt automatisch ueber das bioConductor Netzwerk installiert.
if(!"topGO"%in%rownames(installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
  BiocManager::install("topGO")
}
#topGO muss geladen werden damit wir es später benutzen können.
library(topGO)
#Jetzt laden wir die Daten aus topGODaten.RData. Die Datei sollte im Arbeitsverzeichnis liegen.
getwd() #Hier sind wir gerade.
list.files() #Hier sollte topGODaten.RData dabei sein. Wenn nicht, wechseln Sie den Ordner mit setwd().
load("topGODaten.RData")

#Jetzt bereiten wir ein Objekt vor das die signifikanten Gene so präsentiert, wie topGO es fordert.
trefferVektor <- rep(0,length(alleGene)) #Wir schreiben soviele Nullen, wie es Gene gibt.
names(trefferVektor) <- alleGene
trefferVektor[signifikanteGene] <- 1 #Signifikante Gene werden mit einer 1 gekennzeichnet.
trefferVektor <- as.factor(trefferVektor) #topGO funktioniert nur mit faktoriellen Trefferlisten.
head(trefferVektor)

#Mit den topGO-Funktionen können wir jetzt alle GO-Terme testen. 
#Dabei benutzen wir lapply. Das ist eine Funktion mit der auf einen Input-Vektor oder eine Inputliste
#eine weitere Funktion anwenden kann. Die Funktion definieren wir hier bei FUN=... .Wir wenden sie auf
#den Vektor c("BP","MF","CC") an. Innerhalb der Funktion heisst der Input dann ont.
ergebnisseClassic <- lapply(c("BP","MF","CC"),FUN=function(ont){ 
  GOdata <- new("topGOdata",
                ontology = ont,
                allGenes = trefferVektor,
                annot = annFUN.gene2GO,
                gene2GO = ontologien[[ont]],
                nodeSize=10) #Nur Terme mit mindestens 10 Genen (egal ob sie signifikant sind) werden getestet.
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(object=GOdata, pvalue = resultFisher,topNodes=length(resultFisher@score))
  allRes$pvalue <- resultFisher@score[allRes$GO.ID] #Leider wurden die p-Werte formatiert. Wir ersetzen sie mit den echten Zahlen.
  rownames(allRes) <- allRes[,1] #Wir benutzen die GO.IDs als Zeilennamen.
  return(allRes)
})
head(ergebnisseClassic[[1]])
head(ergebnisseClassic[[2]])
head(ergebnisseClassic[[3]])

ergebnisseClassic[[1]]["GO:0016052",] #Was ist der p-Wert für "carbohydrate catabolic process"?
fisher.test(rbind(c(10,76),c(99,5489)),alternative="greater") #Das ist der gleiche p-Wert.

#Mit dem elim-Algorithmus sieht das etwas anders aus.
ergebnisseElim <- lapply(c("BP","MF","CC"),FUN=function(ont){ 
  GOdata <- new("topGOdata",
                ontology = ont,
                allGenes = trefferVektor,
                annot = annFUN.gene2GO,
                gene2GO = ontologien[[ont]],
                nodeSize=10) #Nur Terme mit mindestens 10 Genen (egal ob sie signifikant sind) werden getestet.
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher") #Wir benutzen jetzt den elim-Algorithmus.
  allRes <- GenTable(object=GOdata, pvalue = resultFisher,topNodes=length(resultFisher@score))
  allRes$pvalue <- resultFisher@score[allRes$GO.ID] #Leider wurden die p-Werte formatiert. Wir ersetzen sie mit den echten Zahlen.
  rownames(allRes) <- allRes[,1] #Wir benutzen die GO.IDs als Zeilennamen.
  return(allRes)
})
head(ergebnisseElim[[1]])
head(ergebnisseElim[[2]])
head(ergebnisseElim[[3]])

ergebnisseElim[[1]]["GO:0016052",] #Der p-Wert für "carbohydrate catabolic process" ist jetzt grösser.
ergebnisseClassic[[1]]["GO:0016052",]
