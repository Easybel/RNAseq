
### this script takes a fasta and blasts it against a db
### if you wish to use a fasta as reference db, then run this command before

#  makeblastdb -in $mydbPath/$dict.fasta -dbtype prot -parse_seqids

# set the variable names
#mydbPath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/Nmeningitidis_MC58/Nmeningitidis_MC58"
mydbPath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/FA1090"
myGenePath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/MS11_ncbi_IR/Ann_Ngo2Nmen/singleGenes"
myOutPath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/MS11_ncbi_IR/Ann_Ngo2Nmen/singleGenesBlast2Nmen"

dict="FA1090"

# blast the single genes
blastn -db $mydbPath/$dict.fasta -query $myGenePath/MS11_singleGenes.fasta -outfmt 7 -out $myOutPath/"Blast_NgonGenes_to_"$dict".tab"

# get only the hit lines
grep -v "^#" $myOutPath/"Blast_NgonGenes_to_"$dict".tab" > $myOutPath/"Blast_NgonGenes_to_"$dict"_editted.tab"
