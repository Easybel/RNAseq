####
### this script takes a fasta and blasts it against a db
####
cd /home/isabel/Documents/Doktorarbeit_Mai2022/software/ncbi-blast-2.13.0+/bin/
# the paths where to find the data
mydbPath="/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/1B_MS11_ncbi"
myGenePath="/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/3_FA1090/ManualLists"
myOutPath="/home/isabel/Documents/Doktorarbeit_Mai2022/P5_Ngo_fromSciebo/dictionaries/3_FA1090/ManualLists/Blast_Phage_2_MS11"

# set the variables 
# The genes that you would like to blast ...
genes="FA1090_singleGenes_nucl"

# ...against this dictionary 
dict="MS11"

# the outputs will have this name:
outName="Blast_"$genes"_2_"$dict

## OPTIONAL ##
### if you wish to use a fasta as reference db, then run this command before to create the reference
#makeblastdb -in $mydbPath/$dict.fasta -dbtype nucl -parse_seqids

# blast the single genes
./blastn -db $mydbPath/$dict.fasta -query $myGenePath/$genes.fasta -outfmt 7 -out $myOutPath/$outName".tab"

# get only the hit lines
grep -v "^#" $myOutPath/$outName".tab" > $myOutPath/$outName"_editted.tab"
