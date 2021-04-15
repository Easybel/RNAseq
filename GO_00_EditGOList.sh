# this script edit the output from the geneontology.org page and gives a slimer file out

#GoPath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/Nmeningitidis_MC58/"
#name="Nmen_GO"
GoPath="/home/isabel/Dokumente/ExpEvol/Ngo/dictionaries/FA1090/"
name="FA1090_GO"

awk -v OFS='\t' '{split($0,data,"\t"); print(data[2],data[3],data[3],data[10],data[13])}' $GoPath/$name".txt" > $GoPath/$name"_editted.txt"


