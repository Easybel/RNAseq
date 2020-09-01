# Here, the fasta will be indexed for the Star mapping.
# Define the paths to the dictionary (this is also the output path) and its name. 
# Also provide the path to STAR and 

myDictPath="/projects/ag-advol/dictionaries/MS11_ncbi_IR" #!

# name of dictionary
dict="MS11"
ann="MS11"

#Define folders where software is installed
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64"

cd $StarFold
./STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $myDictPath --genomeFastaFiles $myDictPath/$dict".fasta" \
--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile $myDictPath/$ann".gff3" --sjdbOverhang 149
