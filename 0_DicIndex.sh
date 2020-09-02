# Here, the fasta will be indexed for the Star mapping.
# Define the paths to the dictionary (this is also the output path) and its name. 
# Also provide the path to STAR and 

myDictPath="/projects/ag-advol/dictionaries/MS11_ncbi_IR" #!

# name of dictionary
dict="MS11"  #!
ann="MS11" #!

#Define folders where software is installed
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64" #!

#The basic options to generate genome indices with STAR are as follows:
#--runThreadN NumberOfThreads
#--runMode genomeGenerate
#--genomeDir /path/to/genomeDir
#--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
#--sjdbGTFfile /path/to/annotations.gtf
#--sjdbOverhang ReadLength-1
# !!--genomeSAindexNbases min(14, log2(GenomeLength)/2 - 1) (For small genomes, this parameter must be scaled down)!!

cd $StarFold
./STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $myDictPath --genomeFastaFiles $myDictPath/$dict".fasta" \
--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile $myDictPath/$ann".gff3" --sjdbOverhang 149 
