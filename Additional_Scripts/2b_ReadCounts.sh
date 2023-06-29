#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH -N 1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user irathman@uni-koeln.de 
#SBATCH --time=05:00:00
#SBATCH --array=1-6
i=$SLURM_ARRAY_TASK_ID

####### This script gives the output .. _stats.txt that contains information on how many alignments in the bam file 
# fall into which category
#######

#Define input and output paths
#Make directories if necessary
myDataTrim="/projects/ag-advol/RNAseq/Trim"
myDictPath="/projects/ag-advol/dictionaries/MS11_ncbi_IR"
myDataPath="/scratch/easy/Ngo/RNA"

#Define folders where software is installed
samFold="/home/irathman/sw/samtools-1.8"
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64"
featureCounts="/home/irathman/sw/subread-2.0.1-source/bin"

# here you map against: .fasta
dict="MS11"
ID=$(ls -1 $myDataTrim | grep "Ngo" | grep "_1P.fastq" | sed -n ''$i'p' | cut -d"_" -f1,2,3)
IDout=$(echo $ID)

###########

touch $myDataPath/$IDout"_stats.txt"
cd $samFold
echo "This script tells you how many alignments in the bam file fall into which category" >  $myDataPath/$IDout"_stats.txt"
echo "1.count all reads with samtools view -c" >>  $myDataPath/$IDout"_stats.txt"
./samtools view -c $myDataPath/$IDout"_Aligned.out.sam" >> $myDataPath/$IDout"_stats.txt"

echo "2. count all primary aligned reads with samtools view -c -F 260" >> $myDataPath/$IDout"_stats.txt"
./samtools view -c -F 260 $myDataPath/$IDout"_Aligned.out.sam" >> $myDataPath/$IDout"_stats.txt"

echo "count all unmapped reads with samtools view -c -f 4" >> $myDataPath/$IDout"_stats.txt"
./samtools view -c -f 4 $myDataPath/$IDout"_Aligned.out.sam" >> $myDataPath/$IDout"_stats.txt"

echo "count all not primary alignments with samtools view -c -f 256" >>  $myDataPath/$IDout"_stats.txt"
./samtools view -c -f 256 $myDataPath/$IDout"_Aligned.out.sam" >>  $myDataPath/$IDout"_stats.txt"

echo "count all reads paired and mapped in proper pairs with samtools view -c -f 3" >>  $myDataPath/$IDout"_stats.txt"
./samtools view -c -f 3 $myDataPath/$IDout"_Aligned.out.sam" >>  $myDataPath/$IDout"_stats.txt"

echo "count only the uniquely mapped ones with the help of the NH flag" >> $myDataPath/$IDout"_stats.txt"
grep -v "^@" $myDataPath/$IDout"_Aligned.out.sam" | awk '{split($12,multi,":"); if(multi[3]==1) print $0}' | wc -l  >> $myDataPath/$IDout"_stats.txt"


exit 0
