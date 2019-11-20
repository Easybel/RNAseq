#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --array=1-11

i=$SLURM_ARRAY_TASK_ID

#Define input and output paths
#Make directories if necessary

date="190319"
myDataRaw="/projects/ag-advol/melih/raw3"

mkdir /scratch/turkmc/FastQC_outputs_before_$date
FastQC="/scratch/turkmc/FastQC_outputs_before_$date" 

mkdir /scratch/turkmc/Trimm_results_$date
myDir="/scratch/turkmc/Trimm_results_$date"

mkdir /scratch/turkmc/Star_$date
myDataPath="/scratch/turkmc/Star_$date"

mkdir /scratch/turkmc/FastQC_outputs_after_$date
FastQC_afterTrimm="/scratch/turkmc/FastQC_outputs_after_$date" 

mkdir /$myDataPath/csvOUTPUTS_$date
csvOUTPUTS="$myDataPath/csvOUTPUTS_$date"

#Path of dictionary, put your dictionaries in that folder
dictPath="/projects/ag-advol/melih/Dict/BsubNC_000964wt/"

#Define folders where software is installed
TrimmFold="/home/myuekse0/sw/Trimmomatic-0.36"
FastQCFold="/home/myuekse0/sw/FastQC"
bwaFold="/home/myuekse0/sw/bwa.kit"
samFold="/home/myuekse0/sw/samtools-1.3.1"
picardFold="/home/myuekse0/sw/picard-2.6.0"
GATKFold="/home/myuekse0/sw/GATK3.5"
snpEffFold="/home/myuekse0/sw/snpEff"
bedtoolsFold="/home/myuekse0/sw/bedtools2/bin"
StarFold="/home/myuekse0/sw/STAR-2.6.0a/bin/Linux_x86_64/"
featureCounts="/home/myuekse0/sw/subread-1.6.2-source/bin/"

#for ((i=1; i<7; i++))
#do
#FastQC.sh
#cd $FastQCFold
#ID=$(ls -1 $myDataRaw | grep "N0" | grep "_1.fq" | sed -n ''$i'p' | cut -d"_" -f1)
#./fastqc -o $FastQC/ $myDataRaw/$ID"_1.fq" -t 8
#./fastqc -o $FastQC/ $myDataRaw/$ID"_2.fq" -t 8

#Trimm_Script.sh
cd $TrimmFold
ID=$(ls -1 $myDataRaw | grep "0" | grep "_1.fq.gz" | sed -n ''$i'p' | cut -d"_" -f1)
java -Xmx30G -Xms24G -jar trimmomatic-0.36.jar PE -threads 8 -trimlog $ID.TrimLog $myDataRaw/$ID"_1.fq.gz" $myDataRaw/$ID"_2.fq.gz" $myDir/$ID"_1P".fq $myDir/$ID"_1U".fq $myDir/$ID"_2P".fq $myDir/$ID"_2U".fq ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#FastQC_script_after.sh
cd $FastQCFold
ID=$(ls -1 $myDir | grep "0" | grep "_1P.fq" | sed -n ''$i'p' | cut -d"_" -f1)
./fastqc -o $FastQC_afterTrimm/ $myDir/$ID"_1P.fq" -t 8
./fastqc -o $FastQC_afterTrimm/ $myDir/$ID"_2P.fq" -t 8

cd $StarFold
ID=$(ls -1 $myDir | grep "0" | grep "_1P.fq" | sed -n ''$i'p' | cut -d"_" -f1)
./STAR --runThreadN 8 --runMode genomeGenerate --genomeSAindexNbases 4 --genomeDir $dictPath --genomeFastaFiles $dictPath/BsubNC_000964wt.fasta
chmod a+x $dictPath/*
./STAR --runThreadN 8 --genomeDir $dictPath/ --outSAMmultNmax 10 --outFileNamePrefix $myDataPath/$ID"_" --readFilesIn $myDir/$ID"_1P".fq $myDir/$ID"_2P".fq

cd $featureCounts
./featureCounts -p -T 8 -F SAF -a $dictPath/BsubNC_000964wt.saf -o $myDataPath/$ID"_raw.count" $myDataPath/$ID"_Aligned.out".sam
#done
exit 0