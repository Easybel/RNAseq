
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
myRawData="/projects/ag-advol/Isabel/RNARawReads"

myTrimData="/scratch2/easy/DirectsRNA/TrimmResults" 

myQC="/scratch2/easy/DirectsRNA/QCResults"

myDataPath="/scratch2/easy/DirectsRNA"
#mkdir /scratch2/easy/DirectsRNA/TrimmResults/FastQCAfterTrimm

#Path of dictionary, put your dictionaries in that folder
myDictPath="/projects/ag-advol/dictionary/BsubNC_000964wt/"

#Define folders where software is installed
TrimmFold="/home/irathman/sw/Trimmomatic-0.36"
FastQCFold="/home/irathman/sw/FastQC"
bwaFold="/home/irathman/sw/bwa-0.7.17"
samFold="/home/irathman/sw/samtools-1.8"
picardFold="/home/irathman/sw"
GATKFold="/home/irathman/sw/GATK3.5"
snpEffFold="/home/irathman/sw/snpEff"
bedtoolsFold="/home/irathman/sw/bedtools2/bin"
StarFold=
featureCounts=


ID=$(ls -1 $myRawData | grep "N0" | grep "_1.fq" | sed -n ''$i'p' | cut -d"_" -f1)
IDout=$ID""

cd $FastQCFold
./fastqc -o $myQC/ $myDataRaw/$ID"_1.fq" -t 8
./fastqc -o $myQC/ $myDataRaw/$ID"_2.fq" -t 8

#Trimm_Script.sh
cd $TrimmFold
ID=$(ls -1 $myDataRaw | grep "0" | grep "_1.fq.gz" | sed -n ''$i'p' | cut -d"_" -f1)
java -Xmx30G -Xms24G -jar trimmomatic-0.36.jar PE -threads 8 -trimlog $ID.TrimLog $myDataRaw/$ID"_1.fq.gz" $myDataRaw/$ID"_2.fq.gz" $myTrimData/$ID"_1P".fq $myTrimData/$ID"_1U".fq $myTrimData/$ID"_2P".fq $myTrimData/$ID"_2U".fq ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#FastQC_script_after.sh
cd $FastQCFold
ID=$(ls -1 $myDir | grep "0" | grep "_1P.fq" | sed -n ''$i'p' | cut -d"_" -f1)
./fastqc -o $myQC/ $myTrimData/$ID"_1P.fq" -t 8
./fastqc -o $myQC/ $myTrimData/$ID"_2P.fq" -t 8

cd $StarFold
ID=$(ls -1 $myDir | grep "0" | grep "_1P.fq" | sed -n ''$i'p' | cut -d"_" -f1)
./STAR --runThreadN 8 --runMode genomeGenerate --genomeSAindexNbases 4 --genomeDir $myDictPath --genomeFastaFiles $myDictPath/BsubNC_000964wt.fasta
chmod a+x $dictPath/*
./STAR --runThreadN 8 --genomeDir $myDictPath/ --outSAMmultNmax 10 --outFileNamePrefix $myDataPath/$ID"_" --readFilesIn $myTrimData/$ID"_1P".fq $myTrimData/$ID"_2P".fq

cd $featureCounts
./featureCounts -p -T 8 -F SAF -a $myDictPath/BsubNC_000964wt.saf -o $myDataPath/$ID"_raw.count" $myDataPath/$ID"_Aligned.out".sam
#done
exit 0
