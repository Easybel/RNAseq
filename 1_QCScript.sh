#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH -N 1
#SBATCH --time=05:00:00
#SBATCH --array=1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user irathman@uni-koeln.de 
i=$SLURM_ARRAY_TASK_ID

#Define input and output paths
#Make directories if necessary
myDataRaw="/projects/ag-advol/RNAseq/Raw"
myQCRaw="/projects/ag-advol/RNAseq/Raw/QC"
myDataTrim="/projects/ag-advol/RNAseq/Trim"
myQCTrim="/projects/ag-advol/RNAseq/Trim/QC"

#Define folders where software is installed
FastQCFold="/home/irathman/sw/FastQC"
TrimmFold="/home/irathman/sw/Trimmomatic-0.36"

ID=$(ls -1 $myDataRaw | grep "R2_001.fastq.gz" | sed -n ''$i'p' | cut -d"_" -f1,2,3,4)

#FastQC before
cd $FastQCFold
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_R1_001.fastq.gz" -t 8
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_R2_001.fastq.gz" -t 8

#Trimm_Script.sh
IDout=$(echo $ID | cut -d"_" -f1,2)
cd $TrimmFold
java -Xmx30G -Xms24G -jar trimmomatic-0.36.jar PE -threads 8 -trimlog $ID.TrimLog $myDataRaw/$ID"_R1_001.fastq.gz" $myDataRaw/$ID"_R2_001.fastq.gz" $myDataTrim/$IDout"_1P".fastq $myDataTrim/$IDout"_1U".fastq $myDataTrim/$IDout"_2P".fastq $myDataTrim/$IDout"_2U".fastq ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#FastQC_script_after.sh
cd $FastQCFold
./fastqc -o $myQCTrim/ $myDataTrim/$IDout"_1P.fastq" -t 8
./fastqc -o $myQCTrim/ $myDataTrim/$IDout"_2P.fastq" -t 8

exit 0
