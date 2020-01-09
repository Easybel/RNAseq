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

#Define input and output paths
#Make directories if necessary
myDataTrim="/projects/ag-advol/RNAseq/Trim"
myDictPath="/projects/ag-advol/dictionaries/BsubNC_000964wt"
myDataPath="/scratch2/easy/Directs/DirectsRNA"

#Define folders where software is installed
FastQCFold="/home/irathman/sw/FastQC"
TrimmFold="/home/irathman/sw/Trimmomatic-0.36"
bwaFold="/home/irathman/sw/bwa-0.7.17"
samFold="/home/irathman/sw/samtools-1.8"
picardFold="/home/irathman/sw"
GATKFold="/home/irathman/sw/GATK3.5"
bcfFold="/home/irathman/sw/bcftools-1.8"
snpEffFold="/home/irathman/sw/snpEff"
bedtoolsFold="/home/irathman/sw/bedtools2/bin"
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64/"
featureCounts="/home/irathman/sw/subread-1.6.5-source/bin"

# here you map against: .fasta
dict="BsubNC_000964wt"
ID=$(ls -1 $myDataRaw | grep "W" | grep "_1.fastq.gz" | sed -n ''$i'p' | cut -d"_" -f1)
IDout=$ID"_lm2Bsub"

cd $StarFold
ID=$(ls -1 $myDir | grep "0" | grep "_1P.fq" | sed -n ''$i'p' | cut -d"_" -f1)
./STAR --runThreadN 8 --runMode genomeGenerate --genomeSAindexNbases 4 --genomeDir $myDictPath --genomeFastaFiles $myDictPath/$dict
chmod a+x $myDictPath/*
./STAR --runThreadN 8 --genomeDir $myDictPath/ --outSAMmultNmax 10 --outFileNamePrefix $myDataPath/$ID"_" --readFilesIn $myDataRaw/$ID"_1P".fastq $myDataRaw/$ID"_2P".fastq

cd $featureCounts
./featureCounts -p -T 8 -F SAF -a $myDictPath/BsubNC_000964wt.saf -o $myDataPath/$ID"_raw.count" $myDataPath/$ID"_Aligned.out".sam
#done
exit 0
