#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -N 1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user irathman@uni-koeln.de 
#SBATCH --time=05:00:00
#SBATCH --array=1-6
i=$SLURM_ARRAY_TASK_ID

#Define input and output paths
#Make directories if necessary

myDataPath="/scratch/easy/NgoTrans_Run3/Results"

#Define folders where software is installed
samFold="/home/irathman/sw/samtools-1.8"
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64/"
bedFold="/home/irathman/sw/bedtools2/bin" #!

# here you map against: .fasta

ID=$(ls -1 $myDataTrim | grep "_1P.fastq" | sed -n ''$i'p' | cut -d"_" -f1,2,3)
IDout=$(echo $ID)

#Coverage
module unload gnu/4.4.7
module load gnu/5.1.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_coverage.txt"


exit 0

