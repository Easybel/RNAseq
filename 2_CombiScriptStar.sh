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
#Define input and output paths
#Make directories if necessary
myDataTrim="/projects/ag-advol/RNAseq/Trim"
myDictPath="/projects/ag-advol/dictionaries/MS11_ncbi_IR"
myDataPath="/scratch/easy/Ngo/RNA"

#Define folders where software is installed
samFold="/home/irathman/sw/samtools-1.8"
StarFold="/home/irathman/sw/STAR-2.5.3a/bin/Linux_x86_64/"
featureCounts="/home/irathman/sw/subread-2.0.1-source/bin"

# here you map against: .fasta
dict="MS11"
ID=$(ls -1 $myDataTrim | grep "_1P.fastq" | sed -n ''$i'p' | cut -d"_" -f1,2,3,4)
IDout=$(echo $ID)

####### Mapping the reads to the reference with STAR
# basic options:
#--runThreadN NumberOfThreads
#--genomeDir /path/to/genomeDir
#--readFilesIn /path/to/read1 [/path/to/read2]
#--outFileNamePrefix /path/to/output/dir/prefix
# Multimappers
# --outSAMmultNmax limits  the  number  of  output  alignments  (SAM  lines)  for multimappers ??
# Chimeric/ circular reads
# --chimSegmentMin positive value ??how to set??
# --chimOutType WithinBAM
cd $StarFold
./STAR --runThreadN 8 --genomeDir $myDictPath/ --outSAMmultNmax 10 --outFileNamePrefix $myDataPath/$IDout"_" \ 
--readFilesIn $myDataTrim/$ID"_1P.fastq" $myDataTrim/$ID"_2P.fastq" 
echo "done with mapping"

##### data cpnverted
#..sam is converted here to bam, sorted and indexed, so that it can be oppened with IGVviewer
cd $samFold
./samtools view -b $myDataPath/$IDout"_Aligned.out.sam" --threads 8 -T $myDictPath/$dict.fasta -o $myDataPath/$IDout".bam"
./samtools sort $myDataPath/$IDout".bam" --threads 8 --reference $myDictPath/$dict.fasta -o $myDataPath/$IDout"_sort.bam"
./samtools index -b $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_sort.bam.bai"

###### FeatureCounts: Read quantification
# Here the featurecount is done with subread package -> output: count table
# Input data (i) aligned reads in sam/bam, (ii) list of genomic features in either GTF, GFF or simplified annotation format (SAF)
# Options:
# -t which feature Type to count, e.g. "exon", "CDS"
# -p : for paired reads and -f for counting fragments instead of reads
# -M count multimappers each as 1 count or with --fraction 
# -O count multi-overlapping reads each as 1 count or with --fraction; --minOverlap give the number of bp that overlap has to be long

cd $featureCounts
# change the output name depending on settings
IDout2=$IDout"_O_minO20"
./featureCounts -p -O --minOverlap 20 -t "CDS" -T 8 -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout2"_CDS_raw.count"
./featureCounts -p -O --minOverlap 20 -t "exon" -T 8 -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout2"_exon_raw.count"
#done
exit 0
