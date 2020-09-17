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


# Read quantification
# Here the featurecount is done with subread package -> output: count table
# Input data (i) aligned reads in sam/bam, (ii) list of genomic features in either GTF, GFF or simplified annotation format (SAF)

cd $featureCounts

# default setting: this is what I assume to be default settings
IDout1=$IDout
./featureCounts -p -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout1"_CDS_raw.count"
./featureCounts -p -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout1"_exon_raw.count"

# fragment option: does it make a difference with the default settings
IDout2=$IDout"_f"
./featureCounts -p -f -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout2"_CDS_raw.count"
./featureCounts -p -f -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout2"_exon_raw.count"

# count multimapper 
IDout3=$IDout"_M"
./featureCounts -p -M -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout3"_CDS_raw.count"
./featureCounts -p -M -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \

# count multimapper as fraction 
IDout4=$IDout"_Mfrac"
./featureCounts -p -M --fraction -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout4"_CDS_raw.count"
./featureCounts -p -M --fraction -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout4"_exon_raw.count"

# all settings also with multimap fraction, overlap and overlap fraction
IDout5=$IDout"_MfracOfrac"
./featureCounts -p -M --fraction -O --fraction -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout5"_CDS_raw.count"
./featureCounts -p -M --fraction -O --fraction -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout5"_exon_raw.count"

# all settings also with multimap fraction, overlap and overlap fraction
IDout5=$IDout"_MfracOfrac_minO10"
./featureCounts -p -M --fraction -O --fraction --minOverlap 10 -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout5"_CDS_raw.count"
./featureCounts -p -M --fraction -O --fraction --minOverlap 10 -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout5"_exon_raw.count"

# all settings also with multimap fraction, overlap and overlap fraction
IDout6=$IDout"_MfracOfrac_minO5"
./featureCounts -p -M --fraction -O --fraction --minOverlap 5 -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout6"_CDS_raw.count"
./featureCounts -p -M --fraction -O --fraction --minOverlap 5 -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout6"_exon_raw.count"


# all settings also with multimap fraction, overlap and overlap fraction
IDout7=$IDout"_MfracOfrac_minO20"
./featureCounts -p -M --fraction -O --fraction --minOverlap 20 -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout7"_CDS_raw.count"
./featureCounts -p -M --fraction -O --fraction --minOverlap 20 -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout7"_exon_raw.count"


# all settings also with multimap fraction, overlap and overlap fraction
IDout8=$IDout"_O_minO20"
./featureCounts -p -O --minOverlap 20 -T 8 -t "CDS" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout8"_CDS_raw.count"
./featureCounts -p -O --minOverlap 20 -T 8 -t "exon" -F GTF -a $myDictPath/$dict".gtf" $myDataPath/$IDout"_Aligned.out".sam \
-o $myDataPath/$IDout8"_exon_raw.count"


exit 0
                                    
