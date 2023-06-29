# RNAseq

![alt text](RNApipeline.png)


... What we have to do


-**Preparing a fasta reference**
  - First we need a DNA reference (fasta) of the organism that we are investigating. 
    - we can optionally adjust the fasta to our lab strain (change SNPs, Indels)
    - Create an index with the script: **0_DicIndex.sh**
    - We need a list of features, (genes in this case), for which we count the reads that map to it. In this case we have it in the gtf format
