# RNAseq

![alt text](RNApipeline.png)


... What we have to do


-**Preparing a fasta reference**
  - First we need a DNA reference (fasta) of the organism that we are investigating. 
    - we can optionally adjust the fasta to our lab strain (change SNPs, Indels)
    - Create an index with the script: **0_DicIndex.sh**
    - We need a list of features, (genes in this case), for which we count the reads that map to it

- **Prepare the raw data**
  - **fastq**: Raw data is provided in the .fastq format, where we have paired-end reads (forward _1 & reverse _2), compressed in the .fastq.gz format.  
  - First, run a Quality control and trimm the data with: **1_ QCScript.sh**
    - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) checks the quality of reads before and after trimming. A while ago, we had the problem that the rRNA depletion was not working int he company. This step can help to identify this problem.
    - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) trimms bad quality reads and removes adapters, that are specified by the company
       

- **Reads mappen**
  - Die getrimmten reads werden mit [Star](https://github.com/alexdobin/STAR) auf die Referenz gemappt
  - mit featurecounts (subread) wird die Anzahl der reads bestimmt, die auf die feature entfallen
  
- **Count quantification/ analysis**
  - Was bedeuten die counts pro feature jetzt? --> das wird in R (matlab?) analysiert 
  
  
- **Other software**
  - **IGV**: http://software.broadinstitute.org/software/igv/

