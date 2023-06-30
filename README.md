# RNAseq

![alt text](RNApipeline.png)


... What we have to do


-**Preparing a fasta reference**
  - First we need a DNA reference (fasta) of the organism that we are investigating. 
    - we can optionally adjust the fasta to our lab strain (change SNPs, Indels)
    - Create an index with the script: [0_DicIndex.sh](https://github.com/Easybel/RNAseq/blob/master/0_DicIndex.sh)
    - We need a list of features, (genes in this case), for which we count the reads that map to it. In this case we have it in the gtf format

- **Prepare the raw data**
  - **fastq**: Raw data is provided in the .fastq format, where we have paired-end reads (forward _1 & reverse _2), compressed in the .fastq.gz format.  
  - First, run a Quality control and trimm the data with: [1_QCScript.sh](https://github.com/Easybel/RNAseq/blob/master/1_QCScript.sh)
    - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) checks the quality of reads before and after trimming. A while ago, we had the problem that the rRNA depletion was not working int he company. This step can help to identify this problem.
    - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) trimms bad quality reads and removes adapters, that are specified by the company
       

- **Mapping the reads**
  - With the script [2_CombiScriptStar.sh](https://github.com/Easybel/RNAseq/blob/master/2_CombiScriptStar.sh) The trimmed reads are mapped with [Star](https://github.com/alexdobin/STAR) against the fasta reference
    - keep in mind that we have to decide what to do with multimapping reads!    
  - with featurecounts (subread), we count the number of reads that map to a specific feature, specified for the reference in the gtf file
  
- **Count quantification/ analysis**
  1. we perform a pre-analysis with [4a_PreAnalysis.R](https://github.com/Easybel/RNAseq/blob/master/4a_PreAnalysis.R) to find out, if the data clusters nicely in a biological way (as opposed to clustering by day). The principle component analysis that was used for this analysis was also implemented by Isabel in the python script in [4a_PreAnalysis_PCA_IR.ipynb](https://github.com/Easybel/RNAseq/blob/master/4a_PreAnalysis_PCA_IR.ipynb). Here, it is explained step by step what has to be done. Further explanations can be found in her phD thesis. 
  
  
- **Other software**
  - **IGV**: http://software.broadinstitute.org/software/igv/

