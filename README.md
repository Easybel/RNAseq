# RNAseq

![alt text](RNApipeline.png)

## TheMaster Plan
... What we have to do

Here I collected all the things that I could think of that we should work on until Viera leaves - and of course also afterwards :)
  
- **Rohdaten bearbeiten**
  - **fastq**: Daten sind im .fastq Format gegeben, da wir paired-end data haben immer mit forward _ 1 und reverse _ 2 run. Die Daten sind komprimiert und liegen als .fastq.gz vor.  
  - Quality control und trimmen: **1_ QCScript.sh**
    - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) benutzen, um die Qualität der Rohdaten zu checken.
    - [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) benutzen, um Adapter und zu kurze/ zu schlechte reads zu trimmen.
    - (check die Qualität nochmal - verbessert?)
 
-**Fasta Referenz erzeugen**
  - ...
  - Indexing

- **Reads mappen**
  - Die getrimmten reads werden mit [Star](https://github.com/alexdobin/STAR) auf die Referenz gemappt. 
  
  
- **Other software**
  - **IGV**: http://software.broadinstitute.org/software/igv/

