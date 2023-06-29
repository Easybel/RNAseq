# This script takes the path where the FastQC results lie and returns 
# the overrepresented sequences from the report in a fasta format

myDataIn="/projects/ag-advol/RNAseq/Raw/QC"
myDataOut="/projects/ag-advol/RNAseq/Raw/Out"
myScriptPath="/home/irathman/scripts/python_scripts"

module unload python/2.7.5
module load python/3.4.3

# inside the unzipped folder the data is usually called
report="fastqc_data.txt"

cd $myDataIn

for f in Ngo*1_fastqc.zip; do

ID=$(echo $f | cut -d"." -f1)
IDout=$(echo $f | cut -d"_" -f1,2,3)

echo $ID
echo $IDout

rm -r $myDataIn/$ID
unzip $myDataIn/$ID".zip"
echo "done unzipping"

# here the sequences are 
cd $myDataIn/$ID
sed -n '/Overrepresented/,/END_MODULE/p' $report | grep -v "^#" | grep -v "^>>" | awk '{print $1}' > $myDataOut/$IDout"_OverSeq.txt"

cd $myDataOut
IDout2=$IDout"_OverSeq"

python3 $myScriptPath/Seqs2Fasta.py "$IDout2"

rm -r $myDataIn/$ID

## part where blast is used

module unload intel/19.0
module load  intel/16.0_gnu_5.1.0

blastn –db nt –query $myDataIn/Ngo_T1ctl_S2_OverSeq.fasta –out $myDataOut/$IDout2"_results.out" -remote

done
