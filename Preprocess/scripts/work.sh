
echo "sample file"
ln -s ~/Projects/Breast/scRNA/Data/Megatron_data/COH060/FT-SE5541/ ./VG_CAMA1_D11_fastq
sbatch betsy_sample.sh
#deeper sequence
ln -s /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH061/FT-SE5576/ VG_CAMA1_D11b_fastq
sbatch betsy_sample.sh

echo "Merge fastq dirs from two runs into one folder"
mkdir VG_CAMA1_D11_ALL_fastq
cd VG_CAMA1_D11_ALL_fastq/
ln -s ../VG_CAMA1_D11b_fastq/FT-SA4266* ./
ln -s ../VG_CAMA1_D11_fastq/FT-S* ./
cp ../VG_CAMA1_D11_fastq/barcodeAssociationTable.txt ./
cat ../VG_CAMA1_D11b_fastq/barcodeAssociationTable.txt >> barcodeAssociationTable.txt 
cp VG_CAMA1_D11b_sample.txt VG_CAMA1_D11_ALL_sample.txt
cat VG_CAMA1_D11_sample.txt >> VG_CAMA1_D11_ALL_sample.txt


echo "preprocess"
sbatch betsy_scRNA_10x_QC.sh

