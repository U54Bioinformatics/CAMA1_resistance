
echo "zinbwave normalized count"
ln -s ../Normalize_count/VG_CAMA1_D11_ALL_10x.zinbwave.normalized.txt ./

echo "run ssGSEA with BETSY using zinbwave normalized count"
sbatch --array 1 betsy_scRNA_ssGSEA.sh

echo "rename files"
cp VG_CAMA1_D11_ALL_10x.ssGSEA/scores.txt VG_CAMA1_D11_ALL_10x.zinbwave.normalized.ssGSEA.scores.txt
cp VG_CAMA1_D11_ALL_10x.ssGSEA/expressed_genes.txt VG_CAMA1_D11_ALL_10x.zinbwave.normalized.ssGSEA.expressed_genes.txt
# Folder VG_CAMA1_D11_ALL_10x.all_filtered_cells has results for all filtered cells after removing low-quality cells and doublets (n=23232)
# Folder VG_CAMA1_D11_ALL_10x.all_filtered_with_marker_cells has results for filtered cells that express mVenus or mCherry. So the resistent and sensitive cells can be seperated. This is the final data Vince and Rena are using (n=16116). 
