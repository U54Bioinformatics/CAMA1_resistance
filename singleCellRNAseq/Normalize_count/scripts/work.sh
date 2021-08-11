echo "Prepare count"
# all filtered cells
ln -s ../Seurat/VG_CAMA1_D11_ALL_10x_gene_symbols.raw.counts.txt ./
ln -s ../Seurat/VG_CAMA1_D11_ALL_10x_gene_symbols.CPM.txt ./
cp ../Seurat/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt .
ln -s VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt VG_CAMA1_D11_ALL_10x.metadata.txt
# cell with marker genes
ln -s ../Plot_gene_violin_seurat/VG_CAMA1_D11_ALL_10x.expression_counts.raw.txt ./
ln -s ../Plot_gene_violin_seurat/VG_CAMA1_D11_ALL_10x.extract.metadata.txt ./
ln -s VG_CAMA1_D11_ALL_10x.expression_counts.raw.txt VG_CAMA1_D11_ALL_10x.txt
ln -s VG_CAMA1_D11_ALL_10x.extract.metadata.txt VG_CAMA1_D11_ALL_10x.metadata.txt

echo "Run zinbwave normalization"
sbatch Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-12 Normalize_count_by_zinbwave.02_normalization.sh
sbatch Normalize_count_by_zinbwave.03_merge_output.sh
