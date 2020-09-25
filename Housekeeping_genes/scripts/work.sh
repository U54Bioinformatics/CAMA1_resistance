
echo "Identify cell expressing marker genes"
# add Gene.ID to the first if there is no one
cp ../Seurat/results/VG_CAMA1_D11_ALL_10x_gene_symbols.filtered.counts.txt ./
ln -s ../Seurat/results/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.txt ./
sbatch --array 1-3 run_slurm.sh

echo "all cell with markers"
# output: VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt, cell meta with marker gene labeled (Marker_genes    Marker_groups); Filter out "nomarker" in column "Marker_genes" to obtain the final cells (n=16116)
bash run_marker_anno.sh


