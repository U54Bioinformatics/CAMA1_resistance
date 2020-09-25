echo "VG_CAMA1_D11_ALL, combined fastq sequence from two runs"
ln -s ../Preprocess/results/VG_CAMA1_D11_ALL_count01.noaggr.txt VG_CAMA1_D11_ALL_10x_count01.noaggr.counts.txt
ln -s ../Preprocess/results/VG_CAMA1_D11_ALL_meta/cell_metadata.txt VG_CAMA1_D11_ALL_10x_cell_metadata.txt

echo "preliminary filter and create raw obj: min_cell=3, filter gene expressed in less than 3 cells;"
sbatch --array 1 scRNA_Seurat_00origial_patient_obj.sh

echo "do doublet prediction in ../Doublet and link the results back here"
ln -s ../Doublet/results/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.doublet_free.* ./

echo "run Seurat and plot"
# input are seurat obj and meta, raw and doublet_free
# output 1: VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.rds, seurat obj on filtered cell by removing low-quality cells and doublet 
# output 2: VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.pdf, seurat plot of above seurat obj
# figure folder:
# VG_CAMA1_D11_ALL_10x_seurat_raw_png, plot of raw cells, first step filter, removing gene expressed in less than 3 cells
# VG_CAMA1_D11_ALL_10x_seurat_doubletfree_png, plot of doublet free cells, second step filter removing doublet
# VG_CAMA1_D11_ALL_10x_seurat_filtered_png, plot of filtered cell, third step filter removing anything outside cutoff range
# VG_CAMA1_D11_ALL_10x_seurat_figures_png, similar to VG_CAMA1_D11_ALL_10x_seurat_filtered_png, ploting on same cells 
sbatch --array 1 scRNA_Seurat_05individual_from_individual.sh

echo "removing cells without confident marker genes and plot"
cp ../Housekeeping_genes/output/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt ./
sbatch --array 1 scRNA_Seurat_05individual_from_individual_marker_genes.sh

