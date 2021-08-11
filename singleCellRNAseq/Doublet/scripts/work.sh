
echo "Use seurat to do an initial filter (keep almost all cells) in folder Seurat"
ln -s ../Seurat/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.raw.rds ./

echo "run scrublet to label doublet and prepare files for Seurat anlaysis"
# input: VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.raw.rds
# output 1: VG_CAMA1_D11_ALL_10x.scrublet_out.csv, doublet prediction: use "scrublet_doublet_call2" which use 0.25 as cutoff for scrublet 
# output 2: Seurat_2kgenes_vst_cc.doublet_free.cell_metadata.txt, meta file of doublet free cells
# output 3: Seurat_2kgenes_vst_cc.doublet_free.rds, seurat obj of doublet free calls 
sbatch scrublet_QC.sh

