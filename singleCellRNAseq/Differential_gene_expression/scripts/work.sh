cp ../../../FELINE/FELINE_merged_011_028_premrna/Seurat/scRNA_Seurat_05individual_from_individual* ./
cp ../../../FELINE/FEL040P701/Seurat/scRNA_Seurat_00origial_patient_obj.* ./
ln -s ../Expression_Tables/VG_CAMA1_D11_10x_c* ./
ln -s ../Expression_Tables/VG_CAMA1_D11b_10x_c* ./
ln -s ../Expression_Tables/VG_CAMA1_D11_ALL_10x_c* ./

echo "VG_CAMA1_D11" > patients.list
echo "VG_CAMA1_D11b" > patients.list
echo "VG_CAMA1_D11_ALL" > patients.list

echo "preliminary filter and create raw obj"
sbatch --array 1 scRNA_Seurat_00origial_patient_obj.sh

echo "remove doublet"
ln -s ../QC/VG_CAMA1_D11_ALL_10x_Seurat_2kgenes_vst_cc.doublet_free.* ./

echo "cluster"
sbatch --array 1 scRNA_Seurat_05individual_from_individual.sh

echo "marker genes labeled"
ln -s ../Housekeeping_genes/VG_CAMA1_D11_10x_cell_metadata.UMAPcluster.marker_genes.txt ./
ln -s ../Housekeeping_genes/VG_CAMA1_D11b_10x_cell_metadata.UMAPcluster.marker_genes.txt ./
ln -s ../Housekeeping_genes/VG_CAMA1_D11_ALL_filtered/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt ./
sbatch --array 1 scRNA_Seurat_05individual_from_individual.sh

echo "compare DE genes in different comparsion"
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_top10_genes_Marker_genes.gene_diff.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_top10_genes_GV014mV_vs_GV014mC.gene_diff.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_top10_genes_GV014mV_vs_GV013mV.gene_diff.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_top10_genes_GV014mC_vs_GV015mC.gene_diff.txt

echo "compare DE genes between scRNA vs. BulkRNA: only 26 shared DE genes?"
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt ../../Vince_resistence_BulkRNA/Plot_DE_heatmap/CAMA_1_DMSO.vs.CAMA_1_riboR_DMSO.fdr_0.05.fc_2.txt

echo "rerun MAST using cells with 1.5k-2.5k genes, similar with all cells for DE (but a lot better for ssGSEA)"
sbatch --array 1 scRNA_Seurat_05individual_from_individual.sh
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV014mC_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_seurat_filtered_marker_genes_MAST_all_cells/VG_CAMA1_D11_10x_top10_genes_GV014mC_vs_GV015mC.gene_diff.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV014mV_vs_GV013mV.gene_diff.txt VG_CAMA1_D11_10x_seurat_filtered_marker_genes_MAST_all_cells/VG_CAMA1_D11_10x_top10_genes_GV014mV_vs_GV013mV.gene_diff.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt VG_CAMA1_D11_10x_seurat_filtered_marker_genes_MAST_all_cells/VG_CAMA1_D11_10x_top10_genes_GV013mV_vs_GV015mC.gene_diff.txt

echo "Merge two run "
cp ~/software/Project_scripts/scRNA/Seurat/Merged_patients_v1/scRNA_Seurat_01cluster_merge_object.R ./

echo "upregulated gene in GV014mC and GV014mV"
awk '$3>0' VG_CAMA1_D11b_10x_top10_genes_GV014mC_vs_GV013mV.gene_diff.txt > VG_CAMA1_D11b_10x_top10_genes_GV014mC_vs_GV013mV.gene_diff.GV014mC_up.txt
awk '$3>0' VG_CAMA1_D11b_10x_top10_genes_GV014mV_vs_GV013mV.gene_diff.txt > VG_CAMA1_D11b_10x_top10_genes_GV014mV_vs_GV013mV.gene_diff.GV014mV_up.txt

echo "pathway enrichment from DE genes"
sbatch scRNA_Seurat_05individual_from_individual_DE_genes.sh

