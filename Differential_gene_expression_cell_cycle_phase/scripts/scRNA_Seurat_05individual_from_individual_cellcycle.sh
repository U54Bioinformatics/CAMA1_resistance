#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=scRNA_Seurat_05individual_from_individual_cellcycle.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
#export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_MAST/lib/R/library

export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=3
fi

echo "CPU: $CPU"
echo "N: $N"

cc=G1
sample=VG_CAMA1_D11_ALL_$cc
platform=10x
#sample=`cat patients.list | head -n $N | tail -n 1`
prefix=$sample\_$platform

if [ ! -e $sample\_$platform\_seurat_raw_png ]; then
   echo "Plot raw, no filtered"
   #mkdir $sample\_$platform\_seurat_raw_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_raw.R $sample $platform
   #mkdir $sample\_$platform\_seurat_doubletfree_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_doubletfree.R $sample $platform
fi

if [ ! -e $sample\_$platform\_seurat_figures_png ]; then
   echo "Running seurat"
   #mkdir $sample\_$platform\_seurat_figures_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_before_singleR.R $prefix
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_singleR_anno.R $prefix
fi

#to rerun singleR
#/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_singleR_anno.R $prefix

#clusters need to be edited based on UMAP and cell type annotation, set to 1 and edit;
#bioinformatics.stackexchange.com/questions/4297/resolution-parameter-in-seurats-findclusters-function-for-larger-cell-numbers
if [ 0 ]; then
   echo "Editing clusters"
   #mkdir $sample\_$platform\_seurat_edit_cluster_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_edit_cluster.R $sample $platform
   #mv $sample\_$platform\_Seurat_2kgenes_vst_cc.rds $sample\_$platform\_Seurat_2kgenes_vst_cc.before_edit.rds
   #cp $sample\_$platform\_Seurat_2kgenes_vst_cc.edit_cluster.rds $sample\_$platform\_Seurat_2kgenes_vst_cc.rds
   #mv $sample\_$platform\_cell_metadata.UMAPcluster.txt $sample\_$platform\_cell_metadata.UMAPcluster.before_edit.txt
   #cp $sample\_$platform\_cell_metadata.UMAPcluster.edit_cluster.txt $sample\_$platform\_cell_metadata.UMAPcluster.txt
fi

if [ ! -e $sample\_$platform\_seurat_filtered_png ]; then
   echo "Plot filtered"
   #mkdir $sample\_$platform\_seurat_filtered_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_filtered.R $sample $platform
fi

if [ ! -e $sample\_$platform\_seurat_filtered_marker_genes_png ]; then
   echo "Plot filtered marker genes"
   mkdir $sample\_$platform\_seurat_filtered_marker_genes_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_filtered_test.R $sample $platform
   /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_filtered_with_markers_CC.R $sample $platform $cc
fi

if [ ! -e $sample\_$platform\_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt ]; then
   echo "Postprocess after singleR"
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_sample_cluster_matrix.R $prefix
   #python scRNA_Seurat_05individual_from_individual_sum_cell_type.py --cluster_anno $prefix\_SingleR_cluster_table_Blueprint.txt --sample_cluster $prefix\_cell_metadata.UMAPcluster.SingleR_anno_cluster.sample_matrix.txt
   #python scRNA_Seurat_05individual_from_individual_add_cell_type_anno.py --cluster_anno $prefix\_SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster $prefix\_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt 
fi


#if [ ! -e $sample\_$platform\_seurat_after_singler_png ]; then
   echo "Plot after singleR, filtered and cell types annotated"
   #mkdir $sample\_$platform\_seurat_after_singler_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "umap" 0
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "umap" 1
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "umap" 2
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "tsne" 0
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "tsne" 1
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform "tsne" 2
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png.R $sample $platform
#fi

if [ ! -e $sample\_$platform\_seurat_singler_infercnv_figures_png ]; then
   echo "Plot after singleR and infercnv"
   #mkdir $sample\_$platform\_$patient\_seurat_singler_infercnv_figures_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_after_infercnv.R $sample $platform $patient
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

