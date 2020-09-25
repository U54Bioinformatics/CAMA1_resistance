cp ~/software/Project_scripts/scRNA/ssGSEA_ANOVA/One_way_ANOVA_normalcell_slurm.sh ./
cp ~/software/Project_scripts/scRNA/ssGSEA_ANOVA/One_way_ANOVA_normalcell.py ./
cp ~/software/Project_scripts/scRNA/ssGSEA_ANOVA/multi-process.pl ./
ln -s ../ssGSEA/VG_CAMA1_D11_10x.ssGSEA/scores.txt ./VG_CAMA1_D11_10x.ssGSEA.scores.txt
awk '{print $1"\t"$5"\t"$21"\t"$22"\t"$23}' ../Seurat/VG_CAMA1_D11_10x_cell_metadata.UMAPcluster.marker_genes.txt > VG_CAMA1_D11_10x.seurat.cell_type.anno.txt

echo "anova pathway analysis"
sbatch --array 1 One_way_ANOVA_normalcell_slurm.sh

echo "anova pathway analysis with nfeature between 1.5 to 2.5 k"
head -n 1 VG_CAMA1_D11_10x.seurat.cell_type.anno.txt > VG_CAMA1_D11_10x.seurat.cell_type.anno.1.5_2.5k.txt
awk '$4>1500 && $4<2500' ../Seurat/VG_CAMA1_D11_10x_cell_metadata.UMAPcluster.marker_genes.txt | awk '{print $1"\t"$5"\t"$21"\t"$22"\t"$23}' >> VG_CAMA1_D11_10x.seurat.cell_type.anno.1.5_2.5k.txt
#run 
sbatch --array 1 One_way_ANOVA_normalcell_slurm.sh
#compare 1.5 to 2.5 k cells vs. all cells
python ~/software/bin/listdiff.py VG_CAMA1_D11_ANOVA_dir_1.test.AdjustP.fdr0.05_diff0.05.txt VG_CAMA1_D11_ANOVA.all_cell_with_markers/VG_CAMA1_D11_ANOVA_dir_1.test.AdjustP.fdr0.05_diff0.05.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_ANOVA_dir_2.test.AdjustP.fdr0.05_diff0.05.txt VG_CAMA1_D11_ANOVA.all_cell_with_markers/VG_CAMA1_D11_ANOVA_dir_2.test.AdjustP.fdr0.05_diff0.05.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_ANOVA_dir_4.test.AdjustP.fdr0.05_diff0.05.txt VG_CAMA1_D11_ANOVA.all_cell_with_markers/VG_CAMA1_D11_ANOVA_dir_4.test.AdjustP.fdr0.05_diff0.05.txt
python ~/software/bin/listdiff.py VG_CAMA1_D11_ANOVA_dir_5.test.AdjustP.fdr0.05_diff0.05.txt VG_CAMA1_D11_ANOVA.all_cell_with_markers/VG_CAMA1_D11_ANOVA_dir_5.test.AdjustP.fdr0.05_diff0.05.txt
#count cells
grep "GV013_mvenus" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.txt
grep "GV014_mvenus" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.txt
grep "GV014_mcherry" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.txt
grep "GV015_mcherry" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.txt
grep "GV013_mvenus" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.all_cell.txt
grep "GV014_mvenus" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.all_cell.txt
grep "GV014_mcherry" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.all_cell.txt
grep "GV015_mcherry" -c VG_CAMA1_D11_10x.seurat.cell_type.anno.all_cell.txt


echo "update deeper sequence"
ln -s ../ssGSEA/VG_CAMA1_D11b_10x.ssGSEA/scores.txt VG_CAMA1_D11b_10x.ssGSEA.scores.txt
ln -s ../Seurat/VG_CAMA1_D11b_10x_cell_metadata.UMAPcluster.marker_genes.txt ./

#extract pathway score for ploting or test
export R_LIBS=/home/jichen/software/BETSY/install/envs/limma/lib/R/library/
python Prepare_ssGSEA_data_per_pathway_CAMA1.py --input VG_CAMA1_D11b_10x.ssGSEA.scores.txt --anno VG_CAMA1_D11b_10x_cell_metadata.UMAPcluster.marker_genes.txt --pathway pathways.list > log 2>&1 &
cat DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP.data.violin.R | R --slave
convert -density 300 -quality 100 DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP.pdf DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP.png

echo "merged run1 and run2"
ln -s ../ssGSEA/VG_CAMA1_D11_ALL_10x.ssGSEA/scores.txt VG_CAMA1_D11_ALL_10x.ssGSEA.scores.txt
awk '{print $1"\t"$5"\t"$22"\t"$23"\t"$24}' ../Seurat/VG_CAMA1_D11_ALL_10x_cell_metadata.UMAPcluster.marker_genes.txt > VG_CAMA1_D11_ALL_10x.seurat.cell_type.anno.txt

