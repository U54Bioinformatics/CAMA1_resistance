#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=scRNA_Seurat_05individual_from_individual_DE_genes.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_enricher/lib/R/library


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

#comparison1
#prefix=VG_CAMA1_D11_ALL_10x_top10_genes_GV014mV_vs_GV013mV
#comparison=GV014mVenus_vs_GV013mVenus
#foldchange=0.2
#comparison2
#prefix=VG_CAMA1_D11b_10x_top10_genes_GV014mC_vs_GV015mC
#comparison=GV014mCherry_vs_GV015mCherry
#foldchange=0.2
#comparison3
#prefix=VG_CAMA1_D11b_10x_top10_genes_GV014mC_vs_GV013mV
#comparison=GV014mCherry_vs_GV013mVenus
#foldchange=0.2
#comparison4
#prefix=VG_CAMA1_D11_ALL_10x_top10_genes_GV013mV_vs_GV015mC
#comparison=GV013mVenus_vs_GV015mCherry
#flip the comarision for NES figure
#comparison=GV015mCherry_vs_GV013mVenus
#foldchange=0.2
#comparison5
prefix=VG_CAMA1_D11_ALL_10x_top10_genes_GV014mV_vs_GV014mC
#comparison=GV014mVenus_vs_GV014mCherry
#flip the comarision for NES figure
comparison=GV014mCherry_vs_GV014mVenus
foldchange=0.2

echo "Filter significant genes and Vocanol plot"
/home/jichen/software/BETSY/install/envs/scRNA_enricher/bin/Rscript scRNA_Seurat_05individual_from_individual_DE_genes_plot.R $prefix $comparison $foldchange
convert -density 300 -quality 100 $prefix.VolcanoPlot.pdf $prefix.VolcanoPlot.png
python ~/software/bin/txt2xlsx.py --input $prefix\.gene_diff.txt
python ~/software/bin/txt2xlsx.py --input $prefix\.gene_diff.significant.txt
echo "Enrichment test: hypergeometric test and GSEA"
/home/jichen/software/BETSY/install/envs/scRNA_enricher/bin/Rscript scRNA_Seurat_05individual_from_individual_DE_genes_enricher.R $prefix $comparison $foldchange
convert -density 300 -quality 100 $prefix\.Hallmark_NES.pdf $prefix\.Hallmark_NES.png
python ~/software/bin/txt2xlsx.py --input $prefix\.c2.GSEA.txt
python ~/software/bin/txt2xlsx.py --input $prefix\.hallmark.GSEA.txt
python ~/software/bin/txt2xlsx.py --input $prefix\.c7.GSEA.txt

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

