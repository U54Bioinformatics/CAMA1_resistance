#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=scRNA_Seurat_05individual_from_individual.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

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

#sample=FEL011016
platform=10x
sample=`cat patients.list | head -n $N | tail -n 1`
prefix=$sample\_$platform

if [ ! -e $sample\_$platform\_seurat_filtered_marker_genes_png ]; then
   echo "Plot filtered marker genes"
   mkdir $sample\_$platform\_seurat_filtered_marker_genes_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_filtered_test.R $sample $platform
   /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_featureplots_png_filtered_with_markers.R $sample $platform
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

