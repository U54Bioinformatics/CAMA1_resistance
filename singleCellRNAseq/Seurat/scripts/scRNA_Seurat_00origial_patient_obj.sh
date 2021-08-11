#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=3:00:00
#SBATCH --output=scRNA_Seurat_00origial_patient_obj.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

###########################################
#input file:
#10x count file:      FEL011_10x_count01.noaggr.counts.txt
#Betsy cell metafile: FEL011_10x_cell_metadata.txt
     
#output file:
#seurat project file: FEL011_10x_Seurat_2kgenes_vst_cc.rds
###########################################

#sample=FEL012016
#sample=FEL023P121
platform=10x
sample=`cat patients.list | head -n $N | tail -n 1`

if [ ! -e $patient\_Seurat_2kgenes_vst_cc.raw.rds ]; then
   echo "Create seurat obj for patient: $patient"
   if [ $platform == '10x' ]; then
      /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_00origial_patient_obj.R $sample $platform
   fi
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

