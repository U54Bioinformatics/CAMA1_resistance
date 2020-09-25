#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=scrublet_QC.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh



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

export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/lib/R/library
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_QC/lib/python3.6/site-packages/

#sample=FEL025P142
#sample=FEL012P101
platform=10x
sample=`cat patients.list | head -n $N | tail -n 1`

echo "extract market matrix from seurat obj"
/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/bin/Rscript scrublet_extract_mm.R $sample GV013 $platform
/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/bin/Rscript scrublet_extract_mm.R $sample GV014 $platform
/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/bin/Rscript scrublet_extract_mm.R $sample GV015 $platform

echo "predict doublet with scrublet"
/home/jichen/software/BETSY/install/envs/scRNA_QC/bin/python scrublet_doublet.py $sample\_GV013
/home/jichen/software/BETSY/install/envs/scRNA_QC/bin/python scrublet_doublet.py $sample\_GV014
/home/jichen/software/BETSY/install/envs/scRNA_QC/bin/python scrublet_doublet.py $sample\_GV015

echo "merge doublet annotation csv files"
python scrublet_merge_meta.py --list patients.list --output $sample

echo "write doublet free RDS and meta"
/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/bin/Rscript scrublet_save_rds.R $sample\_$platform

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

