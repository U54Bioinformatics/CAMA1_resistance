#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --output=One_way_ANOVA_normalcell_slurm.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

start=`date +%s`

export R_LIBS=~/software/BETSY/install/lib/R/library/

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


patient=`cat patients.list | head -n $N | tail -n 1`
score=$patient\_10x.ssGSEA.scores.txt
cluster=$patient\_10x.seurat.cell_type.anno.txt
#1:#comparing GV014_mvenus_vs_GV013_mvenus
#2:#comparing GV013_mvenus_vs_GV015_mcherry
#3:#mvenus_vs_mcherry
#4:#comparing GV014_mcherry_vs_GV015_mcherry
#5:#comparing GV014_mcherry_vs_GV014_mvenus
analysis=4

echo "One way ANOVA"
python One_way_ANOVA_normalcell_CAMA1.py --input $score --cluster $cluster --output $patient\_ANOVA_dir --analysis $analysis
perl multi-process.pl -cpu $CPU $patient\_ANOVA_dir_$analysis\.run.sh
python One_way_ANOVA_normalcell_CAMA1.py --input $score --cluster $cluster --output $patient\_ANOVA_dir --analysis $analysis --summary 1
cat $patient\_ANOVA_dir_$analysis\.test.AdjustP.R | R --slave


echo "Violin plots"
python One_way_ANOVA_cluster_top_plot.py --input $patient\_ANOVA_dir_$analysis\.test.AdjustP.txt --fdr 0.05 --diff 0.1 --top 50
python One_way_ANOVA_cluster_top_plot.py --input $patient\_ANOVA_dir_$analysis\.test.AdjustP.txt --fdr 0.05 --diff 0.05 --top 50
python One_way_ANOVA_cluster_top_plot.py --input $patient\_ANOVA_dir_$analysis\.test.AdjustP.txt --fdr 0.05 --diff 0 --top 50
#mv $patient\_* FEL011016_ssGSEA_ANOVA_seurat_cluster/

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

