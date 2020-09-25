#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=3:00:00
#SBATCH --output=run_slurm.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export PATH=$PATH:/home/jichen/software/BETSY/install/envs/scRNA/bin/:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
#export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library:/home/jichen/software/BETSY/install/envs/scRNAcnv/lib/R/library/

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

#cut into one sequence one file
#perl ~/software/bin/fastaDeal.pl --cuts 1 ~/Projects/Database/Reference/hg19/ucsc.hg19.fasta > hg19.chr1.fasta
#cp ucsc.hg19.fasta.cut/ucsc.hg19.fasta.01 hg19.chrM.fasta
#GENOME=/home/jichen/Projects/Breast/WGS/MDS_Kahn/Reference/hg19
GENOME=/home/jichen/Projects/Database/Reference
Database=/home/jichen/Projects/Database

#prefix=VG_CAMA1_D11.expression_counts
#prefix=VG_CAMA1_D11_10x_gene_symbols.filtered.counts
prefix=VG_CAMA1_D11_ALL_10x_gene_symbols.filtered.counts
prefix2=VG_CAMA1_D11_ALL
sample=`cat samples.list | head -n $N | tail -n 1`
if true; then
python Extract_sample_count.py --count $prefix\.txt --sample VG_CAMA1_D11.$sample\.sample.txt --project $sample
python Housekeeping_genes_expression.py --input $prefix\.$sample\.txt --cutoff 1 > $prefix\.$sample\.hkgene_cutoff1.txt
echo "cutoff=1"
grep "mC" $prefix\.$sample\.hkgene_cutoff1.txt
grep "mV" $prefix\.$sample\.hkgene_cutoff1.txt
python Housekeeping_genes_expression.py --input $prefix\.$sample\.txt --cutoff 2 > $prefix\.$sample\.hkgene_cutoff2.txt
echo "cutoff=2"
grep "mC" $prefix\.$sample\.hkgene_cutoff2.txt
grep "mV" $prefix\.$sample\.hkgene_cutoff2.txt
Rscript VG_CAMA1_D11_count01.counts.marker_cells.R $prefix $sample
fi
#Rscript VG_CAMA1_D11_count01.counts.marker_cells_violin.R $prefix
Rscript VG_CAMA1_D11_ALL_count01.counts.marker_cells_violin.R $prefix $prefix2

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

