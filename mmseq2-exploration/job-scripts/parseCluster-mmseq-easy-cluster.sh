#!/bin/sh
#PBS -N mmseq_clusterParams
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=16gb

set -e
set -o pipefail

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate biohack2025_haplo

echo "start at $(date)"

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/Haploblock_Clusters_ElixirBH25/mmseq2-exploration
export PARAM_SET=0

export RESULTS_FOLDER=$BASE_DIR/results/params${PARAM_SET}/

read p1 p2 p3 < <(tail -n +2 $RESULTS_FOLDER/params${PARAM_SET}.tsv)

echo -e "MIN_SEQ_ID\tCOV_FRACTION\tCOV_MODE\tSUBFOLDER_NAME\tNUM_CLUSTER" > $RESULTS_FOLDER/params${PARAM_SET}_cluster.tsv

for d in $RESULTS_FOLDER/TNF_*; do
    n=$(grep -c "^>" $d/*_rep_seq.fasta)
    echo -e "${p1}\t${p2}\t${p3}\t${d}\t${n}" >> $RESULTS_FOLDER/params${PARAM_SET}_cluster.tsv
done

echo "end at $(date)"
echo "DONE!"