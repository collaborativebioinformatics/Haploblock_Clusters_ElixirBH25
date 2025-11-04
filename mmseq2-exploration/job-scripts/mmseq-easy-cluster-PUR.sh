#!/bin/sh
#PBS -N mmseq_easyClust_PUR
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=16gb

set -e
set -o pipefail

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate biohack2025_haplo

echo "start at $(date)"

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/Haploblock_Clusters_ElixirBH25/mmseq2-exploration

export POP_ID=PUR
export INPUT_FASTA=$BASE_DIR/data/TNF/data/${POP_ID}/haploblock_phased_seq_TNFa/haploblock_phased_seq_merged/chr6_region_31480875-31598421.fa
export RESULTS=TNF_${POP_ID}_chr6_region_31480875-31598421

export TEMP_FOLDER=$EPHEMERAL/BioHack2025_tempfiles/$RESULTS
export MIN_SEQ_ID=0.9936 #derived from jedrzej script
export COV_FRACTION=0.9942 # 1 - 682/HAPLO_LEN - 0.994198016096 rounded off to 0.9942
export COV_MODE=0
export PARAM_SET=2

echo "min_seq_id ${MIN_SEQ_ID}"
echo "C ${COV_FRACTION}"
echo "cov-mode ${COV_MODE}"

mkdir -p $TEMP_FOLDER

mkdir -p $BASE_DIR/results/params${PARAM_SET}/$RESULTS
cd $BASE_DIR/results/params${PARAM_SET}/$RESULTS

mmseqs easy-cluster $INPUT_FASTA $RESULTS $TEMP_FOLDER \
            --min-seq-id $MIN_SEQ_ID -c $COV_FRACTION --cov-mode $COV_MODE

echo "end at $(date)"
echo "DONE!"