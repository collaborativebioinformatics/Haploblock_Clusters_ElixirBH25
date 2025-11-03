#!/bin/sh
#PBS -N mmseq_easyClust_PUR
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=16gb

set -e
set -o pipefail

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate biohack2025_haplo

echo "start at $(date)"

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/mmseq2-exploration/

export INPUT_FASTA=$BASE_DIR/data/TNF/data/PUR/haploblock_phased_seq_TNFa/haploblock_phased_seq_merged/chr6_region_31480875-31598421.fa
export RESULTS=TNF_PUR_chr6_region_31480875-31598421
export TEMP_FOLDER=$EPHEMERAL/BioHack2025_tempfiles
export MIN_SEQ_ID=0.99
export COV_FRACTION=0.8
export COV_MODE=0

echo "min_seq_id ${MIN_SEQ_ID}"
echo "C ${COV_FRACTION}"
echo "cov-mode ${COV_MODE}"

mkdir -p $BASE_DIR/results/$RESULTS
cd $BASE_DIR/results/$RESULTS

mmseqs easy-cluster $INPUT_FASTA $RESULTS $TEMP_FOLDER \
            --min-seq-id $MIN_SEQ_ID -c $COV_FRACTION --cov-mode $COV_MODE

echo "end at $(date)"
echo "DONE!"