#!/bin/bash
#$ -N brbseq_pipeline
#$ -cwd
#$ -o ../logs/brbseq_pipeline.out
#$ -e ../logs/brbseq_pipeline.err
#$ -l h_rt=24:00:00
#$ -l h_vmem=32G
#$ -pe smp 8


# Activate Conda env
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate brb_seq

# Run pipeline
bash brb_seq_pipeline.sh