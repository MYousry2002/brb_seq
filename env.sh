#!/bin/bash

# Create environment
conda init
conda env create -f env.yml

cd . 
conda init
conda activate brb_seq

# Path to BRBseqTools JAR
export BRBSEQ_JAR_PATH="$(pwd)/BRBseqTools.jar"