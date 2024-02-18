#!/usr/bin/bash

pwd="$(pwd)/"
RNAfold_path=$1

# generate sequences 
Rscript 01_Generate_seq.R $pwd $RNAfold_path

# predict forward and reverse sequences
source activate fragility_model
python predict_seq.py
conda deactivate

# analyse predictions
Rscript 02_analyse_predictions.R $pwd