#!/usr/bin/bash

# generate sequences 
Rscript 01_Generate_seq.R

# predict forward and reverse sequences
source /home/imm/hert6114/anaconda3/bin/activate myenv
python predict_seq.py
conda deactivate

# analyse predictions
Rscript 02_analyse_predictions.R