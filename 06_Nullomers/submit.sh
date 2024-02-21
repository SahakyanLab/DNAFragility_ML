#!/usr/bin/bash

pwd="$(pwd)/"

# generate sequences 
Rscript 01_Nullomer_and_shuffle.R $pwd
Rscript 02_SNP_nullomer.R $pwd
# Rscript 02_SNP_nullomer_v2.R $pwd

source /home/imm/hert6114/anaconda3/bin/activate myenv

# predict fragility of SNPs leading to nullomers
python predict_seq.py --SV_type SNP

# predict random sequences composed of nullomers vs. reference genome
python predict_seq.py --SV_type Not_SNP

conda deactivate

# analyse predictions
Rscript 03_analyse_predictions.R $pwd
Rscript 04_analyse_SNP_nullomer.R $pwd