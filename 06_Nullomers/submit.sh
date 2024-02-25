#!/usr/bin/bash

pwd="$(pwd)/"
repeats=10000

# generate sequences 
Rscript 00_Nullomer_and_genome.R $pwd $repeats
Rscript 01_Nullomer_and_shuffle.R $pwd $repeats
Rscript 02_SNP_nullomer.R $pwd $repeats
# Rscript 02_SNP_nullomer_v2.R $pwd $repeats

source /home/imm/hert6114/anaconda3/bin/activate myenv

# predict fragility of SNPs leading to nullomers
python predict_seq.py --SV_type SNP

# predict random sequences composed of nullomers vs. reference genome
python predict_seq.py --SV_type Genome

# predict random sequences composed of nullomers vs. shuffled version
python predict_seq.py --SV_type Shuffle

conda deactivate

# analyse predictions
Rscript 03_analyse_shuffle_predictions.R $pwd
Rscript 04_analyse_SNP_nullomer.R $pwd

get file names
cd ../data/
find nullomer_and_genome/ -type f -name "*final_pred*" -print0 | xargs -0 -n 1 -P 8 sh -c 'echo "$0" >> nullomer_and_genome.csv'
mv nullomer_and_genome.csv ../06_Nullomers/data/
cd ../06_Nullomers/

# get file names
cd ../data/
find nullomer_and_shuffle/ -type f -name "*final_pred*" -print0 | xargs -0 -n 1 -P 8 sh -c 'echo "$0" >> nullomer_and_shuffle.csv'
mv nullomer_and_shuffle.csv ../06_Nullomers/data/
cd ../06_Nullomers/

for shuffle in "TRUE" "FALSE"
do
    Rscript 03_analyse_genome_predictions.R $pwd $shuffle
done