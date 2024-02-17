#!/bin/bash

pwd="$(pwd)/"
RNAfold_path=$1

# extract SV and SNP sequences
Rscript 00_Process_sequences.R $pwd

# extract features
Rscript 00_Process_SNP.R $pwd $RNAfold_path

for batch in {1..40}
do
   Rscript 01_Get_Features_SV.R $pwd $sub_batch $RNAfold_path
done

# run predictions
bash predict_SVs.sh

# get file names
cd ../data/
find deltafragility/ -type f -name "*final_pred*" -print0 | xargs -0 -n 1 -P 8 sh -c 'echo "$0" >> filenames.csv'
mv filenames.csv ../05_DeltaFragility/data/
cd ../05_DeltaFragility/
Rscript 02_Process_SV_predictions.R $pwd

# analyse results
Rscript 02_analyse_predictions_SNPs.R $pwd

for bw in 1960 10000 20000
do
   Rscript 03_analyse_predictions_SVs.R $pwd $bw
done