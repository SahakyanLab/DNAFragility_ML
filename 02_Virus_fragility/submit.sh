#!/usr/bin/bash

# process sequences and extract features
Rscript 00_Process_sequences.R

# get file names
cd ../data/
find human_viruses/ -type f -name "*final_pred*" -print0 | xargs -0 -n 1 -P 8 sh -c 'echo "$0" >> filenames.csv'
mv filenames.csv ../02_Virus_fragility/data/
cd ../02_Virus_fragility/

# run predictions
bash predict_viruses.sh

# analyse results
Rscript 02_analyse_predictions.R