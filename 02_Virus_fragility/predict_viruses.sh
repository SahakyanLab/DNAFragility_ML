#!/usr/bin/bash

source activate fragility_model

for rand in "" "--random"
do
    python predict_viruses.py $rand
done

conda deactivate