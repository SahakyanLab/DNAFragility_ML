#!/usr/bin/bash

source activate fragility_model

# # predict randomly generated DNA virus sequences
# python predict_viruses.py --random True

# predict DNA virus sequences
python predict_viruses.py --random False

conda deactivate