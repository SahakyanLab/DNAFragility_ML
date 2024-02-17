#!/usr/bin/bash

source /home/imm/hert6114/anaconda3/bin/activate myenv

# # predict randomly generated DNA virus sequences
# python predict_viruses.py --random True

# predict DNA virus sequences
python predict_viruses.py --random False

conda deactivate